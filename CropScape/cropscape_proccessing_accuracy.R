library(rgdal)
library(raster)
library(sf)
library(ncdf4)
library(dplyr)
library(lubridate)
library(imputeTS)
library(filesstrings)
library(data.table)
library(tools)
library(stringi)
library(raster)
library(parallel)
library(rlist)
library(exactextractr)
library(parallel)
library(data.table)
library(tictoc)
library(snow)
library(tidyr)
library(ggplot2)
library(extrafont)
library(randomcoloR)
options(scipen = 999)

###############################################################################
#++++++++++++++++++++++ Set the working directory  ++++++++++++++++++++++++++ #
##############################################################################
#RM temp files
tmp_dir <- tempdir()
list.files(tmp_dir)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), 
       recursive = TRUE)
source("D:/Projects/Project_Mother/Scripts/kings.R")
source("D:/Projects/Project_Mother/Scripts/map.R")

########## !!!!!!! You need to change the below for different extractions ######
# this is where the raw raster data is sitting
dir <- "[YOUR IMAGERY DIRECTORY HERE]"
setwd(dir)

## Read in the files if not uploading from R data (this will change w/ every extraction)
scenes = list.files(path=setwd(dir),
                    pattern=glob2rx("*cdls*.tif"), recursive = TRUE)

## Make sure to sub out the years you don't need
scenes.crs= scenes[1]

################################################################################
# +++++++++++++++++++ Read in Shapefiles and set up code++++++++++++++++++++++ #
################################################################################
HRU <- sf::st_read("[YOUR SHAPEFILE HERE]")
HRU = sf::st_transform(HRU, crs = proj4string(brick(scenes.crs)))

area = 900 ## 900 m2 to acres
store = exact_extract(raster::stack(scenes), HRU)
names = HRU$NAME
names(store) = names
map.dat = Map(cbind, store, ID = names)
df.dat = list.rbind(map.dat)
df.grouped = df.dat %>% reshape2::melt(id = c("coverage_fraction", "ID"))
df.grouped.sum = df.grouped %>% group_by(ID, variable, value) %>%
  dplyr::summarize(total = sum(coverage_fraction * area))

df.grouped.sum$Year = str_remove(df.grouped.sum$variable, "_30m_cdls_clip")
df.grouped.sum$Year = str_remove(df.grouped.sum$Year, "X")
df.grouped.sum$ID = as.character(as.factor(df.grouped.sum$ID))

################################################################################
# +++++++++++++++++++++++++++++++ Accuracy Assessment +++++++++++++++++++++++++#
################################################################################

WD_ST = read.csv("[LINK TO LOOKUP FILE BETWEEN STATE & UNIT]")
WD_ST$NAME = as.character(as.factor(WD_ST$NAME))
WD_ST$State = as.character(as.factor(WD_ST$State))

Accuracy = read.csv("[LINK TO ACCURACY ASSESSMENT FILES")
CS_lu = read.csv("[LINK TO CROPSCAPE LOOKUP")[,1:3]
Accuracy = Accuracy %>% reshape2::melt(id = c("X"))
names(Accuracy) = c("State", "Year", "Accuracy")
Accuracy$Year = str_remove(Accuracy$Year, "X")
Accuracy$State = as.character(as.factor(Accuracy$State))

## join the processed data with the State
df.grouped.sum = left_join(df.grouped.sum, WD_ST, 
                           by = c("ID" = "NAME"))

## Now join by state and year for accuracy assessment
df.grouped.acc = left_join(df.grouped.sum, Accuracy, 
                           by = c("State", "Year"))
df.grouped.acc$bound = (100 - df.grouped.acc$Accuracy)/100
df.grouped.acc = left_join(df.grouped.acc, CS_lu, 
                           by = c("value" = "Value"))
df.grouped.acc$area_lb = df.grouped.acc$total -(df.grouped.acc$total * df.grouped.acc$bound)
df.grouped.acc$area_ub = df.grouped.acc$total +(df.grouped.acc$total * df.grouped.acc$bound)

## Now look at volatility in planting
df.sum.var = df.grouped.acc %>% ungroup() %>% 
  dplyr::select("ID", "value", "Year", "State",
                "Accuracy", "bound", "Category",
                "Site.ID","System",
                "Use", "total", "area_lb", "area_ub")

df.sum.agg = df.sum.var %>% group_by(ID, Year, value, 
                                     Category, Use, System, Site.ID) %>%
  dplyr::summarize(Area.adj = sum(area_ub))
names(df.sum.agg) = c("ID", "Year", "value",
                      "crop", "Use", "System", "Site.ID", "Area.adj")

df.sum.agg = read.csv("[BIAS CORRECTED OUTPUT FILES HERE")
################################################################################
# ++++++++++++++++++++++++++++ Uncertainty Assessment +++++++++++++++++++++++++#
################################################################################
CS_BOR_LU = read.csv("[LOOKUP FILE BETWEEN USBR AND CROPSCAPE]",
                     stringsAsFactors = FALSE)
CS_BOR_LU = CS_BOR_LU %>% reshape2::melt(id = c("Value", "CropScape.Category", "Use")) 
Crop_LU_1 = read.csv("[LOOKUP FILE BETWEEN USBR AND CROPSCAPE #1]")[,-1]
CS_BOR_LU = left_join(CS_BOR_LU, Crop_LU_1, by = c("value"))
CS_BOR_LU = na.omit(CS_BOR_LU)
CS_BOR_LU = as.data.frame(CS_BOR_LU)

CS_BOR_LU = CS_BOR_LU %>% mutate(variable =  plyr::revalue(variable,
                                                           c("Costilla.Livestock.Co.Op" = "Costilla.Livestock.Co-Op")))
CS_BOR_LU = CS_BOR_LU %>% dplyr::filter(Use == "Ag") %>%
  dplyr::select("Value", "CropScape.Category", "Use", "variable", "Crop")
names(CS_BOR_LU) = c("value", "crop", "Use", "Site.ID", "BoR.crop")
CS_BOR_LU$value = as.numeric(as.character(CS_BOR_LU$value))

## Now join the BoR categories with CropScape
data.CS2BR = left_join(df.sum.agg, CS_BOR_LU, 
                       by = c("Site.ID", "Use", "value"))

## Filter out the Ag uses
data.CS2BR = data.CS2BR%>% dplyr::filter(Use == "Ag")
QA.QC = data.CS2BR%>% dplyr::filter(is.na(BoR.crop))

## Now look at the variability in land use
m2_2_acres = 0.000247105
data.CS2BR = data.CS2BR %>% dplyr::select(ID, Year, value, Use,
                                          crop.x, System, 
                                          Site.ID, Area.adj, BoR.crop)

data.CS2BR.sum = data.CS2BR %>% group_by(ID, System,Year, BoR.crop,
                                         Site.ID) %>%
  dplyr::summarise(Area.adj = sum(Area.adj))


write.csv(data.CS2BR.sum,
          "[LINK TO LANDUSE UNCERTAINTY HERE]")
