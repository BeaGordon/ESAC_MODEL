library(tidyr)
library(lubridate)
library(rlist)
library(reshape2)
library(dplyr)
library(hydroGOF)
library(lfstat)
library(plyr)
library(imputeTS)
library(ggplot2)
library(stringr)
library(extrafont)
library(usmap)
library(cluster)
library(factoextra)
library(wesanderson)
library(ggpubr)
library(BurStMisc)
library(purrr)
library(lwgeom)
library(sf)
library(ggpmisc)
library(polynom)
library(zoo)
library(smwrBase)

###############################################################################
# ++++++++++++++++++++ READ in REQ FUNCTIONS++++++++++++++++++++++++++++++ #
###############################################################################
hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}

# source("D:/Projects/Project_Mother/Scripts/kings.R")
# source("D:/Projects/Project_Mother/Scripts/map.R")

###############################################################################
# ++++++++++++++++++++ READ in HISTORICAL BF ++++++++++++++++++++++++++++++ #
###############################################################################

## 1. Set the directory to the Monthly_data file obtained from 
dir = "[READ IN YOUR USBR DATA OBTAINED FOR EACH SITE IN ITS OWN FOLDER]"
setwd(dir)

## 2. Create a list of all files matching the GCM name in the SWE directory
list = list.files(path = dir, 
                  pattern=glob2rx("*[PERIOD OF INTEREST]*"), # Input GCM name here
                  recursive = TRUE)
lapply.dat = lapply(list, read.csv, skip = 4)# Read in the list as .csvs

names = sub(".*0_", "", list)
names = str_remove( names, ".csv")
lapply.dat = Map(cbind, lapply.dat, ID = names)

## 3. Reassign column names to remove 'XHUC2'
daily.data = list.rbind(lapply.dat) # create a dataframe of all lists
header.data = colnames(daily.data)

## 4. Complete BOR Col names for each crop
headers = lapply(list, read.csv, skip = 2)
headers.sub = data.frame(t(headers[[1]][1,]))
headers.sub = apply(headers.sub, 2, function(x) gsub("^$|^ $", NA, x))
headers.sub = headers.sub %>% as.data.frame() %>% mutate_all(na_if,"")
sub.lu = na.omit(headers.sub)
sub.lu$number = seq(1, nrow(sub.lu), 1)
headers.sub.join = left_join(headers.sub, sub.lu, by = c("X1"))
headers.sub.c = na_locf(headers.sub.join$number)
headers.new = as.data.frame(headers.sub.c)
names(headers.new) = "number"
headers.new = left_join(headers.new, sub.lu)
headers.new$X2 = headers.sub.c
headers.new$X1 = str_remove(headers.new$X1, "Crop: ")
headers.new$X1 = gsub('[0-9]+', '', headers.new$X1)
headers.new.t = t(headers.new)
bor.colnames = headers.new.t[2,]
colnames(daily.data) = paste(bor.colnames, header.data, sep = "_")
colnames(daily.data) = gsub(" ", "", colnames(daily.data), fixed = TRUE)
names(daily.data)[ncol(daily.data)] <- "ID"

## Now melt after creating a date column
daily.data.sub = daily.data%>%
  mutate(dates = make_date(crops_Year,crops_Mo, crops_Dy))
daily.data.sub = daily.data.sub[,-c(1,3:4)]

## Now format by melting
daily.data.sub.melt = daily.data.sub %>% reshape2::melt(id = c("ID", "dates", "crops_DoY"))
daily.data.sub.melt$variable = gsub("[[:digit:]]+","",daily.data.sub.melt$variable)
daily.data.sub.melt$variable = sub("--", "_", daily.data.sub.melt$variable)
daily.data.sub.melt$variable = sub("-", "_", daily.data.sub.melt$variable)

daily.data.sub.melt$crop = sub("_.*", "", daily.data.sub.melt$variable)
daily.data.sub.melt$ET_var = sub(".*_", "", daily.data.sub.melt$variable)
daily.data.sub.melt$ET_var = gsub('[[:punct:] ]+','',daily.data.sub.melt$ET_var)

## Here is where you need to group by the 14 day total
## Set up your Julian dates 
julian = seq(1,365,1)
julian.breaks = rep(seq(1,365,14), each = 14)
julian.breaks = as.data.frame(julian.breaks)
julian.breaks = julian.breaks[1:365,]
julian.breaks = as.data.frame(julian.breaks)
julian.breaks$waterdayjulian = julian
julian.breaks[365,1]= 351

## Get your water day julian
daily.data.sub.melt$dates = as.Date(daily.data.sub.melt$dates)
daily.data.sub.melt$waterdayjulian = hydro.day.new(daily.data.sub.melt$dates)
daily.data.sub.melt$WY = waterYear(daily.data.sub.melt$dates)

## Now join based on your julian breaks
daily.data.sub.melt = left_join(daily.data.sub.melt, julian.breaks, by = c("waterdayjulian"))
daily.data.sub.melt$twoweeks = ifelse(daily.data.sub.melt$waterdayjulian == 366, 351, 
                                      daily.data.sub.melt$julian.breaks)

## Now get an average2020 in DoWY
daily.data.ave.act = daily.data.sub.melt %>% ungroup() %>% 
  dplyr::filter(ET_var == c("ETact"))%>% ungroup() 

daily.data.ave.act$value = ifelse(daily.data.ave.act$crop == c("Orchards"),
                                  daily.data.ave.act$value/2, daily.data.ave.act$value)
daily.data.ave.act =daily.data.ave.act %>%
  group_by(ID, twoweeks, ET_var, crop,WY) %>%
  dplyr::summarize(ET.act.sum = sum(value)) 

write.csv(daily.data.ave.act, 
          "[LINK TO STORE OF THE PROCESSED FILE]")

daily.data.ave.sum = daily.data.ave.act %>% ungroup() %>% 
  group_by(twoweeks, ET_var, crop) %>%
  dplyr::summarize(mean.ETact = median(ET.act.sum),
                   upper.bound.ETact = quantile(ET.act.sum, 0.75), 
                   lower.bound.ETact = quantile(ET.act.sum, 0.25))

write.csv(daily.data.ave.sum, 
          "[LINK TO STORE OF THE 14 DAY AGGREGATED PROCESSED FILE]")

