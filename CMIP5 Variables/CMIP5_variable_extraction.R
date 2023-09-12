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

###############################################################################
# ++++++++++++++++++++ READ in REQ FUNCTIONS++++++++++++++++++++++++++++++ #
###############################################################################
hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}

##############################################################################################################
# +++++++++++++++++++++++++++++++++ Set the working directory  +++++++++++++++++++++++++++++++++++++++++++++ #
##############################################################################################################
#RM temp files
tmp_dir <- tempdir()
list.files(tmp_dir)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

##############################################################################################################
# +++++++++++++++++++++++++++++++++++++++ Read in csvs and set up code++++++++++++++++++++++++++++++++++++++ #
##############################################################################################################
# this is where the  data is sitting 
dir <- "[LINK TO DIRECTORY FOR RAW FILES]"
setwd(dir)

## Read in the files if not uploading from R data
list = list.files(path=setwd(dir),
                  pattern=glob2rx("*[VARNAME]*.csv"), recursive = TRUE)

## Read in a .csv
lapply.dat = lapply(list, read.csv)

## Set up your names
model.name = sub("/.*", "", list)
lapply.dat = Map(cbind, lapply.dat, Model = model.name)
data= list.rbind(lapply.dat) # create a dataframe of all lists

data$dates = as.Date(data$dates)
data$DoY = format(data$dates, "%j")

write.csv(data, "[LINK TO PROCESSED DIRECTORY FILE")




