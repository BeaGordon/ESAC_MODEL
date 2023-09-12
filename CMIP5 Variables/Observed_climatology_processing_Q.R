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
library(utils)
library(R.utils)
library(ncdf4)
library(ncdf4.helpers)
library(raster)
library(tidync)
library(hydroGOF)

###############################################################################
# +++++++++++++++ulate your modeled and observed climatology ++++++++++++++++ #
###############################################################################
hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}
Site_LU = read.csv("[LINK TO MODS_OBS_LU.csv from the LOOKUP TABLE FILE]")
##++++++++++++++++++++++++++++++Modeled Climatology +++++++++++++++++++++++++##
## 1. read in your data
data = data.table::fread("[LINK TO PROCESSED STREAMFLOW MANUALLY CONDUCTED OR PROVIDED IN RASTER FOLDER]")
data = data[data$dates > as.Date("[BEGIN DATE]") & data$dates < as.Date("[END DATE]"), ] ## sub out the time frame
data$dates = as.Date(data$dates)
data$WY = water_year(data$dates, "usgs") ## add USGS WY
data$val = ifelse(data$val < 0, 0, data$val)

## 2. Get the two week values
julian = seq(1,365,1)
julian.breaks = rep(seq(1,365,14), each = 14)
julian.breaks = as.data.frame(julian.breaks)
julian.breaks = julian.breaks[1:365,]
julian.breaks = as.data.frame(julian.breaks)
julian.breaks$waterdayjulian = julian
julian.breaks[365,1]= 351

## Get your water day julian
data$dates = as.Date(data$dates, format = "%Y-%m-%d")
data$waterdayjulian = hydro.day.new(data$dates)
data$WY = water_year(data$dates, "usgs")

## Now join based on your julian breaks
data.j = left_join(data, julian.breaks, by = c("waterdayjulian"))
data.j$twoweeks = ifelse(data.j$waterdayjulian == 366, 351, 
                         data.j$julian.breaks)

## 5. Get the Median values
## Get the m3 total by converting cfs to cfd
cfs_2_cms = 0.0283168
cms_2_cmd = 86400

data.tw.sum = data.j %>% 
  dplyr::group_by(System,Point,twoweeks) %>% ## Change to Two Weeks if doing daily
  dplyr::summarize(Q_cm_obs = median(val*cfs_2_cms*cms_2_cmd, na.rm = TRUE))

# data.tw.med = data.j %>% 
#   dplyr::group_by(System, Point,twoweeks, WY) %>% ## Change to Two Weeks if doing daily
#   dplyr::summarize(Q_cm_obs = sum(val*cfs_2_cms*cms_2_cmd, na.rm = TRUE))

write.csv(data.tw.sum, 
          "[LINK TO OBSERVED CLIMATOLOGY FOLDER")
