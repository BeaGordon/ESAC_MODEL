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
options(scipen = 999)

###############################################################################
# +++++++++++++++ulate your modeled and observed climatology ++++++++++++++++ #
###############################################################################
hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}
##++++++++++++++++++++++++++++++Observed Climatology +++++++++++++++++++++++++##
## 1. read in your data
data = data.table::fread("[LINK TO HISTORICAL SWE DATA]")
data$date = as.Date(data$date)

## Convert to long format
data = data[,-1] %>% reshape2::melt(id = c("date"))
names(data) = c("dates", "Site", "SWE_mm")

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
## Get your system data
CA.area = read.csv("[LINK TO CA_AREA.CSV FROM LOOKUP FOLDER]") 

data.j = data.j %>%
  left_join(CA.area, by = c("Site" = "site")) %>%
  na.omit()

## Convert to volume
mm_2_m = 1/1000
data.j$SWE_m3 = data.j$SWE_mm * mm_2_m * data.j$Shape_Area

data.tw.sum = data.j %>% 
  group_by(System, dates) %>%
  dplyr::mutate(SWE_m3_sys = sum(SWE_m3)) %>%
  dplyr::select(System, twoweeks, SWE_m3_sys) %>%
  ungroup() %>% 
  group_by(System, twoweeks) %>%
  dplyr::summarize(SWE_m3 = median(SWE_m3_sys))

data.tw.sum = data.j %>% 
  dplyr::group_by(Site,twoweeks) %>% ## Change to Two Weeks if doing daily
  dplyr::summarize(SWE_obs = median(SWE_mm, na.rm = TRUE))

write.csv(data.tw.sum, 
          "[LINK TO PROCESSED SWE CLIMATOLOGY]")
