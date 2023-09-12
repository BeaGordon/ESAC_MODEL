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

###############################################################################
# ++++++++++++++++++++ READ in HISTORICAL BF ++++++++++++++++++++++++++++++ #
###############################################################################

## 1. Set the directory to the Monthly_data file obtained from 
## Preprocessing_CADMUS code
dir = "[LINK TO EXTRACTED USBR FILES HERE]" 
setwd(dir)

## 2. Create a list of all files matching the GCM name in the SWE directory
list = list.files(path = dir, 
                  pattern=glob2rx("*.csv"), # Input GCM name here
                  recursive = TRUE)
list = list[!str_detect(list,pattern="[PATTERN HERE]")]
list = list[!str_detect(list,pattern="[PATTERN HERE]")]

lapply.dat = lapply(list, read.csv, skip = 4, stringsAsFactors = FALSE)# Read in the list as .csvs
keep = c("Year","DoY","Mo","Dy","PMETo")
lapply.dat = lapply(lapply.dat, function(x) subset(x, select = intersect(keep, colnames(x))))

## Create your variables/list names
## For scenarios
Scenario = sub(".*0_", "", list)
Scenario = str_remove( Scenario, ".csv")

## For sites
Site = list
Site = gsub("/", "_", Site, fixed=TRUE)
Site = word(Site, 1, sep = "_")

## For period
Period = list
Period = str_replace(Period, "S0", "baseline")
Period = str_replace(Period, "708_CA9452", "708CA9452")
Period = str_replace(Period, "W-14_WY6555", "W-14WY6555")
Period = str_replace(Period, "14060007_U-2_UT7026/U-2_UT7026", 
                    "14060007U-2UT7026/U-2UT7026")
Period = str_replace(Period, "/14020005_C-5_CO6306/C-5_CO6306", 
                    "/14020005C-5CO6306/C-5CO6306")
Period = sub("^[^_]*_", "", Period)
Period = str_remove(Period, "_S[0-9]+.csv")
Period = str_remove(Period, ".csv")

## Now assign to your DF
lapply.dat = Map(cbind, lapply.dat, Period = Period, Site = Site, Scenario = Scenario)

## 3. Reassign column names to remove 'XHUC2'
daily.data = list.rbind(lapply.dat) # create a dataframe of all lists

## 4. Now melt after creating a date column
daily.data.sub = daily.data%>%
  mutate(dates = make_date(Year,Mo, Dy))
daily.data.sub = daily.data.sub[,-c(1,3:4)]
## Now get the water day of the year
daily.data.sub$waterdayjulian = hydro.day.new(daily.data.sub$dates)

## 5. Now format by melting
daily.data.sub.melt = daily.data.sub %>% 
  reshape2::melt(id = c("DoY", "waterdayjulian", "dates", "Period", "Site", "Scenario"))

## Here is where you add your julian day grouping to match streamflow
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

###############################################################################
# ++++++++++++++++++++ Analyze the res demand, combine +++++++++++++++++++++++ #
###############################################################################
## 1. First you want to look at your upper and lower bounds for PET demand
data.ave.ETo = daily.data.sub.melt %>% ungroup() %>% select(-DoY) %>%
  group_by(Site, Period, variable, twoweeks, WY,Scenario) %>%
  dplyr::summarize(value = sum(value)) %>% ungroup() %>%
  group_by(Site, Period, variable, twoweeks) %>%
  dplyr::summarize(median.ETo = median(value), 
                   upper.bound.ETo = quantile(value, 0.75), 
                   lower.bound.ETo = quantile(value, 0.25))

## 7. Now join with the lookup table
WU_2_Res = read.csv("[LINK TO WU_2_RES LU TABLE IN LOOKUP FOLDER]",
                    stringsAsFactors = FALSE)
data.ETo = left_join(data.ave.ETo, WU_2_Res, by = c("Site"= "Site.Res")) %>%
  na.omit()

## 8. Now join with the reservoir area
Reservoir.Area = read.csv("[LINK TO RESERVOIR AREA IN RASTERS & DATA FOLDER]",
                          stringsAsFactors = FALSE)
data.ETo.Area = left_join(data.ETo,Reservoir.Area, by = c("Reservoir"))

## 2. Add your reservoir area + PET
## declare your constants
mm_2_m =0.00328084
km2_m2 = 1000000
k_mid = 0.5625
k_low  = 0.5625
k_high = 0.5625
cm = 1000000

## Calculate
data.ETo.calc = data.ETo.Area %>% 
  group_by(Reservoir, System, Period, twoweeks,CAP_MAX) %>%
  dplyr::summarize(Lower.bound.Use = (AREA_SKM*km2_m2)  * k_low*(lower.bound.ETo*mm_2_m),
                   Upper.bound.Use = (AREA_SKM*km2_m2)  * k_high*(upper.bound.ETo*mm_2_m),
                   Med.Use = (AREA_SKM*km2_m2) *k_mid* (median.ETo*mm_2_m))

## 3. Now obtain the cumulative sum
data.ETo.sys = data.ETo.calc %>% 
  group_by(Reservoir, System, Period,CAP_MAX) %>%
  dplyr::mutate(Lower.bound.Use.CS = cumsum(Lower.bound.Use),
                Upper.bound.Use.CS = cumsum(Upper.bound.Use),
                Med.Use.CS = cumsum(Med.Use),
                Maximum_capacity = CAP_MAX*cm) %>% na.omit() %>%
  group_by(System, Period, twoweeks) %>%
  dplyr::summarize(twoweeks = twoweeks,
                   Maximum_capacity = sum(Maximum_capacity),
                   LowLU.LowWU = sum(Lower.bound.Use.CS),
                   HighLU.HighWU = sum(Upper.bound.Use.CS),
                   MedLU.MedWU = sum(Med.Use.CS)) %>% unique()
  
## 4. Now format
data.ETo.sys$Period = ifelse(data.ETo.sys$Period == c("2020"), "2020-2050", 
                     ifelse(data.ETo.sys$Period == c("2050"), "2050-2080", 
                            ifelse(data.ETo.sys$Period == c("2080"), "2080-2100",
                                   "baseline")))
data.ETo.sys = data.ETo.sys %>% ungroup()%>%
  dplyr::select(System, Period, twoweeks, Maximum_capacity, LowLU.LowWU, HighLU.HighWU, MedLU.MedWU) %>%
  reshape2::melt(id = c("System", "Period", "twoweeks", "Maximum_capacity"))

names(data.ETo.sys) = c("System", "Period", "twoweeks", "Maximum_capacity", "Scenario",
                        "ETo_m3")

## 5. Update for consistency with other demand
data.ETo.sys$Scenario = str_replace(data.ETo.sys$Scenario, "HighLU.HighWU", "High Water Use")
data.ETo.sys$Scenario = str_replace(data.ETo.sys$Scenario, "MedLU.MedWU", "Business-as-Usual")
data.ETo.sys$Scenario = str_replace(data.ETo.sys$Scenario, "LowLU.LowWU", "Low Water Use")

## 5. Filter out the parts that don't make sense 
data.ETo.sys$screen = ifelse(data.ETo.sys$Scenario == c("High Water Use") & 
                               data.ETo.sys$Period == c("baseline"), 0,
                                    ifelse(data.ETo.sys$Scenario == c("Low Water Use") & 
                                             data.ETo.sys$Period == c("baseline"), 0, 1))

data.ETo.sys$Scenario = ifelse(data.ETo.sys$Scenario == c("Business-as-Usual") & 
                                        data.ETo.sys$Period == c("baseline"), "Historical",
                                      data.ETo.sys$Scenario)

data.ETo.sys$Period = ifelse(data.ETo.sys$Scenario == c("Historical") & 
                               data.ETo.sys$Period == c("baseline"), "1950-1999",
                             data.ETo.sys$Period)

data.ETo.sys = data.ETo.sys %>% dplyr::filter(screen !=0)

data.capacity = data.ETo.sys %>% dplyr::select(System, Period, twoweeks, Maximum_capacity)
write.csv(data.capacity, 
          "[LINK TO HISTORICAL RESERVOIR EVAPORATION]")

## 6. Adjust categories 
data.ETo.sys = data.ETo.sys %>% dplyr::select(System, Period, twoweeks, Scenario, ETo_m3)

write.csv(data.ETo.sys, 
          "[LINK TO SYSTEM ETO RESULTS FOLDER")

data.historical = data.ETo.sys %>% dplyr::filter(Period == c("1950-1999")) %>%
  dplyr::filter(Scenario == c("Business-as-Usual"))

write.csv(data.historical, 
          "[LINK TO HISTORICAL ETO RESULTS FOLDER]")
