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
options(scipen=999)

hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}

source("D:/Projects/Project_Mother/Scripts/kings.R")
source("D:/Projects/Project_Mother/Scripts/map.R")

###############################################################################
# ++++++++++++++++++ Upload your data for processing +++++++++++++++++++++++++#
###############################################################################
## SKIP THIS IF YOU'RE CONFIDENT IN DATA
# read in your data
data45 = data.table::fread("[LINK TO PROCESSED 4.5 STREAMFLOW MANUALLY CONSTRUCTED OR OBTAINED FROM DATA FOLDER]")[,-1]
names(data45) = c("dates", "site", "Q_cms_mod45", "GCM", "Dam_BoB")
data85 = data.table::fread("[LINK TO PROCESSED 8.5 STREAMFLOW MANUALLY CONSTRUCTED OR OBTAINED FROM DATA FOLDER]")[,-1]
names(data85) = c("dates", "site", "Q_cms_mod85", "GCM", "Dam_BoB")
data.RCP = left_join(data45, data85, by= c("dates", "site", "GCM", "Dam_BoB"))

## 1. Break your data up into the periods of interest: 1) 2020-2050; 2) 2050-2080;
## 3) 2080-2100
data.RCP = data.RCP[data.RCP$dates >= as.Date("2020-10-01"),]
data.RCP$Period = ifelse(data.RCP$dates >= as.Date("2020-10-01") & 
                           data.RCP$dates <= as.Date("2050-09-30"), "2020-2050",
                         ifelse(data.RCP$dates >= as.Date("2050-10-01") & 
                                  data.RCP$dates <= as.Date("2080-09-30"), "2050-2080",
                                "2080-2100"))

## 1a. Obtain the monthly sum totals to match the observed/simulated climatology
data.RCP$WY = water_year(data.RCP$dates, "usgs")
data.RCP$month = month(data.RCP$dates)

## Sanity check
check = data.RCP %>% group_by(site, GCM) %>% tally() 

## Set up your Julian dates 
julian = seq(1,365,1)
julian.breaks = rep(seq(1,365,14), each = 14)
julian.breaks = as.data.frame(julian.breaks)
julian.breaks = julian.breaks[1:365,]
julian.breaks = as.data.frame(julian.breaks)
julian.breaks$waterdayjulian = julian
julian.breaks[365,1]= 351

## Get your water day julian
data.RCP$dates = as.Date(data.RCP$dates)
data.RCP$waterdayjulian = hydro.day.new(data.RCP$dates)

## Now join based on your julian breaks
data.RCP = left_join(data.RCP, julian.breaks, by = c("waterdayjulian"))
data.RCP$twoweeks = ifelse(data.RCP$waterdayjulian == 366, 351, 
                           data.RCP$julian.breaks)

## Store to save yourself data procesing once all of this is finalized
write.csv(data.RCP,
          "[LINK TO STORAGE LOCATION FOR ALL DATA")

###############################################################################
# +++++++++++++++ Second Time running the code start +++++++++++++++++++++++++#
###############################################################################
## If you've already run the code then just use
data.RCP = data.table::fread("[LINK TO PROCESSED DATA MANUALLY CONSTRUCTED OR FROM DATA FOLDER")

###############################################################################
# ++++++++++++Now following the SCF methodology ++++++++++++++++++++++++++++++#
###############################################################################
## 1. Get the median value for each 2 week period 
## a. Units conversion for cms to cmd
cms2cmd = 86400

## Calculate the monthly median in cubic meters for each period for each site
data.RCP.sub = data.RCP %>% 
  dplyr::group_by(GCM, site,twoweeks, Period) %>%
  dplyr::summarize(Q_cm_mod45 = median(Q_cms_mod45*cms2cmd),
                   Q_cm_mod85 = median(Q_cms_mod85*cms2cmd))

## Sanity check
check = data.RCP.sub %>% group_by(site, GCM) %>% tally()

## 1b. Convert your data into long format to make the next couple steps a bit
## easier
data.RCP.sub.m = data.RCP.sub %>% 
  reshape2::melt(id = c("GCM", "site", "twoweeks","Period"))

data.RCP.sub.m$RCP = ifelse(data.RCP.sub.m$variable == c("Q_cm_mod45"),
                            "RCP_45", "RCP_85") 

data.RCP.sub.m = data.RCP.sub.m %>% dplyr::select(-c("variable"))
names(data.RCP.sub.m) = c("GCM", "site", "twoweeks", "Period", "Q_cm_mod_proj", "RCP")

## 2. Join with historical climatology
## Read in your modeled historical
## First add your lookup table to reference across obs and mod
obs.mod.lu = read.csv("[LINK TO LOOKUP TABLE MODS_OBS_LU IN LOOKUP FOLDER]",
                      stringsAsFactors = FALSE)
data.RCP.sub.m = left_join(data.RCP.sub.m, obs.mod.lu, by = c("site")) ## join by site

obs.clim = data.table::fread("[LINK TO OBSERVED CLIMATOLOGY OBTAINED FROM *_PROCESSING.r]")
mod.clim = data.table::fread("[LINK TO SIMULATED CLIMATOLOGY OBTAINED FROM *_PROCESSING.r]")
check = mod.clim %>% group_by(Point, GCM) %>% tally()

## b. Join the dataframes
data.full = left_join(data.RCP.sub.m, mod.clim, by = c("twoweeks", "GCM", "System", "Point"))%>%
  left_join(obs.clim, by= c("twoweeks", "System", "Point"))
check = data.full %>% group_by(Point, GCM, RCP) %>% tally()

## 3. Multiplicative Change factor
data.full$CF_m = (data.full$Q_cm_mod_proj/data.full$Q_cm_mod)*data.full$Q_cm_obs

## 4. Clean up your data to speed things up and check
rm(list=ls()[! ls() %in% c("data.full", "hydro.day.new", "theme_kings", "theme_map")])

###############################################################################
# ++++++++++++ Here you will get your system metrics +++++++++++++++++++++++++#
###############################################################################
## 1. First break out your data into the cumulative sums for each period,
## RCP, Point, and GCM per two weeks
data.full.point = data.full %>% #dplyr::filter(Point == "Bridgeport") %>%
  dplyr::select(twoweeks, GCM, System, Point,Q_cm_obs, Period, RCP, CF_m) %>%
  group_by(GCM, Period, RCP,System, Point) %>%
  dplyr::mutate(Q_cm_obs_running = cumsum(Q_cm_obs),
                Q_cm_CF_m_running = cumsum(CF_m)) %>%
  slice(which.max(twoweeks))%>%
  dplyr::select(twoweeks, GCM, System, Point, Period, RCP, Q_cm_CF_m_running) %>%
  na.omit()
names(data.full.point) = c("twoweeks", "GCM", "System", "Point", "Period", "RCP", "Q_cm_point")

## 2. Now estimate the contribution for aggregate total for each system
data.full.system = data.full.point %>% #dplyr::filter(Point == "Bridgeport") %>%
  group_by(twoweeks, GCM, Period, RCP,System) %>%
  dplyr::mutate(Q_cm_total = sum(Q_cm_point))

## 3. Now get your point contribution for each system
data.full.system$frac_contribution = data.full.system$Q_cm_point/data.full.system$Q_cm_total
data.full.system = data.full.system %>% dplyr::filter(Period == "2080-2100") %>% ungroup()%>%
  dplyr::select(GCM, RCP, System, Point, frac_contribution)

## 4. Clean up your data to speed things up and check
rm(list=ls()[! ls() %in% c("data.full", "hydro.day.new", "theme_kings", "theme_map",
                           "data.full.system")])

###############################################################################
# ++++++++++++ Now look at your fractional contribution ++++++++++++++++++++++#
###############################################################################
## 1. First break out your data into the cumulative sums for each period,
## RCP, Point, and GCM per two weeks
data.full.sort = data.full %>% #dplyr::filter(Point == "Bridgeport") %>%
  dplyr::select(twoweeks, GCM, System, Point,Q_cm_obs, Period, RCP, CF_m) %>%
  group_by(GCM, Period, RCP,System, Point) %>%
  dplyr::mutate(Q_cm_obs_running = cumsum(Q_cm_obs),
                Q_cm_CF_m_running = cumsum(CF_m)) 

## 2. Join that record by your fractional weights obtained from above
data.full.sort.j= left_join(data.full.sort, data.full.system, 
                            by =c("GCM", "RCP", "System", "Point"))%>%
  dplyr::mutate(Q_cm_obs_running.f = Q_cm_obs_running*frac_contribution,
                Q_cm_CF_m_running.f =Q_cm_CF_m_running*frac_contribution)

## 3. Now merge at the system level for your scenario selection
data.full.sys = data.full.sort.j %>% na.omit() %>% ungroup()%>%
  group_by(GCM, Period, RCP, System, twoweeks) %>% 
  dplyr::summarize(Q_cm_obs_running.sys = sum(Q_cm_obs_running.f),
                   Q_cm_CF_m_running.sys = sum(Q_cm_CF_m_running.f))

## 4. Get your timing metrics for plotting
## Q50 & Annual
data.Q50.cor = data.full.sys %>%
  group_by(GCM, Period, RCP,System) %>%
  dplyr::mutate(Q_cm_CF_m_annual = max(Q_cm_CF_m_running.sys),
                Q_cm_obs_annual = max(Q_cm_obs_running.sys),
                Q_percent_m= Q_cm_CF_m_running.sys/Q_cm_CF_m_annual,
                Q_percent_obs= Q_cm_obs_running.sys/Q_cm_obs_annual) %>%
  dplyr::select("twoweeks", "GCM", "System", "Period", "RCP", 
                "Q_percent_m", "Q_percent_obs") %>%
  reshape2::melt(id = c("RCP", "GCM", "twoweeks", "System", "Period"))%>%
  dplyr::filter(value >= 0.5)%>%
  group_by(RCP, GCM,System, Period,variable) %>% 
  dplyr::slice(which.min(twoweeks)) 

## 5a. Spread your DF
data.Q50.spread = data.Q50.cor %>% tidyr::spread(key = "variable", value = "twoweeks") %>%
  group_by(GCM, Period, RCP,System) %>%
  dplyr::summarize(across(everything(), ~ first(na.omit(.)))) %>%
  group_by(GCM, Period, RCP,System) %>%
  dplyr::summarise(Timing_change_m = Q_percent_m-Q_percent_obs)
  
## 6. Get your timing metrics for plotting
## Q50 & Annual
data.Q.cor = data.full.sys %>%
  group_by(GCM, Period, RCP,System) %>%
  dplyr::mutate(Q_cm_CF_m_annual = max(Q_cm_CF_m_running.sys),
                Q_cm_obs_annual = max(Q_cm_obs_running.sys)) %>%
  dplyr::select("twoweeks", "GCM", "System", "Period", "RCP", "Q_cm_CF_m_annual",
                "Q_cm_obs_annual") %>%
  reshape2::melt(id = c("RCP", "System",  "GCM", "twoweeks", "Period"))%>%
  group_by(RCP, GCM,System, Period,variable) %>% 
  dplyr::slice(which.max(twoweeks)) 


## 6a. Spread your DF
data.Q.spread = data.Q.cor %>% tidyr::spread(key = "variable", value = "value") %>%
  group_by(GCM, Period, RCP,System) %>%
  dplyr::summarize(across(everything(), ~ first(na.omit(.)))) %>%
  group_by(GCM, Period, RCP,System) %>%
  dplyr::summarise(Amount_change_m = ((Q_cm_CF_m_annual-Q_cm_obs_annual)/Q_cm_obs_annual)*100)


## 7. Now join your data and look at plotting and get your quantiles
data.Q.joined = left_join(data.Q50.spread, data.Q.spread, 
                          by = c("RCP", "GCM", "Period", "System")) %>%
  dplyr::select("System", "RCP", "GCM", "Period", "Timing_change_m", "Amount_change_m")

data.Q.EOC = data.Q.joined %>% dplyr::filter(Period == c("2080-2100"))
q = c(.25, .5, .75)
quantiles = data.Q.EOC %>% ungroup() %>% dplyr::select(-c("System")) %>%
  reshape2::melt(id = c("RCP", "GCM", "Period"))%>%  
group_by(variable)%>%
  dplyr::summarize(quant25 = quantile(value, probs = q[1]), 
            quant50 = quantile(value, probs = q[2]),
            quant75 = quantile(value, probs = q[3]))

## 8. Now do your scenario selection
quantiles.m = quantiles %>% reshape2::melt(id = c("variable")) 
names(quantiles.m) = c("metric", "variable", "value")
quantiles.m = quantiles.m %>% tidyr::unite("key", c("metric", "variable")) %>%
  tidyr::spread(key = "key", value = "value")

## 9. Try out your selection 
scenario.selection = data.Q.EOC
scenario.selection$Timing_change_m_quant50 = quantiles.m$Timing_change_m_quant50
scenario.selection$Amount_change_m_quant50 = quantiles.m$Amount_change_m_quant50
scenario.selection$simple_results = ifelse(scenario.selection$Timing_change_m > scenario.selection$Timing_change_m_quant50 & 
                                             scenario.selection$Amount_change_m <= scenario.selection$Amount_change_m_quant50, 
                                           "Less/Neutral",
                                           ifelse(scenario.selection$Timing_change_m > scenario.selection$Timing_change_m_quant50 & 
                                                    scenario.selection$Amount_change_m > scenario.selection$Amount_change_m_quant50, 
                                                  "More/Neutral",
                                                  ifelse(scenario.selection$Timing_change_m <= scenario.selection$Timing_change_m_quant50 & 
                                                           scenario.selection$Amount_change_m <= scenario.selection$Amount_change_m_quant50, 
                                                         "Less/Earlier", "More/Earlier")))


## 10. Now plot (likely extract an example)
windows()
scenario.plot = ggplot(data = scenario.selection) +geom_point(aes(x= Amount_change_m, y = Timing_change_m, 
                                                  fill = simple_results),
                                      shape = 21, size = 6, alpha = 0.75)+
  geom_hline(aes(yintercept = Timing_change_m_quant50),size = 1, linetype = "dashed",
             color = "grey50")+
  geom_vline(aes(xintercept = Amount_change_m_quant50),size = 1, linetype = "dashed",
             color = "grey50")+
  scale_y_continuous(limits = c(-75, 75), n.breaks = 3)+
  scale_x_continuous(limits = c(-300,300), n.breaks =5 )+
  ylab("Change in DOQ50 Timing, 14 days")+xlab("Change in Annual Streamflow, %")+
  scale_fill_manual("Future Streamflow Scenario", values = c("darkred", "brown1", "darkblue", "dodgerblue3"))+
  theme_minimal(base_size =30)+facet_wrap(~System)+theme(legend.position = "bottom",
                                                         legend.direction = "horizontal",
                                                         legend.title = element_text(size = 28))
ggsave(plot = scenario.plot, 
       "[LINK TO SCENARIO DEVELOPMENT NOT USED IN MAIN MANUSCRIPT]", 
       width=21, height=10.5, dpi=300)

## 11. Now write your scenario selection to a specific file
scenario.selection.ex = scenario.selection %>% dplyr::select(System, RCP, GCM, simple_results) %>%
  ungroup()

write.csv(scenario.selection.ex,
          "[LINK TO SCENARIO DEVELOPMENT NOT USED IN MAIN MANUSCRIPT]")

## 11a. Quick table of results
scenario.selection.table = scenario.selection.ex 
scenario.selection.table$count = rep(1)

scenario.selection.table.risk = scenario.selection.table %>% 
  group_by(System, simple_results) %>%
  dplyr::summarize(scenario.sum = sum(count),
                   scenario.risk = (scenario.sum/64)*100) %>% ungroup() %>%
  complete(simple_results, System) %>%
  na_replace(0)

table.spread = scenario.selection.table.risk %>% group_by(System) %>%
  tidyr::spread(key = simple_results, value = scenario.risk)%>%
  dplyr::summarize(across(everything(), ~ first(na.omit(.)))) %>%
  dplyr::select(-c("scenario.sum"))

write.csv(table.spread,
          "[LINK TO SCENARIO DEVELOPMENT NOT USED IN MAIN MANUSCRIPT]")

scenario.selection.pairs = scenario.selection.table %>% 
  group_by(System, simple_results, GCM) %>%
  dplyr::summarize(scenario.sum = sum(count),
                   scenario.risk = (scenario.sum/(32))*100) %>% ungroup() %>%
  complete(simple_results, GCM) %>%
  na_replace(0)

table.spread = scenario.selection.pairs %>% group_by(GCM) %>%
  tidyr::spread(key = simple_results, value = scenario.sum)%>%
  dplyr::summarize(across(everything(), ~ first(na.omit(.)))) %>%
  dplyr::select(-c("scenario.risk"))

test = table.spread %>% reshape2::melt(id = "GCM") %>%
  group_by(GCM) %>%
  dplyr::summarize(count = sum(value))

ggplot()+geom_col(data = scenario.selection.pairs, aes(x = simple_results, y = scenario.sum, 
                                           group = GCM, fill = GCM), position = "dodge",
                  stat = "identity")+
  theme_kings(base_size = 22) + xlab("Future Streamflow Scenario")+
  ylab("Total Count across Systems")+theme(legend.position = "none")

write.csv(table.spread,
          "[LINK TO SCENARIO DEVELOPMENT NOT USED IN MAIN MANUSCRIPT]")

## 12. Join, sort, and plot your cumulative scenarios
target = c("Bridger", "Costilla", "Kern", "Little Wood", 
           "Paonia", "Price", "Shoshone", "Success", "Sun", "Umatilla", 
           "Walker", "Wind", "Yakima")

data.scenario.sys = data.full.sort %>% ungroup()%>%
  dplyr::filter(System %in% target) %>%
  group_by(twoweeks, GCM, Period, RCP, System) %>% 
  dplyr::summarize(Q_cm_obs_running.sys = sum(Q_cm_obs_running, na.rm = TRUE),
                   Q_cm_CF_m_running.sys = sum(Q_cm_CF_m_running, na.rm = TRUE))

## Join with your scenarios from above
data.scenario.initial = left_join(data.scenario.sys, scenario.selection.ex, 
                                  by = c("System", "RCP", "GCM"))

## Now get the median value per scenario
data.scenario = data.scenario.initial %>% #dplyr::filter(Point == "Bridgeport") %>%
  group_by(twoweeks, System, simple_results, Period.x) %>%
  dplyr::summarize(Q_cm_obs = median(Q_cm_obs_running.sys),
                Q_cm_sim_scen = median(Q_cm_CF_m_running.sys)) 

write.csv(data.scenario, "[LINK TO SCENARIO DEVELOPMENT NOT USED IN MAIN MANUSCRIPT]")

ggplot(data = data.scenario[data.scenario$System == "Walker",]) + 
  geom_line(aes(x= twoweeks, y = Q_cm_sim_scen,
                                           color = simple_results, group = simple_results),
                     size = 1) +
  geom_line(aes(x= twoweeks, y = Q_cm_obs,
                                  group = simple_results),
            color = "black", size = 1, linetype = "dashed")+
  theme_kings(base_size = 30) + facet_wrap(.~System+Period.x, ncol = 3)+
  scale_color_manual("Future Streamflow Scenario", values = c("darkred", "brown1", "darkblue", "dodgerblue3"))+
  ylab("Median cumulative streamflow (m3)")+xlab("Day of the Water Year")