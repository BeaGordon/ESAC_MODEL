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
library(usmap)
library(RColorBrewer)
library(ggnewscale)
library(tidyr)
library(maps)
library(stars)

###############################################################################
# ++++++++++ Standard uploads for processing and functions +++++++++++++++++++#
###############################################################################
us = map_data("state")
wus = us %>% 
  dplyr::filter(region %in% c("arizona", "california",
                              "colorado", "idaho", "montana", "nevada",
                              "oregon", "utah", "washington", "wyoming",
                              "new mexico", "south dakota", "north dakota",
                              "nebraska", "kansas", "texas", "oklahoma"))


## Conversions and other things
cm_2_mcm = 1/1000000 
crs1 = crs(regions)

options(scipen=999)

hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}

sites = read.csv("[LINK TO SITE COORDS]")[,-c(1:2,4)]
sites$System = str_replace(sites$System, "Little Wood River", "Little Wood")
sites$System = str_replace(sites$System, " Shoshone", "Shoshone")
sites$System = str_replace(sites$System, "Yakima", "Kittitas")

states = c("Arizona", "California",
           "Colorado", "Idaho", "Montana", "Nevada",
           "Oregon", "Utah", "Washington", "Wyoming","New Mexico","South Dakota", 
           "North Dakota",
           "Nebraska", "Kansas", "Texas", "Oklahoma")

## Here is your elevation raster to make the maps pretty
elev = raster("[LINK TO ELEVATION FILE WITHIN RASTERS FOLDER]")
newproj = "+proj=longlat +datum=WGS84"
elev2 = projectRaster(elev, crs=newproj, res=0.1)

#> returning a list of RasterLayer objects
usa = spData::us_states %>% sf::st_transform(4326)
usa = usa[usa$REGION ==c("West"),] 
usa_elev = crop(elev2[[1]], usa)
rdf = as.data.frame(usa_elev, xy=TRUE) #Convert raster to data.frame
names(rdf)[3] = 'magnitude' #Name value column
rdf$magnitude = ifelse(rdf$magnitude < 500, NA, rdf$magnitude)
rdf$magnitude = ifelse(rdf$x < -100 & rdf$x > -120 & rdf$magnitude < 1300,
                       NA, rdf$magnitude)
rdf = na.omit(rdf)

## Regions for basin backdrop
regions = sf::st_read("[LINK TO BASINS SHAPEFILE FROM RASTERS FOLDER]")
regions =sf::st_transform(regions, "+proj=longlat +datum=WGS84")

## Rivers for map backdrop
rivers = sf::st_read("[LINK TO BASINS SHAPEFILE FROM RASTERS FOLDER]")
rivers =sf::st_transform(rivers, "+proj=longlat +datum=WGS84")

###############################################################################
# +++++++++++++++++++++++++ EXPSOSURE ANALYSIS  ++++++++++++++++++++++++++++++#
###############################################################################
## read in your processed data
data = read.csv("[READ IN PROCESSED DATA OBTAINED FROM ACCOMPANYING SCRIPTS 
                OR PROVIDED AS 'paper_data.csv' in RASTERS FOLDER")

# add color scale information
data = data %>%
  mutate(color =  plyr::revalue(System,
                                c(
                                  "Little Wood" = "dodgerblue4", 
                                  "Kittitas" = "cornflowerblue",
                                  "Umatilla" = "deepskyblue",
                                  # West-orangered
                                  "California " = "thistle4", 
                                  "Great Basin " = "thistle", 
                                  #Southwest-dark yellows
                                  "Bridger" = "orangered2",
                                  "Paonia" = "coral3", 
                                  "Price" = "lightcoral",
                                  #"Souris-Red-Rainy " = "lightsalmon1",
                                  # West north central-light yellows
                                  "Costilla" = "peachpuff", 
                                  #East north central-light greens
                                  "Kaweah" = "darkgoldenrod3",
                                  "Kern" = "lightgoldenrod1",
                                  # Central-medium greens
                                  "Sun" = "olivedrab3", 
                                  "Shoshone" = "lightseagreen", 
                                  "Wind" = "darkturquoise",
                                  ## Great Basin
                                  "Walker" = "slategray2" 
                                )))

myColorScale = unique(as.vector(data$color))
names(myColorScale) = unique(as.vector(data$System))
colScale = scale_fill_manual(name = "System",values = myColorScale)

## First set up your exposure DF
data.fig1 = data

## Next obtain the raw magnitude and timing differences
data.fig1$Amount_dif = (data.fig1$Supply/data.fig1$Supply_obs)
data.fig1$Timing_dif = abs(data.fig1$Supply_50-data.fig1$Supply_50_obs)

## In this mutation: 1) Group by System and Period to obtain a median mag and timing
## difference across all 64 GCM-RCP pairs; 2) use the median values across all systems
## and each period to get the minimum and maximum for re-scaling
data.fig1 = data.fig1 %>%
  group_by(System,Period) %>%
  dplyr::mutate(Timing = median(Timing_dif),
                Mag = median(Amount_dif),
                Supply_med = median(Supply),
                Supply_50_med = median(Supply_50),
                NaturalStorage_med = median(NaturalStorage),
                Demand_med = median(Demand),
                timing_delta = ((Supply_50_med-Supply_50_obs)/Supply_50_obs)*100,
                magnitude_delta = ((Supply_med-Supply_obs)/Supply_obs)*100)%>%
  ungroup() %>%
  dplyr::mutate(A_timing = min(Timing),
                B_timing = max(Timing),
                A_mag = min(Mag),
                B_mag = max(Mag))

## In this mutation: 1) Group by System and Period to obtain a Timing I and Magnitude I
## based on the median Timing and Mag values for all 64 GCM-RCP pairs from above
data.fig1 = data.fig1 %>% group_by(System, Period) %>%
  dplyr::mutate(Timing_I = 1+(10-1)*((Timing - A_timing)/(B_timing-A_timing)),
                Mag_I = 1+(10-1)*((Mag - A_mag)/(B_mag-A_mag)),
                Mag_I_radar = 11-Mag_I) 

## Here you get the exposure metric either as a geometric mean or as just an average
data.fig1$Exp = ((11-data.fig1$Mag_I) * (data.fig1$Timing_I))^(1/2)## Geometric mean

## Here you get the median change in storage across the GCM-RCP pairs
data.fig1 = data.fig1 %>% ungroup()%>% group_by(System, Period) %>%
  dplyr::mutate(storage_med = median((NaturalStorage-NaturalStorage_obs)/BuiltStorage),
                storage_med_snow = median((NaturalStorage-NaturalStorage_obs)/NaturalStorage_obs)) %>%
  left_join(sites, by = c("System")) %>% unique()

## Add in elevation
elevation = read.csv("D:/Projects/Project_Reservoir/Data/Spreadsheets/System Results/Elevation.csv")
median_elevation = median(elevation$Elevation)
elevation$Systems = str_replace(elevation$Systems, "Yakima", "Kittitas")
elevation$Systems = str_replace(elevation$Systems, "Little Wood River", "Little Wood")
elevation$Systems = str_replace(elevation$Systems, " Shoshone", "Shoshone")

## Create a dataframe including elevation and crop data
data.fig1 = left_join(data.fig1, crop.frac, by = c("System", "Period")) %>%
  left_join(elevation, by = c("System" = "Systems")) %>%
  left_join(hist.frac, by = c("System"))

## data for analysis
data.fig1.sub = data.fig1 %>% dplyr::select(System,Exp, Period,
                                            storage_med_snow,
                                            NaturalStorage_obs, BuiltStorage,
                                            Demand_obs, Supply_med, Supply_obs,Elevation,
                                            Demand_med, NaturalStorage_med, Flexible,
                                            Capital, Flexible.hist, Captial.hist,
                                            timing_delta, magnitude_delta) %>% ungroup() %>%
  group_by(Period) %>%
  dplyr::mutate(Snow_med = median(NaturalStorage_med)) %>% 
  ungroup() %>%
  group_by(System, Exp, Period)%>%
  dplyr::mutate(Hist_Supply_Demand = Demand_obs/Supply_obs, 
                Proj_Supply_Demand = Demand_med/Supply_med,
                Hist_Storage = NaturalStorage_obs/BuiltStorage,
                Proj_Storage = NaturalStorage_med/BuiltStorage,
                Hist_Snow_Supply = NaturalStorage_obs/Supply_obs,
                Proj_Snow_Supply = NaturalStorage_med/Supply_med,
                Delta_Supply_Demand = (Proj_Supply_Demand-Hist_Supply_Demand)/Hist_Supply_Demand,
                Delta_Storage = (Proj_Storage-Hist_Storage)/Hist_Storage,
                Delta_Snow_Supply = (Proj_Snow_Supply-Hist_Snow_Supply)/Hist_Snow_Supply,
                Delta_Supply = (Supply_med-Supply_obs)/Supply_obs,
                Delta_Demand = (Demand_med-Demand_obs)/Demand_obs,
                Delta_Snow = (NaturalStorage_med-NaturalStorage_obs)/NaturalStorage_obs)%>%
  dplyr::select(System, Exp, storage_med_snow, 
                Hist_Supply_Demand,Hist_Storage,Hist_Snow_Supply,
                Proj_Snow_Supply,Proj_Storage,Proj_Supply_Demand,
                Delta_Supply_Demand,Delta_Storage,Delta_Snow_Supply,Delta_Supply,
                Delta_Demand,Delta_Snow)%>% unique()

###############################################################################
# +++++++++++++++++++++++++ CODE FOR FIGURE 2  +++++++++++++++++++++++++++++++#
###############################################################################

######################### METRIC BASED MAP OF EXPOSURE ##########################
## set up your custom palette
exposure.pal = c('#67000d','#a50f15','#cb181d','#ef3b2c',
                 '#fb6a4a','#fc9272','#fcbba1','#fee0d2', '#fff5f0')

## set up your labels for plotting
data.fig1$System = ifelse(data.fig1$System == c("Little\nWood"), "Little Wood",
                          data.fig1$System)
unique(data.fig1$POINT_Y)
unique(data.fig1$POINT_X)
unique(data.fig1$System)
#Shift points
data.fig1$POINT_Y = ifelse(data.fig1$System == c("Shoshone"), 45,data.fig1$POINT_Y)

# Shift points
data.fig1$POINT_Y_lab = ifelse(data.fig1$System == c("Bridger"), 40.8,
                               ifelse(data.fig1$System == c("Costilla"), 34.8,
                                      ifelse(data.fig1$System == c("Kaweah"), 35.85,
                                             ifelse(data.fig1$System == c("Kern"), 33.0,
                                                    ifelse(data.fig1$System == c("Kittitas"), 47,
                                                           ifelse(data.fig1$System == c("Little Wood"), 41.5,
                                                                  ifelse(data.fig1$System == c("Paonia"),39.1,
                                                                         ifelse(data.fig1$System == c("Price"), 37.4,
                                                                                ifelse(data.fig1$System == c("Wind"),43.0,
                                                                                       ifelse(data.fig1$System == c("Sun"), 47.5,
                                                                                              ifelse(data.fig1$System == c("Umatilla"), 43.3,
                                                                                                     ifelse(data.fig1$System == c("Walker"), 38.5,
                                                                                                            45.1))))))))))))
data.fig1$POINT_X_lab = ifelse(data.fig1$System == c("Bridger"), -109.6,
                               ifelse(data.fig1$System == c("Costilla"), -112.1,
                                      ifelse(data.fig1$System == c("Kaweah"), -125.8,
                                             ifelse(data.fig1$System == c("Kern"), -124.5280,
                                                    ifelse(data.fig1$System == c("Kittitas"), -126.9,
                                                           ifelse(data.fig1$System == c("Little Wood"), -119.9,
                                                                  ifelse(data.fig1$System == c("Paonia"),-107.0,
                                                                         ifelse(data.fig1$System == c("Price"), -116.3,
                                                                                ifelse(data.fig1$System == c("Wind"),-108.3,
                                                                                       ifelse(data.fig1$System == c("Sun"), -112.2486,
                                                                                              ifelse(data.fig1$System == c("Umatilla"), -126.3,
                                                                                                     ifelse(data.fig1$System == c("Walker"), -126.3,
                                                                                                            -109))))))))))))


data.fig1 = data.fig1 %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(
                                    "Little Wood" = "LW-ID",
                                    "Sun" = "Su-MT",
                                    "Shoshone" = "Sh-WY", 
                                    "Wind" = "Wi-WY",
                                    "Bridger" = "Br-WY", 
                                    "Price" = "Pr-UT",
                                    "Paonia" = "Pa-CO", 
                                    "Costilla" = "Co-NM",

                                    "Kern" = "Ke-CA", 
                                    "Kaweah" = "Ka-CA",
                                    "Walker" = "Wa-NV", 
                                    "Umatilla" = "Um-OR",
                                    "Kittitas" = "Ki-WA"
                                    #Southwest-dark yellows
                                    #"Souris-Red-Rainy " = "lightsalmon1",
                                    # West north central-light yellows
                                    #East north central-light greens
                                    # Central-medium greens
                                  )))

########################## CODE FOR FIGURE 2A-C ###############################
pal = c("#01665e","#8c510a","#c7eae5","#5ab4ac" )
exposure.map = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), 
                                      fill = "white",
                                      color = "grey77", size = 0.35, alpha = 0.4)+
  geom_sf(data = regions, fill = "#f5f5f5", color = "grey50", alpha = 0.5) +
  geom_raster(data = rdf, aes(x=x, y=y, fill = magnitude), alpha = 0.4) +
  scale_fill_gradientn(colors = rev(pal))+
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_label(data = data.fig1,aes(x=POINT_X_lab, y=POINT_Y_lab,
                                  label = Systems), hjust = 0,vjust = -0.55, family = "Helvetica", 
             size = 7.3, fill = "#f7fcfd", label.size = 0, label.padding = unit(0, "lines"))+
  geom_point(data = data.fig1,
             aes(x=POINT_X, y=POINT_Y, fill =Exp, size =storage_med_snow*100), 
             shape = 21, color = "black", alpha = 1, stroke = 1.05)+
  scale_fill_stepsn("Exposure",
                    colors=rev(exposure.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-.25,-.5,-.75,-1.0))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+  
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top")) +
  facet_wrap(~Period, nrow = 1)+theme(panel.grid= element_blank(), 
                                      legend.position = "none")+
  theme(strip.text.x = element_text(
    face = "bold",
    margin = unit(rep(25, 4), "pt")
  ))

## Make a custom scale for your fill
exposure.map.fill = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                           color = "grey85", size = 0.25, alpha = 0.55)+
  scale_fill_gradient2(low = "white", high = "grey20")+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig1,
             aes(x=POINT_X, y=POINT_Y, fill =Exp), size = 5,
             shape = 21, color = "black", alpha = 1, stroke = 1)+
  scale_fill_stepsn("Relative Exposure",
                    colors=rev(exposure.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous("DOQ50 Timing [Days]",
                        range = c(5,15))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(size = 28, family  = "Helvetica"))+
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0))+
  guides(fill = guide_colorsteps(frame.colour = "black", ticks=TRUE,
                                 ticks.colour = "black",
                                 ticks.linewidth = 1,
                                 title.position = "top"))

exposure.legend = get_legend(exposure.map.fill+theme(legend.position = "right",
                                                     legend.title.align = 0.5))

## Make a custom scale for your size
exposure.map.size = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                           color = "grey85", size = 0.25, alpha = 0.55)+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig1,
             aes(x=POINT_X, y=POINT_Y, size =storage_med_snow*100),
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Exposure",
                    colors=rev(exposure.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", S[Snow])),
                        range = c(15,5), breaks = c(-15,-30,-45,-60),
                        labels = c("-15%", "-30%","-45%", "-60%"))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2.5, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top",
                             label.position = "bottom")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0))

exposure.legend.size = get_legend(exposure.map.size+theme(legend.position = "right",
                                                          legend.title.align = 0.5))

########################## CODE FOR FIGURE 2D-F ###############################
data.fig1.radar = data.fig1 %>% ungroup() %>%
  dplyr::select(System, Period,  timing_delta, magnitude_delta) %>%
  unique() %>% group_by(System, Period) %>%
  reshape2::melt(id = c("System", "Period"))


data.fig1.radar$variable = str_replace(data.fig1.radar$variable, "timing_delta", " \u0394 in Timing ")
data.fig1.radar$variable = str_replace(data.fig1.radar$variable, "magnitude_delta", "\u0394 in Magnitude")

data.fig1.radar = data.fig1.radar %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(
                                    "Little Wood" = "LW-ID         ",
                                    "Sun" = "         Su-MT",
                                    "Shoshone" = "         Sh-WY", 
                                    "Wind" = "         Wi-WY",
                                    "Bridger" = "         Br-WY", 
                                    "Price" = "         Pr-UT",
                                    "Paonia" = "         Pa-CO", 
                                    "Costilla" = "Co-NM",

                                    "Kern" = "Ke-CA         ", 
                                    "Kaweah" = "Ka-CA         ",
                                    "Walker" = "Wa-NV         ", 
                                    "Umatilla" = "Um-OR         ",
                                    "Kittitas" = "Ki-WA         "
                                    #Southwest-dark yellows
                                    #"Souris-Red-Rainy " = "lightsalmon1",
                                    # West north central-light yellows
                                    #East north central-light greens
                                    # Central-medium greens
                                  )))

data.fig1.radar$Systems = factor(data.fig1.radar$Systems, 
                                  c(
                                    "LW-ID         ",
                                    "         Su-MT",
                                    "         Sh-WY", 
                                    "         Wi-WY",
                                    "         Br-WY", 
                                    "         Pr-UT",
                                    "         Pa-CO", 
                                    "Co-NM",

                                    "Ke-CA         ", 
                                    "Ka-CA         ",
                                    "Wa-NV         ", 
                                    "Um-OR         ",
                                    "Ki-WA         "
                                    #Southwest-dark yellows
                                    #"Souris-Red-Rainy " = "lightsalmon1",
                                    # West north central-light yellows
                                    #East north central-light greens
                                    # Central-medium greens
                                  ))

df1 = data.fig1.radar %>% 
  dplyr::arrange((levels = Systems))

## Set up your radar plot for comparison
coord_radar = function (theta = "x", start = 0, direction = 1) 
{
  theta = match.arg(theta, c("x", "y"))
  r = if (theta == "x") 
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

## Now plot
exp.compass = ggplot(df1, aes(x=Systems, y=value, color=variable,group = variable),
                     alpha = 0.5) +
  geom_hline(aes(yintercept = 0), color = "red", alpha = 0.75)+
  geom_polygon(aes(color = variable, linetype = variable), fill = NA, size =1.75) +
  theme_minimal(base_size = 30)+
  coord_polar(start = 5.7, clip = "off")+
  scale_y_continuous(lim = c(-40,10))+
  scale_fill_manual("",
                    values = c("grey50","black" ))+
  scale_color_manual("Components",
                     values = c("grey50","black" ))+
  facet_wrap(Period ~ ., labeller = label_wrap_gen(width=14), nrow =1) +
  ylab("")+xlab("")+
  guides(colour=guide_legend(nrow=1, byrow=TRUE)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c("bottom"))+theme(legend.direction="horizontal")+
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
  theme(axis.text.x = element_text(colour="black", size=21,  hjust = -4),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "grey80", size = 0.75))+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing = unit(9, "lines"))+
  theme(legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+
  theme(legend.position = "none")+
  theme(plot.background = element_rect(color = 0,
                                       size = 0),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) +#
  
  annotate(geom ="label", x = 1.5, alpha = 0.7,
           y = seq(-40, 10, 10), 
           label = c("-40%", "-30%", "-20%", "-10%", "0%", "10%"),
           size = unit(8, "pt"), fill = "white",label.size=NA) +
  
  # remove original y-axis text / ticks & black border
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 

## Now make your legend plot
exp.compass.legend.plot = ggplot(df1, aes(x=Systems, y=value, color=variable,group = variable),
                                 alpha = 0.5) +
  geom_line(aes(color = variable, linetype = variable), fill = NA,  size =1.5) +
  theme_minimal(base_size = 30)+
  coord_polar(start = 5.7, clip = "off")+
  scale_y_continuous(lim = c(-40,10))+
  scale_fill_manual("",
                    values = c("grey50","black" ))+
  scale_color_manual("Components of Exposure",
                     values = c("grey50","black" ))+
  facet_wrap(Period ~ ., labeller = label_wrap_gen(width=14), nrow =1) +
  ylab("")+xlab("")+
  guides(colour=guide_legend(nrow=1, byrow=TRUE)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c("bottom"))+theme(legend.direction="horizontal")+
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
  theme(axis.text.x = element_text(colour="black", size=21, face = "bold", hjust = -4),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "grey90"))+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing = unit(8, "lines"))+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(size = 14, family  = "Helvetica"))+
  theme(plot.background = element_rect(color = 0,
                                       size = 0),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))

override.linetype = c("solid", "dashed")

exp.compass.legend.plot =exp.compass.legend.plot + guides(colour = guide_legend(override.aes = list(linetype = override.linetype),
                                                                                 label.position = "bottom",
                                                                                 title.position = "top", title.hjust = 0.5))
exp.compass.legend.plot = exp.compass.legend.plot + scale_linetype(guide = FALSE)+
  theme(legend.title=element_text(size=22))

exp.compass.legend = get_legend(exp.compass.legend.plot+theme(legend.position = "right",
                                                              legend.title.align = 0.5))

#################### PLOT FIGURE 2 & SAVE W/ LABELS ############################
plot = cowplot::ggdraw() +
  cowplot::draw_plot(exposure.map,x=0-0.055,  y = 0.5-0.009,width = 1.1, height = 0.55) +
  cowplot::draw_plot(exp.compass,x=0,  y = 0.17-0.035,width = 1, height = 0.4) +
  cowplot::draw_plot(exposure.legend,x=-0.25,  y = -0.3-0.05,width = 1.5, height = 0.48)+
  cowplot::draw_plot(exposure.legend.size,x=0-0.55,  y = -0.29-0.05,width = 1.5, height = 0.48)+
  cowplot::draw_plot(exp.compass.legend,x=0.05,  y = -0.13-0.049,width = 1.5, height = 0.48)+
  cowplot::draw_text("A", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.1-0.05, y = 0.93)+
  cowplot::draw_text("B", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.4-0.035, y = 0.93)+
  cowplot::draw_text("C", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.7-0.02, y = 0.93)+
  cowplot::draw_text("D", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.1-0.05, y = 0.535+0.003)+
  cowplot::draw_text("E", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.4-0.035, y = 0.535+0.003)+
  cowplot::draw_text("F", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.7-0.02, y = 0.535+0.003)

ggsave(plot = plot, "[LINK TO EXPOSURE FIGURE]", width=21, height=13, dpi=300)

###############################################################################
# +++++++++++++++++++++++++ SENSITIVITY ANALYSIS  ++++++++++++++++++++++++++++#
###############################################################################
data.fig3 = data.fig2 %>%
  dplyr::filter(Scenario == c("Business-as-Usual")) ## Look at your median demand

## Define a W projected and a W_baseline for dW
data.fig3$W_proj = data.fig3$Supply_med/data.fig3$Demand_med
data.fig3$W_base = data.fig3$Supply_obs/data.fig3$Demand_obs

## For consistency with Exposure, first obtain median values
data.fig3 = data.fig3 %>%
  group_by(System, Period)%>%
  dplyr::mutate(dw = W_proj-W_base, ## there's only one W_base per sys
                W = W_proj,
                dX = Exp,
                denom =W/W_base,
                num = (dw/dX),## this is technically correct
                S = num/denom,
                W_base,
                demand_delta = ((Demand_med-Demand_obs)/Demand_obs)*100,
                supply_delta = ((Supply_med-Supply_obs)/Supply_obs)*100,
                dw_perc = (dw/W_base)*100)%>%
  ungroup() %>%
  dplyr::mutate(A_S = min(S),
                B_S = max(S))

## Now you rescale your Sensitivity values where negative should be HIGHER Sens
data.fig3 = data.fig3 %>%group_by(System, Period) %>%
  dplyr::mutate(Sen_I = 11-(1+(10-1)*((S - A_S)/(B_S-A_S))))

data.fig3$Sensitivity = data.fig3$Sen_I

###############################################################################
# +++++++++++++++++++++++++ CODE FOR FIGURE 3  ++++++++++++++++++++++++++++++#
###############################################################################
## Set up your sensitivity palette
sensitivity.pal = c('#08306b', '#08519c', '#2171b5','#4292c6','#6baed6',
                    '#9ecae1','#c6dbef','#deebf7','#f7fbff')

############################### CODE FOR FIGURE 3A-C ##########################
sensitivity.map = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), 
                                         fill = "white",
                                         color = "grey77", size = 0.25, alpha = 0.4)+
  geom_sf(data = regions, fill = "#f5f5f5", color = "grey50", alpha = 0.5) +
  geom_raster(data = rdf, aes(x=x, y=y, fill = magnitude), alpha = 0.4) +
  scale_fill_gradientn(colors = rev(pal))+
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_label(data = data.fig3,aes(x=POINT_X_lab, y=POINT_Y_lab,
                                  label = Systems), hjust = 0,vjust = -0.55, family = "Helvetica", 
             size = 7.3, fill = "#f7fcfd", label.size = 0, label.padding = unit(0, "lines"))+
  geom_point(data = data.fig3,
             aes(x=POINT_X, y=POINT_Y, fill =Sensitivity, size =storage_med_snow*100), 
             shape = 21, color = "black", alpha = 1, stroke = 1.05)+
  scale_fill_stepsn("Sensitivity",
                    colors=rev(sensitivity.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-.25,-.5,-.75,-1.0))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+  
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top")) +
  facet_wrap(~Period, nrow = 1)+theme(panel.grid= element_blank(), 
                                      legend.position = "none")+
  theme(strip.text.x = element_text(
    face = "bold",
    margin = unit(rep(25, 4), "pt")
  ))

## Make a custom scale for your fill
sensitivity.map.fill = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                              color = "grey85", size = 0.25, alpha = 0.55)+
  scale_fill_gradient2(low = "white", high = "grey20")+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig3,
             aes(x=POINT_X, y=POINT_Y, fill =Sensitivity), size = 5,
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Relative Sensitivity",
                    colors=rev(sensitivity.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("% \U0394 in ", S[Snow])),
                        range = c(15,5), breaks = c(-15,-30,-45,-60))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+  
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0))+
  guides(fill = guide_colorsteps(frame.colour = "black", ticks=TRUE,
                                 ticks.colour = "black",
                                 ticks.linewidth = 1,
                                 title.position = "top"))

sensitivity.legend = get_legend(sensitivity.map.fill+theme(legend.position = "right", 
                                                           legend.title.align = 0.5))

## Make a custom scale for your size
sensitivity.map.size = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                              color = "grey85", size = 0.25, alpha = 0.55)+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig3,
             aes(x=POINT_X, y=POINT_Y, size =storage_med_snow*100), 
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Sensitivity",
                    colors=rev(sensitivity.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("% \U0394 in ", S[Snow])),
                        range = c(15,5), breaks = c(-15,-30,-45,-60))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+  
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top",
                             label.position = "bottom")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0))

sensitivity.legend.size = get_legend(sensitivity.map.size+theme(legend.position = "right", 
                                                                legend.title.align = 0.5))

######################### CODE FOR FIGURE 3F-D #################################
data.fig3.radar = data.fig3 %>% ungroup() %>%
  dplyr::select(System, Period,  dw_perc, demand_delta) %>%
  unique() %>% group_by(System, Period) %>%
  reshape2::melt(id = c("System", "Period"))

data.fig3.radar$variable = str_replace(data.fig3.radar$variable, "dw_perc", "    \u0394 in W    ")
data.fig3.radar$variable = str_replace(data.fig3.radar$variable, "demand_delta", "\u0394 in Demand")

data.fig3.radar = data.fig3.radar %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(
                                    "Little Wood" = "LW-ID         ",
                                    "Sun" = "         Su-MT",
                                    "Shoshone" = "         Sh-WY", 
                                    "Wind" = "         Wi-WY",
                                    "Bridger" = "         Br-WY", 
                                    "Price" = "         Pr-UT",
                                    "Paonia" = "         Pa-CO", 
                                    "Costilla" = "Co-NM",
                                    "Kern" = "Ke-CA         ", 
                                    "Kaweah" = "Ka-CA         ",
                                    "Walker" = "Wa-NV         ", 
                                    "Umatilla" = "Um-OR         ",
                                    "Kittitas" = "Ki-WA         "
                              
                                  )))

data.fig3.radar$Systems = factor(data.fig3.radar$Systems, 
                                  c(
                                    "LW-ID         ",
                                    "         Su-MT",
                                    "         Sh-WY", 
                                    "         Wi-WY",
                                    "         Br-WY", 
                                    "         Pr-UT",
                                    "         Pa-CO", 
                                    "Co-NM",

                                    "Ke-CA         ", 
                                    "Ka-CA         ",
                                    "Wa-NV         ", 
                                    "Um-OR         ",
                                    "Ki-WA         "
                                    
                                  ))


df3 = data.fig3.radar %>% 
  dplyr::arrange((levels = Systems))

## Set up your radar plot for comparison
coord_radar = function (theta = "x", start = 0, direction = 1) 
{
  theta = match.arg(theta, c("x", "y"))
  r = if (theta == "x") 
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}
## Now plot
sen.compass = ggplot(df3, aes(x=Systems, y=value, color=variable,group = variable),
                     alpha = 0.5) +
  geom_hline(aes(yintercept = 1), color = "red", alpha = 0.75)+
  geom_polygon(aes(color = variable, linetype = variable), fill = NA, size =1.75) +
  theme_minimal(base_size = 30)+
  coord_polar(start = 5.7, clip = "off")+
  scale_fill_manual("",
                    values = c("grey50","black" ))+
  scale_color_manual("Components",
                     values = c("grey50","black" ))+
  facet_wrap(Period ~ ., labeller = label_wrap_gen(width=14), nrow =1) +
  ylab("")+xlab("")+
  guides(colour=guide_legend(nrow=1, byrow=TRUE)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c("bottom"))+theme(legend.direction="horizontal")+
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
  theme(axis.text.x = element_text(colour="black", size=21,  hjust = -4),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "grey80", size = 0.75))+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing = unit(9, "lines"))+
  theme(legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+
  theme(legend.position = "none")+
  theme(plot.background = element_rect(color = 0,
                                       size = 0),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) +#
  
  annotate(geom ="label", x = 1.5, alpha = 0.7,
           y = seq(-40, 20, 20), 
           label = c("-40%",  "-20%", "0%", "20%"),
           size = unit(8, "pt"), fill = "white",label.size=NA) +
  
  # remove original y-axis text / ticks & black border
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 

######################### ROSE BASED MAP OF SENSITIITY ##########################
sen.compass.legend.plot = ggplot(df3, aes(x=Systems, y=value, color=variable,group = variable),
                                 alpha = 0.5) +
  geom_line(aes(color = variable, linetype = variable), fill = NA,  size =1.5) +
  theme_minimal(base_size = 30)+
  coord_polar(start = 5.7, clip = "off")+
  scale_fill_manual("",
                    values = c("grey50","black" ))+
  scale_color_manual("Components of Sensitivity",
                     values = c("grey50","black" ))+
  facet_wrap(Period ~ ., labeller = label_wrap_gen(width=14), nrow =1) +
  ylab("")+xlab("")+
  guides(colour=guide_legend(nrow=1, byrow=TRUE)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c("bottom"))+theme(legend.direction="horizontal")+
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
  theme(axis.text.x = element_text(colour="black", size=21, face = "bold", hjust = -4),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "grey90"))+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing = unit(8, "lines"))+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(size = 14, family  = "Helvetica"))+
  theme(plot.background = element_rect(color = 0,
                                       size = 0),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))

override.linetype = c("solid", "dashed")

sen.compass.legend.plot =sen.compass.legend.plot + 
  guides(colour = guide_legend(override.aes = list(linetype = override.linetype),
                               label.position = "bottom",
                               title.position = "top", title.hjust = 0.5))

sen.compass.legend.plot = sen.compass.legend.plot + scale_linetype(guide = FALSE)+
  theme(legend.title=element_text(size=24))

sen.compass.legend = get_legend(sen.compass.legend.plot+theme(legend.position = "right",
                                                              legend.title.align = 0.5))


############################# FIGURE 3 W/ LABELS ##############################
sensitivity.plot = cowplot::ggdraw() +
  cowplot::draw_plot(sensitivity.map,x=0-0.055,  y = 0.5-0.009,width = 1.1, height = 0.55) +
  cowplot::draw_plot(sen.compass,x=0,  y = 0.17-0.035,width = 1, height = 0.4) +
  cowplot::draw_plot(sensitivity.legend,x=-0.25,  y = -0.3-0.05,width = 1.5, height = 0.48)+
  cowplot::draw_plot(sensitivity.legend.size,x=0-0.55,  y = -0.29-0.05,width = 1.5, height = 0.48)+
  cowplot::draw_plot(sen.compass.legend,x=0.05,  y = -0.13-0.049,width = 1.5, height = 0.48)+
  cowplot::draw_text("A", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.1-0.05, y = 0.93)+
  cowplot::draw_text("B", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.4-0.035, y = 0.93)+
  cowplot::draw_text("C", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.7-0.02, y = 0.93)+
  cowplot::draw_text("D", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.1-0.05, y = 0.535+0.003)+
  cowplot::draw_text("E", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.4-0.035, y = 0.535+0.003)+
  cowplot::draw_text("F", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.7-0.02, y = 0.535+0.003)
windows()
ggsave(plot = sensitivity.plot, "[LINK TO SENSITIVITY PLOT LOCATION]", width=21, height=13, dpi=300)


###############################################################################
# +++++++++++++++++++++++++ ADAPTIVE CAPACITY ANALYSIS  ++++++++++++++++++++++#
###############################################################################
## Set up your data
data = read.csv("[LINK TO PAPER RESULTS 'paper_data.csv' IN RASTER FOLDER")
data.fig4 = data

## Now get your demand and supply metrics
data.fig4 = data.fig4 %>% ungroup() %>% 
  group_by(System, Period) %>%
  dplyr::mutate(Supply_med = median(Supply),
                NaturalStorage_med = median(NaturalStorage),
                Demand_med = median(Demand[Scenario == 'Business-as-Usual']),
                Demand_BAS_median = median(Demand[Scenario == 'Business-as-Usual']),
                Demand_LWU_median = median(Demand[Scenario == 'Low Water Use']),
                Demand_mgmt = abs(Demand_LWU_median-Demand_BAS_median)/(Demand_BAS_median),
                N_B_hist = (NaturalStorage_obs/BuiltStorage),
                N_B_proj = (NaturalStorage_med/BuiltStorage),
                Storage = (N_B_hist-N_B_proj)/N_B_hist,
                RI_raw = (Supply_med-Demand_BAS_median),
                RI_mod =ifelse(RI_raw < -100*cm_2_mcm, -100*cm_2_mcm, RI_raw),
                RI_proj = (Supply_med-Demand_BAS_median)/BuiltStorage,
                RI_hist = (Supply_obs-Demand_obs)/BuiltStorage,
                RI = (RI_proj-RI_hist),
                RI_delta = ((Supply_med-Demand_BAS_median)-(Supply_obs-Demand_obs))/BuiltStorage,
                Demand_delta = (Demand_BAS_median-Demand_LWU_median)/(BuiltStorage))

## Get the minimum and maximum re-scaling values across all systems and periods
data.fig4 = data.fig4 %>%
  ungroup() %>%
  dplyr::mutate(A_Demand = min(Demand_mgmt),
                B_Demand = max(Demand_mgmt),
                A_Storage = min(Storage),
                B_Storage = max(Storage),
                A_RI = min(RI),
                B_RI = max(RI))

## Rescale the conservation and storage indicators
data.fig4 = data.fig4 %>% group_by(System, Period) %>%
  dplyr::mutate(Storage_I = 1+(10-1)*((Storage - A_Storage)/(B_Storage-A_Storage)),
                Demand_I = 1+(10-1)*((Demand_mgmt - A_Demand)/(B_Demand-A_Demand)),
                Storage_I_radar = 11-Storage_I,
                Recharge_I =1+(10-1)*((RI - A_RI)/(B_RI-A_RI)),
                Recharge_I_radar =  Recharge_I,
                RI_delta,
                Demand_delta)

## Combine into an adaptive capacity indicator via geometric mean
data.fig4$AC = (data.fig4$Demand_I * (data.fig4$Recharge_I))^(1/2)

## Get change in snowpack
data.fig4 = data.fig4 %>% ungroup()%>% group_by(System, Period) %>%
  dplyr::mutate(storage_med = median((NaturalStorage-NaturalStorage_obs)/BuiltStorage),
                storage_med_snow = median((NaturalStorage-NaturalStorage_obs)/NaturalStorage_obs)) %>%
  left_join(sites, by = c("System")) %>% unique()

data.fig4 = data.fig4 %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(
                                    "Little Wood" = "LW-ID",
                                    "Sun" = "Su-MT",
                                    "Shoshone" = "Sh-WY", 
                                    "Wind" = "Wi-WY",
                                    "Bridger" = "Br-WY", 
                                    "Price" = "Pr-UT",
                                    "Paonia" = "Pa-CO", 
                                    "Costilla" = "Co-NM",

                                    "Kern" = "Ke-CA", 
                                    "Kaweah" = "Ka-CA",
                                    "Walker" = "Wa-NV", 
                                    "Umatilla" = "Um-OR",
                                    "Kittitas" = "Ki-WA"
                                 
                                  )))
#Shift points
data.fig4$POINT_Y = ifelse(data.fig4$System == c("Shoshone"), 45,data.fig4$POINT_Y)

#Shift points
data.fig4$POINT_Y_lab = ifelse(data.fig4$System == c("Bridger"), 40.8,
                               ifelse(data.fig4$System == c("Costilla"), 34.8,
                                      ifelse(data.fig4$System == c("Kaweah"), 35.85,
                                             ifelse(data.fig4$System == c("Kern"), 33.0,
                                                    ifelse(data.fig4$System == c("Kittitas"), 47,
                                                           ifelse(data.fig4$System == c("Little Wood"), 41.5,
                                                                  ifelse(data.fig4$System == c("Paonia"),39.1,
                                                                         ifelse(data.fig4$System == c("Price"), 37.4,
                                                                                ifelse(data.fig4$System == c("Wind"),43.0,
                                                                                       ifelse(data.fig4$System == c("Sun"), 47.5,
                                                                                              ifelse(data.fig4$System == c("Umatilla"), 43.3,
                                                                                                     ifelse(data.fig4$System == c("Walker"), 38.5,
                                                                                                            45.1))))))))))))
data.fig4$POINT_X_lab = ifelse(data.fig4$System == c("Bridger"), -109.6,
                               ifelse(data.fig4$System == c("Costilla"), -112.1,
                                      ifelse(data.fig4$System == c("Kaweah"), -125.8,
                                             ifelse(data.fig4$System == c("Kern"), -124.5280,
                                                    ifelse(data.fig4$System == c("Kittitas"), -126.9,
                                                           ifelse(data.fig4$System == c("Little Wood"), -119.9,
                                                                  ifelse(data.fig4$System == c("Paonia"),-107.0,
                                                                         ifelse(data.fig4$System == c("Price"), -116.3,
                                                                                ifelse(data.fig4$System == c("Wind"),-108.3,
                                                                                       ifelse(data.fig4$System == c("Sun"), -112.2486,
                                                                                              ifelse(data.fig4$System == c("Umatilla"), -126.3,
                                                                                                     ifelse(data.fig4$System == c("Walker"), -126.3,
                                                                                                            -109))))))))))))


## Add your systems
data.fig4 = data.fig4 %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(
                                    "Little Wood" = "LW-ID",
                                    "Sun" = "Su-MT",
                                    "Shoshone" = "Sh-WY", 
                                    "Wind" = "Wi-WY",
                                    "Bridger" = "Br-WY", 
                                    "Price" = "Pr-UT",
                                    "Paonia" = "Pa-CO", 
                                    "Costilla" = "Co-NM",

                                    "Kern" = "Ke-CA", 
                                    "Kaweah" = "Ka-CA",
                                    "Walker" = "Wa-NV", 
                                    "Umatilla" = "Um-OR",
                                    "Kittitas" = "Ki-WA"
                          
                                  )))

###############################################################################
# ++++++++++++++++++++++++++++++++++++++++++ FIGURE 4  ++++++++++++++++++++++#
###############################################################################
adaptivecapacity.pal = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8',
                         '#807dba','#6a51a3','#54278f','#3f007d')

############################### CODE FOR FIGURE 4 a-c ##########################
adaptivecapacity.map = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), 
                                              fill = "white",
                                              color = "grey77", size = 0.25, alpha = 0.4)+
  geom_sf(data = regions, fill = "#f5f5f5", color = "grey50", alpha = 0.5) +
  geom_raster(data = rdf, aes(x=x, y=y, fill = magnitude), alpha = 0.4) +
  scale_fill_gradientn(colors = rev(pal))+
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_label(data = data.fig4,aes(x=POINT_X_lab, y=POINT_Y_lab,
                                  label = Systems), hjust = 0,vjust = -0.55, family = "Helvetica", 
             size = 7.3, fill = "#f7fcfd", label.size = 0, label.padding = unit(0, "lines"))+
  geom_point(data = data.fig4,
             aes(x=POINT_X, y=POINT_Y, fill =AC, size =storage_med_snow*100), 
             shape = 21, color = "black", alpha = 1, stroke = 1.05)+
  scale_fill_stepsn("Exposure",
                    colors=(adaptivecapacity.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-.25,-.5,-.75,-1.0))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+  
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top")) +
  facet_wrap(~Period, nrow = 1)+theme(panel.grid= element_blank(), 
                                      legend.position = "none")+
  theme(strip.text.x = element_text(
    face = "bold",
    margin = unit(rep(25, 4), "pt")
  ))

## Custom fill legend
adaptivecapacity.map.fill = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
 
  scale_fill_gradient2(low = "white", high = "grey20")+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig4,
             aes(x=POINT_X, y=POINT_Y, fill =AC), 
             shape = 21, color = "black", alpha = 1, stroke = 1)+
  scale_fill_stepsn("Relative Adaptive Capacity",
                    colors=(adaptivecapacity.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous("Standard Deviation in adaptivecapacity (S1-S4)", 
                        range = c(5,15))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+  
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0))+guides(fill = guide_colorsteps(frame.colour = "black", ticks=TRUE,
                                                                            ticks.colour = "black",
                                                                            ticks.linewidth = 1,
                                                                            title.position = "top"))

adaptivecapacity.legend = get_legend(adaptivecapacity.map.fill+theme(legend.position = "right", 
                                                                     legend.title.align = 0.5))
## Custom size legend
adaptivecapacity.map.size = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                                   color = "grey85", size = 0.25, alpha = 0.55)+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig4,
             aes(x=POINT_X, y=POINT_Y, size =storage_med_snow*100), 
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_size_continuous(expression(paste("% \U0394 in ", S[Snow])),
                        range = c(15,5), breaks = c(-15,-30,-45,-60))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top",
                             label.position = "bottom")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0)) 

adaptivecapacity.legend.size = get_legend(adaptivecapacity.map.size+theme(legend.position = "right", 
                                                                          legend.title.align = 0.5))

############################### CODE FOR FIGURE 4 d-f ##########################
data.fig4.radar = data.fig4 %>% 
  dplyr::select(System, Period, RI_delta, Demand_delta,
                storage_med,
                POINT_X, POINT_Y) %>%
  unique() %>% group_by(System, Period, POINT_X, POINT_Y, storage_med) %>%
  reshape2::melt(id = c("System", "Period", "storage_med", "POINT_X", "POINT_Y"))

## Information for radar plots
data.fig4.radar = data.fig4.radar %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(
                                    "Little Wood" = "LW-ID         ",
                                    "Sun" = "         Su-MT",
                                    "Shoshone" = "         Sh-WY", 
                                    "Wind" = "         Wi-WY",
                                    "Bridger" = "         Br-WY", 
                                    "Price" = "         Pr-UT",
                                    "Paonia" = "         Pa-CO", 
                                    "Costilla" = "Co-NM",

                                    "Kern" = "Ke-CA         ", 
                                    "Kaweah" = "Ka-CA         ",
                                    "Walker" = "Wa-NV         ", 
                                    "Umatilla" = "Um-OR         ",
                                    "Kittitas" = "Ki-WA         "
                                  )))

data.fig4.radar$Systems = factor(data.fig4.radar$Systems, 
                                  c(
                                    "LW-ID         ",
                                    "         Su-MT",
                                    "         Sh-WY", 
                                    "         Wi-WY",
                                    "         Br-WY", 
                                    "         Pr-UT",
                                    "         Pa-CO", 
                                    "Co-NM",

                                    "Ke-CA         ", 
                                    "Ka-CA         ",
                                    "Wa-NV         ", 
                                    "Um-OR         ",
                                    "Ki-WA         "
                                  ))

df4 = data.fig4.radar %>% 
  dplyr::arrange((levels = Systems))

## Set up your radar plot for comparison
coord_radar = function (theta = "x", start = 0, direction = 1) {
  theta = match.arg(theta, c("x", "y"))
  r = if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

## Now plot
ac.compass = ggplot(df4, aes(x=Systems, y=value*100, color=variable,group = variable),
                    alpha = 0.5) +
  geom_hline(aes(yintercept = 0), color = "red", alpha = 0.75)+
  geom_polygon(aes(color = variable, linetype = variable), fill = NA, size =1.75) +
  theme_minimal(base_size = 30)+
  coord_polar(start = 5.7, clip = "off")+
  # scale_y_continuous(lim = c(-40,20))+
  scale_fill_manual("",
                    values = c("grey50","black" ))+
  scale_color_manual("Components",
                     values = c("grey50","black" ))+
  facet_wrap(Period ~ ., labeller = label_wrap_gen(width=14), nrow =1) +
  ylab("")+xlab("")+
  guides(colour=guide_legend(nrow=1, byrow=TRUE)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c("bottom"))+theme(legend.direction="horizontal")+
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
  theme(axis.text.x = element_text(colour="black", size=21,  hjust = -4),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "grey80", size = 0.75))+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing = unit(9, "lines"))+
  theme(legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+
  theme(legend.position = "none")+
  theme(plot.background = element_rect(color = 0,
                                       size = 0),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) +#
  
  annotate(geom ="label", x = 1.5, alpha = 0.7,
           y = seq(-150, 50, 50), 
           label = c("-150%",  "-100%", "-50%", "0%", "50%"),
           size = unit(8, "pt"), fill = "white",label.size=NA) +
  
  # remove original y-axis text / ticks & black border
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 


## Custom fill legend
legend_labels = c(expression(frac("Demand Management", S[Built])), 
                   expression(frac("Supply Management", S[Built])))

ac.compass.legend.plot =  ggplot(df4, aes(x=Systems, y=value, color=variable,group = variable),
                                 alpha = 0.5) +
  geom_line(aes(color = variable, linetype = variable), fill = NA,  size =1.5) +
  theme_minimal(base_size = 30)+
  coord_polar(start = 5.7, clip = "off")+
  scale_y_continuous(lim = c(-40,10))+
  scale_fill_manual("",
                    values = c("grey50","black" ))+
  scale_color_manual("Components of Adaptive Capacity",
                     values = c("grey50","black" ),
                     labels = legend_labels)+
  facet_wrap(Period ~ ., labeller = label_wrap_gen(width=14), nrow =1) +
  ylab("")+xlab("")+
  guides(colour=guide_legend(nrow=1, byrow=TRUE)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c("bottom"))+theme(legend.direction="horizontal")+
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
  theme(axis.text.x = element_text(colour="black", size=21, face = "bold", hjust = -4),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "grey90"))+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing = unit(8, "lines"))+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(size = 14, family  = "Helvetica"))+
  theme(plot.background = element_rect(color = 0,
                                       size = 0),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))

override.linetype = c("solid", "dashed")

ac.compass.legend.plot =ac.compass.legend.plot + 
  guides(colour = 
           guide_legend(override.aes = list(linetype = override.linetype),
                        label.position = "bottom",
                        title.position = "top", title.hjust = 0.5))

ac.compass.legend.plot = ac.compass.legend.plot + scale_linetype(guide = FALSE)+
  theme(legend.title=element_text(size=24))

ac.compass.legend = get_legend(ac.compass.legend.plot+theme(legend.position = "right",
                                                            legend.title.align = 0.5))

########################## PLOT FIGURE 4 W/ LABELS #############################
ac.plot = cowplot::ggdraw() +
  cowplot::draw_plot(adaptivecapacity.map,x=0-0.055,  y = 0.5-0.009,width = 1.1, height = 0.55) +
  cowplot::draw_plot(ac.compass,x=0,  y = 0.17-0.018,width = 0.95, height = 0.38) +
  cowplot::draw_plot(adaptivecapacity.legend,x=-0.28,  y = -0.3-0.028+0.005,width = 1.5, height = 0.48)+
  cowplot::draw_plot(adaptivecapacity.legend.size,x=0-0.6,  y = -0.29-0.028+0.005,width = 1.5, height = 0.48)+
  cowplot::draw_plot(ac.compass.legend,x=0.05,  y = -0.13-0.0441+0.005,width = 1.5, height = 0.48)+
  cowplot::draw_text("A", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.1-0.05, y = 0.93)+
  cowplot::draw_text("B", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.4-0.035, y = 0.93)+
  cowplot::draw_text("C", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.7-0.02, y = 0.93)+
  cowplot::draw_text("D", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.1-0.05, y = 0.535+0.003)+
  cowplot::draw_text("E", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.4-0.035, y = 0.535+0.003)+
  cowplot::draw_text("F", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.7-0.02, y = 0.535+0.003)

ggsave(plot = ac.plot, "[LINK TO AC PLOT LOCATION]", width=21, height=13, dpi=300)


###############################################################################
# ++++++++++++++++++ VULNERABILTIY & RESILIENCE ANALYSIS  ++++++++++++++++++++#
###############################################################################

## Select exposure and sensitivity from the first two figures
es.data = data.fig2 %>% dplyr::select(System, Systems, Period, Exp, Sensitivity) %>%
  unique()

## Combine exposure, sensitivity, and adaptive capacity data into a single DF
Vulnerability.data = left_join(es.data, data.fig3, 
                               by = c("Systems", "Period", "System"))

## Calculate Vulnerability w/o adaptive capacity
Vulnerability.data$ESAC_vulnerability = (Vulnerability.data$Exp+Vulnerability.data$Sensitivity)

## Calculate Vulnerability w adaptive capacity
Vulnerability.data$ESAC_vulnerability_AC = (Vulnerability.data$Exp+Vulnerability.data$Sensitivity)-
  Vulnerability.data$AC

## Write a csv
Vulnerability.data.sub = Vulnerability.data %>% dplyr::select(System, Systems, Period,Scenario,
                                                              Exp, Sensitivity, AC,
                                                              storage_med_snow,
                                                              NaturalStorage_obs, BuiltStorage,
                                                              Demand_obs, Supply_med, Supply_obs,Elevation,
                                                              Demand_med, NaturalStorage_med, Flexible,
                                                              Capital, Flexible.hist, Captial.hist,
                                                              ESAC_vulnerability, ESAC_vulnerability_AC) %>% unique() 

## Write your csv 
table2 = Vulnerability.data.sub %>% dplyr::select("Systems", "Period", "Exp", "Sensitivity",
                                                  "AC", "ESAC_vulnerability", "ESAC_vulnerability_AC", 
                                                  "storage_med_snow", "Flexible", "Capital", "Elevation")

table2 = table2 %>%
  ungroup() %>%
  dplyr::mutate(A_vul = min(ESAC_vulnerability_AC),
                B_vul = max(ESAC_vulnerability))

## Rescale the conservation and storage indicators
table2 = table2 %>% 
  dplyr::mutate(Vulnerability = 1+(10-1)*((ESAC_vulnerability - A_vul)/(B_vul-A_vul)),
                Vulnerability_AC =1+(10-1)*((ESAC_vulnerability_AC - A_vul)/(B_vul-A_vul))) %>%
  dplyr::select(-c("ESAC_vulnerability", "ESAC_vulnerability_AC", "A_vul", "B_vul")) %>%
  dplyr::select("Systems", "Period", "Exp", "Sensitivity", "AC", 
                "Vulnerability", "Vulnerability_AC", "storage_med_snow", "Flexible", "Capital",
                "Elevation")

table2$Systems = factor(table2$Systems, 
                         levels = c("LW-ID",
                                    "Su-MT",
                                    "Sh-WY", 
                                    "Wi-WY",
                                    "Br-WY", 
                                    "Pr-UT",
                                    "Pa-CO", 
                                    "Co-NM",
                                    "Ke-CA", 
                                    "Ka-CA",
                                    "Wa-NV", 
                                    "Um-OR",
                                    "Ki-WA"))
table2 = table2 %>% 
  dplyr::arrange((levels = Systems)) %>% unique()


write.csv(table2,"[LINK TO VULNERABILITY RESULTS]")

###############################################################################
# ++++++++++++++++++ CORRELATION PLOTS FOR SI FIGURES ++++++++++++++++++++++++#
###############################################################################
table2 = table2 %>%
  mutate(color =  plyr::revalue(Systems,
                                c(
                                  "LW-ID" = "dodgerblue4", 
                                  "Ki-WA" = "cornflowerblue",
                                  "Um-OR" = "deepskyblue",
                                  "Br-WY" = "orangered2",
                                  "Pa-CO" = "coral3", 
                                  "Pr-UT" = "lightcoral",
                                  "Co-NM" = "peachpuff", 
                                  "Ka-CA" = "darkgoldenrod3",
                                  "Ke-CA" = "lightgoldenrod1",
                                  "Su-MT" = "olivedrab3", 
                                  "Sh-WY" = "lightseagreen", 
                                  "Wi-WY" = "darkturquoise",
                                  "Wa-NV" = "slategray2" 
                                )))

myColorScale = unique(as.vector(table2$color))
names(myColorScale) = unique(as.vector(table2$Systems))
colScale = scale_fill_manual(name = "Systems",values = myColorScale)

## EXPOSURE VS SENSITIVITY BUT CAN BE SWITCHED TO EXAMINE OTHER VARIABLES
ggplot(data = table2)+geom_point(aes(x= Exp, y = Sensitivity, fill = Systems), size = 5, 
                                 shape = 21)+
  theme_minimal(base_size = 22)+facet_wrap(~Period)+colScale

cor(table2$Exp, table2$Sensitivity)

###############################################################################
# ++++++++++++++++++++++++++++++++ CODE FOR FIGURE 5  ++++++++++++++++++++++++#
###############################################################################

data.fig5 = Vulnerability.data %>% ungroup()%>% group_by(System, Period) %>%
  dplyr::select(System, Period, POINT_X, POINT_Y, ESAC_vulnerability,
                ESAC_vulnerability_AC,storage_med, storage_med_snow,
                Flexible.hist) %>%
  reshape2::melt(id = c("System", "Period", "POINT_X", "POINT_Y","storage_med", 
                        "storage_med_snow", "Flexible.hist")) %>%
  ungroup() 

## Now melt for figure
data.fig5$variable_rename = ifelse(data.fig5$variable == c("ESAC_vulnerability_AC"), "Adaptation",
                                   "No Adaptation")

data.fig5 = transform(data.fig5,
                      Name=factor(variable_rename,levels=c("No Adaptation", 
                                                           "Adaptation")))

## Get the minimum and maximum re-scaling values across all systems and periods
data.fig5 = data.fig5 %>%
  ungroup() %>%
  dplyr::mutate(A_vul = min(value),
                B_vul = max(value))

## Rescale the conservation and storage indicators
data.fig5 = data.fig5 %>% 
  dplyr::mutate(Vul_I = 1+(10-1)*((value - A_vul)/(B_vul-A_vul)))

## Change the point labels
#Shift points
data.fig3$POINT_Y = ifelse(data.fig3$System == c("Shoshone"), 45,data.fig3$POINT_Y)

#Shift points
data.fig3$POINT_Y_lab = ifelse(data.fig3$System == c("Bridger"), 40.7,
                               ifelse(data.fig3$System == c("Costilla"), 34.5,
                                      ifelse(data.fig3$System == c("Kaweah"), 35.85,
                                             ifelse(data.fig3$System == c("Kern"), 33.0,
                                                    ifelse(data.fig3$System == c("Kittitas"), 47,
                                                           ifelse(data.fig3$System == c("Little Wood"), 41.3,
                                                                  ifelse(data.fig3$System == c("Paonia"),39.1,
                                                                         ifelse(data.fig3$System == c("Price"), 37.4,
                                                                                ifelse(data.fig3$System == c("Wind"),43.0,
                                                                                       ifelse(data.fig3$System == c("Sun"), 47.5,
                                                                                              ifelse(data.fig3$System == c("Umatilla"), 42.8,
                                                                                                     ifelse(data.fig3$System == c("Walker"), 38.5,
                                                                                                            45.1))))))))))))
data.fig3$POINT_X_lab = ifelse(data.fig3$System == c("Bridger"), -109.4,
                               ifelse(data.fig3$System == c("Costilla"), -113.5,
                                      ifelse(data.fig3$System == c("Kaweah"), -127.0,
                                             ifelse(data.fig3$System == c("Kern"), -126.1,
                                                    ifelse(data.fig3$System == c("Kittitas"), -127.5,
                                                           ifelse(data.fig3$System == c("Little Wood"), -121.3,
                                                                  ifelse(data.fig3$System == c("Paonia"),-107.0,
                                                                         ifelse(data.fig3$System == c("Price"), -117.7,
                                                                                ifelse(data.fig3$System == c("Wind"),-108.3,
                                                                                       ifelse(data.fig3$System == c("Sun"), -112.2486,
                                                                                              ifelse(data.fig3$System == c("Umatilla"), -127.5,
                                                                                                     ifelse(data.fig3$System == c("Walker"), -127.5,
                                                                                                            -109))))))))))))

unique(data.fig3$System)

############################# Figure 5 a-f #####################################
vulnerability.map = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), 
                                           fill = "white",
                                           color = "grey77", size = 0.25, alpha = 0.4)+
  # geom_sf(data = regions,aes(fill = NAME, color = NAME), alpha = 0.2) +
  geom_sf(data = regions, fill = "#f5f5f5", color = "grey50", alpha = 0.5) +
  geom_raster(data = rdf, aes(x=x, y=y, fill = magnitude), alpha = 0.4) +
  scale_fill_gradientn(colors = rev(pal))+
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_label(data = data.fig3,aes(x=POINT_X_lab, y=POINT_Y_lab,
                                  label = Systems), hjust = 0,vjust = -0.55, family = "Helvetica", 
             size = 6.9, fill = "#f7fcfd", label.size = 0, label.padding = unit(0, "lines"))+
  geom_point(data = data.fig5,
             aes(x=POINT_X, y=POINT_Y, fill =Vul_I, size =storage_med_snow*100), 
             shape = 21, color = "black", alpha = 1, stroke = 1.05)+
  scale_fill_stepsn("Relative Resilience or Vulnerability",
                    colors=c("#3288bd", "#74add1", #blues
                             "#66c2a5", "#a6d96a", #greens
                             "#e6f598", "#ffffbf", #yellows
                             "#fee08b",  #oranges
                             "#f46d43", "#d53e4f"), #darkreds
                    n.breaks=13, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-.25,-.5,-.75,-1.0))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica", face = "bold"))+
  facet_grid(Name~Period)+guides(fill = guide_colorsteps(frame.colour = "black", ticks=TRUE,
                                                         ticks.colour = "black",
                                                         ticks.linewidth = 1)) +
  theme(panel.grid= element_blank(),
        legend.position = "none")+ 
  theme(panel.spacing.y = unit(0.2, "lines"))+
  theme(strip.text.x = element_text(
    face = "bold",
    margin = unit(rep(15, 4), "pt")
  ))

## Construct your own vulnerabiltiy map fill
vulnerability.map.fill = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                                color = "grey85", size = 0.25, alpha = 0.55)+
  scale_fill_gradient2(low = "white", high = "grey20")+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig5,
             aes(x=POINT_X, y=POINT_Y, fill =value),
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Relative Vulnerability\nor Resilience\n",
                    colors=c("#3288bd", "#74add1", #blues
                             "#66c2a5", "#a6d96a", #greens
                             "#e6f598", "#ffffbf", #yellows
                             "#fee08b",  #oranges
                             "#f46d43", "#d53e4f"), #darkreds
                    n.breaks=13, limits = c(1,10), show.limits=T)+
  scale_size_continuous("Standard Deviation in adaptivecapacity (S1-S4)", 
                        range = c(5,15))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(4, "cm"),
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"),
        legend.direction = "vertical")+  
  theme(legend.title.align = 0.5,
        legend.direction = "vertical",
        legend.box.just = "center")+  
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0))+
  guides(fill = guide_colorsteps(frame.colour = "black", ticks=TRUE,
                                 ticks.colour = "black",
                                 ticks.linewidth = 1,
                                 title.position = "top")) + theme(legend.title.align=0.5)

vulnerability.legend = get_legend(vulnerability.map.fill)

## Custom map size scale
vulnerability.map.size = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                                color = "grey85", size = 0.25, alpha = 0.55)+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig5,
             aes(x=POINT_X, y=POINT_Y, size = storage_med_snow*100),
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Vulnerability",
                    colors=c("#5e4fa2", "#3288bd", "#74add1", #Low
                             "#e6f598", "#fee08b", "#fdae61", #moderate
                             "#d53e4f", "#b2182b", "#67001f"),#high
                    n.breaks=9, limits = c(-8,20), show.limits=T)+
  scale_size_continuous(expression(paste("% \U0394 in ", S[Snow])),
                        range = c(15,5), breaks = c(-15,-30,-45,-60))+theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"),
        legend.direction = "vertical")+  
  theme(legend.title.align = 0.5,
        legend.direction = "vertical",
        legend.box.just = "center")+
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0))+ theme(legend.title.align=0.5)

vulnerability.legend.size = get_legend(vulnerability.map.size+theme(legend.position = "right", 
                                                                    legend.title.align = 0.5))

############################# Figure 5 g-i #####################################
data.fig5.radar = data.fig5 %>% 
  dplyr::select(System, Period,  Vul_I, variable_rename) %>%
  unique() %>% 
  tidyr::spread(key = variable_rename, value = Vul_I)

names(data.fig5.radar) = c("System", "Period", "Adapation", "NonAdaptation")

data.fig5.radar$Difference = abs(data.fig5.radar$Adapation-data.fig5.radar$NonAdaptation)
data.fig5.radar$Difference_per = (abs(data.fig5.radar$Adapation-data.fig5.radar$NonAdaptation)/data.fig5.radar$NonAdaptation)*100

## Radar plot
data.fig5.radar = data.fig5.radar %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(
                                    "Little Wood" = "LW-ID         ",
                                    "Sun" = "         Su-MT",
                                    "Shoshone" = "         Sh-WY", 
                                    "Wind" = "         Wi-WY",
                                    "Bridger" = "         Br-WY", 
                                    "Price" = "         Pr-UT",
                                    "Paonia" = "         Pa-CO", 
                                    "Costilla" = "Co-NM",

                                    "Kern" = "Ke-CA         ", 
                                    "Kaweah" = "Ka-CA         ",
                                    "Walker" = "Wa-NV         ", 
                                    "Umatilla" = "Um-OR         ",
                                    "Kittitas" = "Ki-WA         "
                                  )))

data.fig5.radar$Systems = factor(data.fig5.radar$Systems, 
                                  c(
                                    "LW-ID         ",
                                    "         Su-MT",
                                    "         Sh-WY", 
                                    "         Wi-WY",
                                    "         Br-WY", 
                                    "         Pr-UT",
                                    "         Pa-CO", 
                                    "Co-NM",

                                    "Ke-CA         ", 
                                    "Ka-CA         ",
                                    "Wa-NV         ", 
                                    "Um-OR         ",
                                    "Ki-WA         "
                                    #Southwest-dark yellows
                                    #"Souris-Red-Rainy " = "lightsalmon1",
                                    # West north central-light yellows
                                    #East north central-light greens
                                    # Central-medium greens
                                  ))

df5 = data.fig5.radar %>% 
  dplyr::arrange((levels = Systems))


Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")

## Set up your radar plot for comparison
coord_radar = function (theta = "x", start = 0, direction = 1) {
  theta = match.arg(theta, c("x", "y"))
  r = if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

## Set up your data and factor as a rename
df5 = df5 %>% dplyr::select(Systems, Period, NonAdaptation, Difference) 
df5 = data.fig5.radar %>% 
  dplyr::arrange((levels = Systems))

## Now plot
vulnerability.rose =ggplot(df5) +
  geom_col(aes(x=Systems, 
               y=Difference, fill = NonAdaptation),
           position = "dodge2", show.legend = TRUE, alpha = 0.9, color = "black")+
  coord_polar(start = 5.7, clip = "off")+
  facet_wrap(Period ~ ., labeller = label_wrap_gen(width=14), nrow =1) +
  ylab("")+xlab("")+  
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c("bottom"))+theme(legend.direction="horizontal")+
  scale_fill_stepsn("Relative Resilience or Vulnerability",
                    colors=c("#3288bd", "#74add1", #blues
                             "#66c2a5", "#a6d96a", #greens
                             "#e6f598", "#ffffbf", #yellows
                             "#fee08b",  #oranges
                             "#f46d43", "#d53e4f"), #darkreds
                    n.breaks=13, limits = c(1,10), show.limits=T)+
  theme_minimal(base_size = 30)+
  ylab("")+xlab("")+
  guides(colour=guide_legend(nrow=1, byrow=TRUE)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c("bottom"))+theme(legend.direction="horizontal")+
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
  theme(axis.text.x = element_text(colour="black", size=21,  hjust = -4),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "grey80", size = 0.75))+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing = unit(9, "lines"))+
  theme(legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+
  theme(legend.position = "none")+
  theme(plot.background = element_rect(color = 0,
                                       size = 0),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) +#
  
  annotate(geom ="label", x = 1.6, alpha = 0.5,
           y = seq(0, 4, 1), 
           label = c("0",  "1", "2", "3", "4"),
           size = unit(8, "pt"), fill = "white",label.size=NA) +
  
  # remove original y-axis text / ticks & black border
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 


################### PLOT FIGURE 5 W LABELS ####################################
vulnerability = cowplot::ggdraw() +
  cowplot::draw_plot(vulnerability.map,x=0-0.15,  y = 0.275-0.005,width = 1.1, height = 0.76) + 
  cowplot::draw_plot(vulnerability.rose,x=0-0.115,  y = 0.0-0.005,width = 1.0, height = 0.28)  +
  cowplot::draw_plot(vulnerability.legend.size,x=0.12,  y = -0.25,width = 1.5, height = 0.48)+
  cowplot::draw_plot(vulnerability.legend,x=0.12,  y = 0.5,width = 1.5, height = 0.48)+
  cowplot::draw_text("A", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.04, y = 0.95)+
  cowplot::draw_text("B", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.3-0.01, y = 0.95)+
  cowplot::draw_text("C", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.55-0.01, y = 0.95)+
  cowplot::draw_text("D", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.04, y = 0.63)+
  cowplot::draw_text("E", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.3-0.01, y = 0.63)+
  cowplot::draw_text("F", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.55-0.01, y = 0.63)+
  cowplot::draw_text("G", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.04, y = 0.30)+
  cowplot::draw_text("H", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.3-0.01, y = 0.30)+
  cowplot::draw_text("I", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.55-0.01, y = 0.30)+
  cowplot::draw_text("Difference", 
                     size = 30, angle = -90,
                     family = "Helvetica", x= 0.77, y = 0.15)

ggsave(plot = vulnerability, "[LINK TO FIGURE 5]", width=21, height=13, dpi=300)

