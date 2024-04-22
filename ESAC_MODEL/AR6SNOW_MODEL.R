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

sites = read.csv("[LINK TO SITE COORDS WITHIN LU TABLE FOLDER")[,-c(1:2,4)]
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
# +++++++++++++++++++++++++ HAZARD ANALYSIS  ++++++++++++++++++++++++++++++#
###############################################################################
## read in your processed data
data = read.csv("D:/Projects/Project_Reservoir/Data/Spreadsheets/System Results/paper_data.csv")

data = read.csv("[READ IN PROCESSED DATA OBTAINED FROM ACCOMPANYING SCRIPTS 
                OR PROVIDED AS 'paper_data.csv' in RASTERS FOLDER")

## First set up your hazard DF
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

## Here you get the hazard metric either as a geometric mean or as just an average
data.fig1$Haz = ((11-data.fig1$Mag_I) * (data.fig1$Timing_I))^(1/2)## Geometric mean

## Here you get the median change in storage across the GCM-RCP pairs
data.fig1 = data.fig1 %>% ungroup()%>% group_by(System, Period) %>%
  dplyr::mutate(storage_med = median((NaturalStorage-NaturalStorage_obs)/BuiltStorage),
                storage_med_snow = median((NaturalStorage-NaturalStorage_obs)/NaturalStorage_obs)) %>%
  left_join(sites, by = c("System")) %>% unique()



###############################################################################
# +++++++++++++++++++++++++ CODE FOR FIGURE 2  +++++++++++++++++++++++++++++++#
###############################################################################

######################### METRIC BASED MAP OF HAZARD ##########################
## set up your custom palette
hazard.pal = c('#67000d','#a50f15','#cb181d','#ef3b2c',
               '#fb6a4a','#fc9272','#fcbba1','#fee0d2', '#fff5f0')

## set up your labels for plotting
data.fig1$System_fac = as.factor(as.character(data.fig1$System))
data.fig1$System = ifelse(data.fig1$System == c("Little\nWood"), "Little Wood",
                          data.fig1$System)

#Shift points
data.fig1$POINT_Y = ifelse(data.fig1$System == c("Shoshone"), 45,data.fig1$POINT_Y)

#Shift points
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
                                  c(# Northwest-blue
                                    "Little Wood" = "LW-ID",
                                    "Sun" = "Su-MT",
                                    "Shoshone" = "Sh-WY", 
                                    "Wind" = "Wi-WY",
                                    "Bridger" = "Br-WY", 
                                    "Price" = "Pr-UT",
                                    "Paonia" = "Pa-CO", 
                                    "Costilla" = "Co-NM",
                                    # West-orangered
                                    "Kern" = "Ke-CA", 
                                    "Kaweah" = "Ka-CA",
                                    "Walker" = "Wa-NV", 
                                    "Umatilla" = "Um-OR",
                                    "Kittitas" = "Ki-WA"
                                    
                                  )))
windows()

hazard.map = ggplot() +
  geom_polygon(data = wus, aes(long, lat, group = group),
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
             aes(x=POINT_X, y=POINT_Y, fill =Haz, size =storage_med*100), 
             shape = 21, color = "black", alpha = 1, stroke = 1.05)+
  scale_fill_stepsn("Hazard",
                    colors=rev(hazard.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-.25,-.5,-.75,-1.0))+
  theme_map(base_size = 30)+
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

hazard.map.fill = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                         color = "grey85", size = 0.25, alpha = 0.55)+
  scale_fill_gradient2(low = "white", high = "grey20")+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig1,
             aes(x=POINT_X, y=POINT_Y, fill =Exp), size = 5,
             shape = 21, color = "black", alpha = 1, stroke = 1)+
  scale_fill_stepsn("Hazard Index",
                    colors=rev(hazard.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous("Standard Deviation in DOQ50 Timing [Days]",
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

hazard.legend = get_legend(hazard.map.fill+theme(legend.position = "right",
                                                 legend.title.align = 0.5))

hazard.map.size = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                         color = "grey85", size = 0.25, alpha = 0.55)+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  
  geom_point(data = data.fig1,
             aes(x=POINT_X, y=POINT_Y, size =storage_med*100),
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Hazard",
                    colors=rev(hazard.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
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

hazard.legend.size = get_legend(hazard.map.size+theme(legend.position = "right",
                                                      legend.title.align = 0.5))

######################### ROSE BASED MAP OF HAZARD ##########################
data.fig1.radar = data.fig1 %>% ungroup() %>%
  dplyr::select(System, Period,  timing_delta, magnitude_delta) %>%
  unique() %>% group_by(System, Period) %>%
  reshape2::melt(id = c("System", "Period"))


data.fig1.radar$variable = str_replace(data.fig1.radar$variable, "timing_delta", " \u0394 in Timing ")
data.fig1.radar$variable = str_replace(data.fig1.radar$variable, "magnitude_delta", "\u0394 in Magnitude")


data.fig1.radar = data.fig1.radar %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(# Northwest-blue
                                    "Little Wood" = "LW-ID         ",
                                    "Sun" = "         Su-MT",
                                    "Shoshone" = "         Sh-WY", 
                                    "Wind" = "         Wi-WY",
                                    "Bridger" = "         Br-WY", 
                                    "Price" = "         Pr-UT",
                                    "Paonia" = "         Pa-CO", # maybe tan 3
                                    "Costilla" = "Co-NM",
                                    # West-orangered
                                    "Kern" = "Ke-CA         ", 
                                    "Kaweah" = "Ka-CA         ",
                                    "Walker" = "Wa-NV         ", 
                                    "Umatilla" = "Um-OR         ",
                                    "Kittitas" = "Ki-WA         "
                                  )))

data.fig1.radar$Systems <- factor(data.fig1.radar$Systems, 
                                  c(# Northwest-blue
                                    "LW-ID         ",
                                    "         Su-MT",
                                    "         Sh-WY", 
                                    "         Wi-WY",
                                    "         Br-WY", 
                                    "         Pr-UT",
                                    "         Pa-CO", # maybe tan 3
                                    "Co-NM",
                                    # West-orangered
                                    "Ke-CA         ", 
                                    "Ka-CA         ",
                                    "Wa-NV         ", 
                                    "Um-OR         ",
                                    "Ki-WA         "
                                  ))

df1 <- data.fig1.radar %>% 
  dplyr::arrange((levels = Systems))
## Set up your radar plot for comparison
coord_radar <- function (theta = "x", start = 0, direction = 1) 
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") 
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

## Now plot
windows()
haz.compass = ggplot(df1, aes(x=Systems, y=value, color=variable,group = variable),
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

######################### ROSE BASED MAP OF EXPOSURE ##########################
windows()
haz.compass.legend.plot = ggplot(df1, aes(x=Systems, y=value, 
                                          color=variable,group = variable),
                                 alpha = 0.5) +
  geom_line(aes(color = variable, linetype = variable), fill = NA,  size =1.5) +
  theme_minimal(base_size = 30)+
  coord_polar(start = 5.7, clip = "off")+
  scale_y_continuous(lim = c(-40,10))+
  scale_fill_manual("",
                    values = c("grey50","black" ))+
  scale_color_manual("Components of Hazard Index",
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

override.linetype <- c("solid", "dashed")

haz.compass.legend.plot <-haz.compass.legend.plot + 
  guides(colour = guide_legend(override.aes = list(linetype = override.linetype),
                               label.position = "bottom",
                               title.position = "top", title.hjust = 0.5))
haz.compass.legend.plot <- haz.compass.legend.plot + scale_linetype(guide = FALSE)+
  theme(legend.title=element_text(size=22))

haz.compass.legend = get_legend(haz.compass.legend.plot+theme(legend.position = "right",
                                                              legend.title.align = 0.5))


windows()
plot = cowplot::ggdraw() +
  cowplot::draw_plot(hazard.map,x=0-0.055,  y = 0.5-0.009,width = 1.1, height = 0.55) +
  cowplot::draw_plot(haz.compass,x=0,  y = 0.17-0.035,width = 1, height = 0.4) +
  cowplot::draw_plot(hazard.legend,x=-0.22,  y = -0.3-0.05,width = 1.5, height = 0.48)+
  cowplot::draw_plot(hazard.legend.size,x=0-0.55,  y = -0.29-0.05,width = 1.5, height = 0.48)+
  cowplot::draw_plot(haz.compass.legend,x=0.08,  y = -0.13-0.049,width = 1.5, height = 0.48)+
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
                     family = "Helvetica", fontface = "bold", x= 0.7-0.02, y = 0.535+0.003)+
  cowplot::draw_text("\U2190 Lower", 
                     size = 22, angle = 0,
                     family = "Helvetica", fontface = "italic", x=0.44-0.02, y = 0.095)+
  cowplot::draw_text("Higher \U2192", 
                     size = 22, angle = 0,
                     family = "Helvetica", fontface = "italic", x= 0.64-0.0, y = 0.095)

ggsave(plot = plot, "[LINK TO STORED HAZARD MAP", width=21, height=13, dpi=300)

###############################################################################
# +++++++++++++++++++++++++ SENSITIVITY ANALYSIS  ++++++++++++++++++++++++++++#
###############################################################################
data.fig3 = data.fig1 %>%
  dplyr::filter(Scenario == c("Business-as-Usual")) ## Look at your median demand

## Define a W projected and a W_baseline for dW
data.fig3$W_proj = data.fig3$Supply_med/data.fig3$Demand_med
data.fig3$W_base = data.fig3$Supply_obs/data.fig3$Demand_obs

## For consistency with Exposure, first obtain median values
data.fig3 = data.fig3 %>%
  group_by(System, Period)%>%
  dplyr::mutate(dw = W_proj-W_base, ## there's only one W_base per sys
                W = W_proj,
                dX = 1,
                denom =W/W_base,
                num = abs(dw),## this is technically correct
                S = dw/W_base) %>%
  ungroup() %>%
  dplyr::mutate(A_S = min(S),
                B_S = max(S))

## Now you rescale your vulnerability values
data.fig3 = data.fig3 %>%group_by(System, Period) %>%
  dplyr::mutate(Vul_I = 11-(1+(10-1)*((S - A_S)/(B_S-A_S))))

data.fig3$Vulernability = data.fig3$Vul_I

###############################################################################
# +++++++++++++++++++++++++ CODE FOR FIGURE 3  ++++++++++++++++++++++++++++++#
###############################################################################
## Set up your vulnerability palette
vulnerability.pal = c('#08306b', '#08519c', '#2171b5','#4292c6','#6baed6',
                    '#9ecae1','#c6dbef','#deebf7','#f7fbff')

############################### CODE FOR FIGURE 3A-C ##########################
vulnerability.map = ggplot() +
  geom_polygon(data = wus, aes(long, lat, group = group),
                                         fill = "white",
                                         color = "grey77", size = 0.25, alpha = 0.4)+
  geom_sf(data = regions, fill = "#f5f5f5", color = "grey50", alpha = 0.5) +
  geom_raster(data = rdf, aes(x=x, y=y, fill = magnitude), alpha = 0.4) +
  scale_fill_gradientn(colors = rev(pal))+
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_label(data = data.fig3,aes(x=POINT_X_lab, y=POINT_Y_lab,
                                  label = Systems), hjust = 0,vjust = -0.55, 
             family = "Helvetica", 
             size = 7.3, fill = "#f7fcfd", label.size = 0, 
             label.padding = unit(0, "lines"))+
  geom_point(data = data.fig3,
             aes(x=POINT_X, y=POINT_Y, fill =Vulernability, size =storage_med_snow*100), 
             shape = 21, color = "black", alpha = 1, stroke = 1.05)+
  scale_fill_stepsn("Vulernability",
                    colors=rev(vulnerability.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-.25,-.5,-.75,-1.0))+
  theme_map(base_size = 30)+
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

vulnerability.map = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group),
                                         fill = "white",
                                         color = "grey77", size = 0.25, alpha = 0.4)+
  geom_sf(data = regions, fill = "#f5f5f5", color = "grey50", alpha = 0.5) +
  geom_raster(data = rdf, aes(x=x, y=y, fill = magnitude), alpha = 0.4) +
  scale_fill_gradientn(colors = rev(pal))+
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_label(data = data.fig2,aes(x=POINT_X_lab, y=POINT_Y_lab,
                                  label = Systems), hjust = 0,vjust = -0.55,
             family = "Helvetica",
             size = 7.3, fill = "#f7fcfd", label.size = 0, label.padding = unit(0, "lines"))+
  geom_point(data = data.fig2,
             aes(x=POINT_X, y=POINT_Y, fill =Vulnerability, size =storage_med*100),
             shape = 21, color = "black", alpha = 1, stroke = 1.05)+
  scale_fill_stepsn("Exposure",
                    colors=rev(vulnerability.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-.25,-.5,-.75,-1.0))+
  theme_map(base_size = 30)+
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

vulnerability.map.fill = ggplot() +geom_polygon(data = wus,
                                                aes(long, lat, group = group), fill = "white",
                                              color = "grey85", size = 0.25, alpha = 0.55)+
  scale_fill_gradient2(low = "white", high = "grey20")+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig2,
             aes(x=POINT_X, y=POINT_Y, fill =Vulnerability), size = 5,
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Vulnerability Index",
                    colors=rev(vulnerability.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("% \U0394 in ", S[Snow])),
                        range = c(15,5), breaks = c(-15,-30,-45,-60))+
  theme_map(base_size = 30)+
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

vulnerability.legend = get_legend(vulnerability.map.fill+theme(legend.position = "right",
                                                           legend.title.align = 0.5))

vulnerability.map.size = ggplot() +geom_polygon(data = wus, 
                                                aes(long, lat, group = group), fill = "white",
                                              color = "grey85", size = 0.25, alpha = 0.55)+

  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig2,
             aes(x=POINT_X, y=POINT_Y, size =storage_med*100),
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Vulnerability",
                    colors=rev(vulnerability.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-15,-30,-45,-60))+
  theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+  
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top",
                             label.position = "bottom")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0))

vulnerability.legend.size = get_legend(vulnerability.map.size+theme(legend.position = "right",
                                                                legend.title.align = 0.5))

######################### ROSE BASED MAP OF Vulnerability ##########################
data.fig2.radar = data.fig2 %>% ungroup() %>
  dplyr::select(System, Period,  supply_delta, demand_delta) %>%
  unique() %>% group_by(System, Period) %>%
  reshape2::melt(id = c("System", "Period"))
data.fig2.radar$variable = as.character(as.factor(data.fig2.radar$variable))
data.fig2.radar <- data.fig2.radar %>%
  dplyr::mutate(variable_1 = str_replace(variable, "supply_delta", "    \u0394 in Supply    "),
                variable_1 = str_replace(variable_1, "demand_delta", "\u0394 in Demand"))
data.fig2.radar$variable_2 = str_replace(data.fig2.radar$variable, "supply_delta", "A")
data.fig2.radar$variable_2 = str_replace(data.fig2.radar$variable, "demand_delta", "B")

data.fig2.radar = data.fig2.radar %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(# Northwest-blue
                                    "Little Wood" = "LW-ID         ",
                                    "Sun" = "         Su-MT",
                                    "Shoshone" = "         Sh-WY",
                                    "Wind" = "         Wi-WY",
                                    "Bridger" = "         Br-WY",
                                    "Price" = "         Pr-UT",
                                    "Paonia" = "         Pa-CO", # maybe tan 3
                                    "Costilla" = "Co-NM",
                                    # West-orangered
                                    "Kern" = "Ke-CA         ",
                                    "Kaweah" = "Ka-CA         ",
                                    "Walker" = "Wa-NV         ",
                                    "Umatilla" = "Um-OR         ",
                                    "Kittitas" = "Ki-WA         "
                                  )))

data.fig2.radar$Systems <- factor(data.fig2.radar$Systems,
                                  c(# Northwest-blue
                                    "LW-ID         ",
                                    "         Su-MT",
                                    "         Sh-WY",
                                    "         Wi-WY",
                                    "         Br-WY",
                                    "         Pr-UT",
                                    "         Pa-CO", # maybe tan 3
                                    "Co-NM",
                                    # West-orangered
                                    "Ke-CA         ",
                                    "Ka-CA         ",
                                    "Wa-NV         ",
                                    "Um-OR         ",
                                    "Ki-WA         "
                                  ))

df2.base <- data.fig2.radar %>%
  dplyr::arrange((levels = Systems))%>% dplyr::filter(variable == c("W_base"))

df2 <- data.fig2.radar %>%
  dplyr::arrange((levels = Systems))

## Set up your radar plot for comparison
coord_radar <- function (theta = "x", start = 0, direction = 1)
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}
## Now plot
windows()
vul.compass = ggplot(df2, aes(x=Systems, y=value*100, color=variable_2,group = variable_2),
                     alpha = 0.5) +
  geom_hline(aes(yintercept = 1), color = "red", alpha = 0.75)+
  geom_polygon(aes(color = variable_2, linetype = variable_2), fill = NA, size =1.75) +
  theme_minimal(base_size = 30)+
  # scale_linetype_manual("", values=c("dashed", "solid")) +  # Set linetypes
  coord_polar(start = 5.7, clip = "off")+
  # scale_y_continuous(lim = c(-40,20))+
  scale_fill_manual("",
                    values = c("grey50","black" ))+
  scale_color_manual("Components",
                     values = c("grey50","black"  ))+
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
# override.linetype <- c("solid", "dashed")

######################### ROSE BASED MAP OF EXPOSURE ##########################
vul.compass.legend.plot = ggplot(df2, aes(x=Systems, y=value*100, 
                                          color=variable_1,group = variable_1),
                                 alpha = 0.5) +
  geom_line(aes(color = variable_1, linetype = variable_1), fill = NA,  size =1.5) +
  theme_minimal(base_size = 30)+
  coord_polar(start = 5.7, clip = "off")+
  scale_fill_manual("",
                    values = c("black","brey50" ))+
  scale_color_manual("Components of Vulnerability Index",
                     values = c("black","grey50" ))+
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

override.linetype <- c("dashed","solid")

vul.compass.legend.plot <-vul.compass.legend.plot + 
  guides(colour = guide_legend(override.aes = list(linetype = override.linetype),
                               label.position = "bottom",
                               title.position = "top", title.hjust = 0.5))

vul.compass.legend.plot <- vul.compass.legend.plot + scale_linetype(guide = FALSE)+
  theme(legend.title=element_text(size=24))

vul.compass.legend = get_legend(vul.compass.legend.plot+theme(legend.position = "right",
                                                              legend.title.align = 0.5))



windows()
vulnerability.plot = cowplot::ggdraw() +
  cowplot::draw_plot(vulnerability.map,x=0-0.055, 
                     y = 0.5-0.009,width = 1.1, height = 0.55) +
  cowplot::draw_plot(vul.compass,x=0,  
                     y = 0.17-0.035,width = 1, height = 0.4) +
  cowplot::draw_plot(vulnerability.legend,x=-0.22-0.03,  
                     y = -0.3-0.05,width = 1.5, height = 0.48)+
  cowplot::draw_plot(vulnerability.legend.size,x=0-0.57,  
                     y = -0.29-0.05,width = 1.5, height = 0.48)+
  cowplot::draw_plot(vul.compass.legend,x=0.08, 
                     y = -0.13-0.049,width = 1.5, height = 0.48)+
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
                     family = "Helvetica", fontface = "bold", x= 0.7-0.02, y = 0.535+0.003)+
  cowplot::draw_text("\U2190 Lower", 
                     size = 22, angle = 0,
                     family = "Helvetica", fontface = "italic", x=0.44-0.02-0.03, y = 0.095)+
  cowplot::draw_text("Higher \U2192", 
                     size = 22, angle = 0,
                     family = "Helvetica", fontface = "italic", x= 0.64-0.0-0.03, y = 0.095)
windows()
ggsave(plot = vulnerability.plot, "[LINK TO STORED FILE]", width=21, height=13, dpi=300)

###############################################################################
# +++++++++++++++++++++++++ ADAPTIVE CAPACITY ANALYSIS  ++++++++++++++++++++++#
###############################################################################
## Set up your data
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
# +++++++++++++++++++++++++ CODE FOR FIGURE 4  ++++++++++++++++++++++++++++++#
###############################################################################
## Set up your adaptive capacity palette
adaptivecapacity.pal = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8',
                         '#807dba','#6a51a3','#54278f','#3f007d')

############################### CODE FOR FIGURE 4 a-c ##########################
adaptivecapacity.map = ggplot() +
  geom_polygon(data = wus, aes(long, lat, group = group),
                                              fill = "white",
                                              color = "grey77", size = 0.25, alpha = 0.4)+
  geom_sf(data = regions, fill = "#f5f5f5", color = "grey50", alpha = 0.5) +
  geom_raster(data = rdf, aes(x=x, y=y, fill = magnitude), alpha = 0.4) +
  scale_fill_gradientn(colors = rev(pal))+
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_label(data = data.fig4,aes(x=POINT_X_lab, y=POINT_Y_lab,
                                  label = Systems), hjust = 0,vjust = -0.55, 
             family = "Helvetica",
             size = 7.3, fill = "#f7fcfd", label.size = 0, label.padding = unit(0, "lines"))+
  geom_point(data = data.fig4,
             aes(x=POINT_X, y=POINT_Y, fill =AC, size =storage_med_snow*100),
             shape = 21, color = "black", alpha = 1, stroke = 1.05)+
  scale_fill_stepsn("Exposure",
                    colors=(adaptivecapacity.pal),
                    n.breaks=10, limits = c(1,10), show.limits=T)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-.25,-.5,-.75,-1.0))+
  theme_map(base_size = 30)+
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

adaptivecapacity.map.fill = ggplot() +
  geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
               color = "grey85", size = 0.25, alpha = 0.55)+
  scale_fill_gradient2(low = "white", high = "grey20")+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig4,
             aes(x=POINT_X, y=POINT_Y, fill =AC),
             shape = 21, color = "black", alpha = 1, stroke = 1)+
  scale_fill_stepsn("Adaptation Index",
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
  theme(strip.text.x = element_text(size=0))+
  guides(fill = guide_colorsteps(frame.colour = "black", ticks=TRUE,
                                 ticks.colour = "black",
                                 ticks.linewidth = 1,
                                 title.position = "top"))

adaptivecapacity.legend = get_legend(adaptivecapacity.map.fill+
                                       theme(legend.position = "right",
                                                                     legend.title.align = 0.5))

adaptivecapacity.map.size = ggplot() +geom_polygon(data = wus, 
                                                   aes(long, lat, group = group), 
                                                   fill = "white",
                                                   color = "grey85", size = 0.25, 
                                                   alpha = 0.55)+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", 
          alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig4,
             aes(x=POINT_X, y=POINT_Y, size =storage_med_snow*100),
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_size_continuous(expression(paste("\U0394 in ", bar(S[Snow]), " [%]")),
                        range = c(15,5), breaks = c(-15,-30,-45,-60))+
  theme_map(base_size = 30)+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(size = 22, family  = "Helvetica"))+
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top",
                             label.position = "bottom")) +
  facet_wrap(~Period)+
  theme(strip.text.x = element_text(size=0))

adaptivecapacity.legend.size = get_legend(adaptivecapacity.map.size+
                                            theme(legend.position = "right",
                                                      legend.title.align = 0.5))
######################### ROSE BASED MAP OF ADAPTIVE CAPACITY ##########################
data.fig4.radar = data.fig4 %>%
  dplyr::select(System, Period, RI_delta, Demand_delta,
                storage_med,
                POINT_X, POINT_Y) %>%
  unique() %>% group_by(System, Period, POINT_X, POINT_Y, storage_med) %>%
  reshape2::melt(id = c("System", "Period", "storage_med", "POINT_X", "POINT_Y"))

data.fig4.radar = data.fig4.radar %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(# Northwest-blue
                                    "Little Wood" = "LW-ID         ",
                                    "Sun" = "         Su-MT",
                                    "Shoshone" = "         Sh-WY",
                                    "Wind" = "         Wi-WY",
                                    "Bridger" = "         Br-WY",
                                    "Price" = "         Pr-UT",
                                    "Paonia" = "         Pa-CO", # maybe tan 3
                                    "Costilla" = "Co-NM",
                                    # West-orangered
                                    "Kern" = "Ke-CA         ",
                                    "Kaweah" = "Ka-CA         ",
                                    "Walker" = "Wa-NV         ",
                                    "Umatilla" = "Um-OR         ",
                                    "Kittitas" = "Ki-WA         "
                                  )))

data.fig4.radar$Systems <- factor(data.fig4.radar$Systems,
                                  c(# Northwest-blue
                                    "LW-ID         ",
                                    "         Su-MT",
                                    "         Sh-WY",
                                    "         Wi-WY",
                                    "         Br-WY",
                                    "         Pr-UT",
                                    "         Pa-CO", # maybe tan 3
                                    "Co-NM",
                                    # West-orangered
                                    "Ke-CA         ",
                                    "Ka-CA         ",
                                    "Wa-NV         ",
                                    "Um-OR         ",
                                    "Ki-WA         "
                                  ))

df3 <- data.fig4.radar %>%
  dplyr::arrange((levels = Systems))

## Set up your radar plot for comparison
coord_radar <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}
df3 <- df3 %>%
  mutate(variable_1 = str_replace(variable, "RI_delta", "A"),
         variable_1 = str_replace(variable, "Demand_delta", "B"))

# data
## Now plot
windows()
ac.compass = ggplot(df3, aes(x=Systems, y=value*100, color=variable_1,group = variable_1),
                    alpha = 0.5) +
  geom_hline(aes(yintercept = 0), color = "red", alpha = 0.75)+
  geom_polygon(aes(color = variable_1, linetype = variable_1), fill = NA, size =1.35) +
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


######################### ROSE BASED MAP OF AC ##########################
legend_labels <- c(expression(frac("Supply Augmentation", S[Built])),
                   expression(frac("Water Conservation", S[Built])))

ac.compass.legend.plot =  ggplot(df3, aes(x=Systems, y=value, color=variable_1,
                                          group = variable_1),
                                 alpha = 0.5) +
  geom_line(aes(color = variable_1, linetype = variable_1), fill = NA,  size =1.5) +
  theme_minimal(base_size = 30)+
  coord_polar(start = 5.7, clip = "off")+
  scale_y_continuous(lim = c(-40,10))+
  scale_fill_manual("",
                    values = c("black","grey50" ))+
  scale_color_manual("Components of Adaptive Capacity Index",
                     values = c("black","grey50" ),
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

override.linetype <- c("dashed", "solid")

ac.compass.legend.plot <-ac.compass.legend.plot + guides(colour = guide_legend(override.aes = list(linetype = override.linetype),
                                                                               label.position = "bottom",
                                                                               title.position = "top", title.hjust = 0.5))

ac.compass.legend.plot <- ac.compass.legend.plot + scale_linetype(guide = FALSE)+
  theme(legend.title=element_text(size=24))

ac.compass.legend = get_legend(ac.compass.legend.plot+theme(legend.position = "right",
                                                            legend.title.align = 0.5))

## Now construct your plot
ac.plot = cowplot::ggdraw() +
  cowplot::draw_plot(adaptivecapacity.map,x=0-0.055,  y = 0.5-0.009,width = 1.1, height = 0.55) +
  cowplot::draw_plot(ac.compass,x=0,  y = 0.17-0.018,width = 0.95, height = 0.38) +
  cowplot::draw_plot(adaptivecapacity.legend,x=-0.28,  y = -0.3-0.028,width = 1.5, height = 0.48)+
  cowplot::draw_plot(adaptivecapacity.legend.size,x=0-0.59,  y = -0.29-0.028,width = 1.5, height = 0.48)+
  cowplot::draw_plot(ac.compass.legend,x=0.07,  y = -0.13-0.0441+0.009,width = 1.5, height = 0.48)+
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
                     family = "Helvetica", fontface = "bold", x= 0.7-0.02, y = 0.535+0.003)+
  cowplot::draw_text("\U2190 Lower", 
                     size = 22, angle = 0,
                     family = "Helvetica", fontface = "italic", x=0.360, y = 0.125)+
  cowplot::draw_text("Higher \U2192", 
                     size = 22, angle = 0,
                     family = "Helvetica", fontface = "italic", x= 0.58-0.0, y = 0.125)

ggsave(plot = ac.plot, "[LINK TO AC PLOT]", width=21, height=13, dpi=300)

###############################################################################
# ++++++++++++++++++ RISK AND THEORETICAL MINIMUM RISK ANALYSIS ++++++++++++++#
###############################################################################

## Select hazard and vulnerability from Figure 2-3
Risk.data.sub = data.fig3 %>% dplyr::select(System, Systems, Period, Haz, Vulernability,
                                            storage_med_snow, POINT_X, POINT_Y) %>%
  unique()

## And adaptive capcity from Figure 4
AC.data = data.fig4 %>% dplyr::select(Systems, System, Period, AC)

## Now join
Risk.data = left_join(Risk.data.sub, AC.data, by = c("Systems", "System", "Period"))

## Calculate Risk w/o adaptive capacity
Risk.data$Risk = (Risk.data$Haz+Risk.data$Vulernability)

## Calculate Risk w adaptive capacity
Risk.data$Risk_Min = (Risk.data$Haz+Risk.data$Vulernability)-
  Risk.data$AC

## Write your csv 
table2 = Risk.data %>% dplyr::select("Systems", "Period", "Haz", "Vulernability",
                                                  "AC", "Risk", "Risk_Min")

table2 = table2 %>%
  ungroup() %>%
  dplyr::mutate(A_risk = min(Risk_Min),
                B_risk = max(Risk))

## Rescale the conservation and storage indicators
table2 = table2 %>% 
  dplyr::mutate(Risk = 1+(10-1)*((Risk - A_risk)/(B_risk-A_risk)),
                Risk_Min =1+(10-1)*((Risk_Min - A_risk)/(B_risk-A_risk))) %>%
  dplyr::select("Systems", "Period", "Haz", "Vulernability", "AC", 
                "Risk", "Risk_Min")

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
# ++++++++++++++++++++++++++++++++ CODE FOR FIGURE 5  ++++++++++++++++++++++++#
###############################################################################

data.fig5 = Risk.data %>% ungroup()%>% group_by(System, Period) %>%
  dplyr::select(System, Period, POINT_X, POINT_Y, Risk_Min,
                Risk, storage_med_snow) %>%
  reshape2::melt(id = c("System", "Period", "POINT_X", "POINT_Y",
                        "storage_med_snow")) %>%
  ungroup() 

## Now melt for figure
data.fig5$variable_rename = ifelse(data.fig5$variable == c("Risk_Min"), "Risk Reduction",
                                   "Risk")

data.fig5 = transform(data.fig5,
                      Name=factor(variable_rename,levels=c("Risk", 
                                                           "Risk Reduction")))

## Get the minimum and maximum re-scaling values across all systems and periods
data.fig5 = data.fig5 %>%
  ungroup() %>%
  dplyr::mutate(A_risk = min(value),
                B_risk = max(value))

## Rescale 
data.fig5 = data.fig5 %>% 
  dplyr::mutate(Risk_I = 1+(10-1)*((value - A_risk)/(B_risk-A_risk)))

## Change the point labels
data.fig5$POINT_Y = ifelse(data.fig5$System == c("Shoshone"), 45,data.fig5$POINT_Y)

#Shift points
data.fig5$POINT_Y_lab = ifelse(data.fig5$System == c("Bridger"), 40.7,
                               ifelse(data.fig5$System == c("Costilla"), 34.5,
                                      ifelse(data.fig5$System == c("Kaweah"), 35.85,
                                             ifelse(data.fig5$System == c("Kern"), 33.0,
                                                    ifelse(data.fig5$System == c("Kittitas"), 47,
                                                           ifelse(data.fig5$System == c("Little Wood"), 41.3,
                                                                  ifelse(data.fig5$System == c("Paonia"),39.1,
                                                                         ifelse(data.fig5$System == c("Price"), 37.4,
                                                                                ifelse(data.fig5$System == c("Wind"),43.0,
                                                                                       ifelse(data.fig5$System == c("Sun"), 47.5,
                                                                                              ifelse(data.fig5$System == c("Umatilla"), 42.8,
                                                                                                     ifelse(data.fig5$System == c("Walker"), 38.5,
                                                                                                            45.1))))))))))))
data.fig5$POINT_X_lab = ifelse(data.fig5$System == c("Bridger"), -109.4,
                               ifelse(data.fig5$System == c("Costilla"), -113.5,
                                      ifelse(data.fig5$System == c("Kaweah"), -127.0,
                                             ifelse(data.fig5$System == c("Kern"), -126.1,
                                                    ifelse(data.fig5$System == c("Kittitas"), -127.5,
                                                           ifelse(data.fig5$System == c("Little Wood"), -121.3,
                                                                  ifelse(data.fig5$System == c("Paonia"),-107.0,
                                                                         ifelse(data.fig5$System == c("Price"), -117.7,
                                                                                ifelse(data.fig5$System == c("Wind"),-108.3,
                                                                                       ifelse(data.fig5$System == c("Sun"), -112.2486,
                                                                                              ifelse(data.fig5$System == c("Umatilla"), -127.5,
                                                                                                     ifelse(data.fig5$System == c("Walker"), -127.5,
                                                                                                            -109))))))))))))


############################# Figure 5 a-f #####################################
Risk.map = ggplot() +
  geom_polygon(data = wus, aes(long, lat, group = group),
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
             aes(x=POINT_X, y=POINT_Y, fill =Risk_I, size =storage_med_snow*100), 
             shape = 21, color = "black", alpha = 1, stroke = 1.05)+
  scale_fill_stepsn("Relative Resilience or Risk",
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

## Construct your own risk map fill
Risk.map.fill = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                                color = "grey85", size = 0.25, alpha = 0.55)+
  scale_fill_gradient2(low = "white", high = "grey20")+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig5,
             aes(x=POINT_X, y=POINT_Y, fill =Risk_I),
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Relative Risk\nor Resilience\n",
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

Risk.legend = get_legend(Risk.map.fill)

## Custom map size scale
Risk.map.size = ggplot() +geom_polygon(data = wus, aes(long, lat, group = group), fill = "white",
                                                color = "grey85", size = 0.25, alpha = 0.55)+
  geom_sf(data = regions, fill = "grey77", color = "grey55", alpha = 0.2) +
  geom_sf(data = rivers, fill = "cornflowerblue", color = "cornflowerblue", alpha = 1)+
  new_scale_fill() +
  geom_point(data = data.fig5,
             aes(x=POINT_X, y=POINT_Y, size = storage_med_snow*100),
             shape = 21, color = "black", alpha = 0.75, stroke = 1)+
  scale_fill_stepsn("Risk",
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

Risk.legend.size = get_legend(Risk.map.size+theme(legend.position = "right", 
                                                                    legend.title.align = 0.5))

############################# Figure 5 g-i #####################################
data.fig5.radar = data.fig5 %>% 
  dplyr::select(System, Period,  Risk_I, variable_rename) %>%
  unique() %>% 
  tidyr::spread(key = variable_rename, value = Risk_I)

names(data.fig5.radar) = c("System", "Period", "Risk for Climate Change", "Reduction in Risk from Adaptive Capacity")

data.fig5.radar$Difference = abs(data.fig5.radar$`Risk for Climate Change`-data.fig5.radar$`Reduction in Risk from Adaptive Capacity`)
data.fig5.radar$Difference_per = (abs(data.fig5.radar$`Risk for Climate Change`-data.fig5.radar$`Reduction in Risk from Adaptive Capacity`)/data.fig5.radar$`Risk for Climate Change`)*100

# Construct your radar plot
data.fig5.radar = data.fig5.radar %>%
  mutate(Systems =  plyr::revalue(System,
                                  c(# Northwest-blue
                                    "Little Wood" = "LW-ID         ",
                                    "Sun" = "         Su-MT",
                                    "Shoshone" = "         Sh-WY",
                                    "Wind" = "         Wi-WY",
                                    "Bridger" = "         Br-WY",
                                    "Price" = "         Pr-UT",
                                    "Paonia" = "         Pa-CO", # maybe tan 3
                                    "Costilla" = "Co-NM",
                                    # West-orangered
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

data.fig5.radar$Systems <- factor(data.fig5.radar$Systems,
                                  c(# Northwest-blue
                                    "LW-ID         ",
                                    "         Su-MT",
                                    "         Sh-WY",
                                    "         Wi-WY",
                                    "         Br-WY",
                                    "         Pr-UT",
                                    "         Pa-CO", # maybe tan 3
                                    "Co-NM",
                                    # West-orangered
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

df5 <- data.fig5.radar %>%
  dplyr::arrange((levels = Systems))


Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")

## Set up your radar plot for comparison
coord_radar <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

## Set up your data and factor as a rename
df5 <- df5 %>% dplyr::select(Systems, Period, 'Risk for Climate Change', Difference)
df5 <- data.fig5.radar %>%
  dplyr::arrange((levels = Systems))

## Now plot
vulnerability.rose =ggplot(df5) +
  geom_col(aes(x=Systems,
               y=Difference, fill = `Risk for Climate Change`),
           position = "dodge2", show.legend = TRUE, alpha = 0.9, color = "black")+
  coord_polar(start = 5.7, clip = "off")+
  facet_wrap(Period ~ ., labeller = label_wrap_gen(width=14), nrow =1) +
  ylab("")+xlab("")+  
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c("bottom"))+theme(legend.direction="horizontal")+
  scale_fill_stepsn("Risk Index",
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

## Construct your plot
vulnerability = cowplot::ggdraw() +
  cowplot::draw_plot(vulnerability.map,x=0-0.15,  y = 0.275-0.005,width = 1.1, height = 0.76) +
  cowplot::draw_plot(vulnerability.rose,x=0-0.115,  y = 0.0-0.005,width = 1.0, height = 0.28)  +
  cowplot::draw_plot(vulnerability.legend.size,x=0.15,  y = -0.25,width = 1.5, height = 0.48)+
  cowplot::draw_plot(vulnerability.legend,x=0.15,  y = 0.5,width = 1.5, height = 0.48)+
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
                     font = "Helvetica", fontface = "bold", x= 0.55-0.01, y = 0.30)+
  cowplot::draw_label(expression(paste("Risk for")),
                      size = 30, angle = -90,
                      fontfamily = "Helvetica", x= 0.795, y = 0.8)  +
  cowplot::draw_label(expression(paste("\U0394 in ", bar(S[Snow]))),
                      size = 30, angle = -90,
                      fontfamily = "Helvetica", x= 0.77, y = 0.8)+
  cowplot::draw_label(expression(paste("Potential Minimum Risk for")),
                      size = 30, angle = -90,
                      fontfamily = "Helvetica", x= 0.795, y = 0.47)  +
  cowplot::draw_label(expression(paste("\U0394 in ", bar(S[Snow]))),
                      size = 30, angle = -90,
                      fontfamily = "Helvetica", x= 0.77, y = 0.47)  +
  cowplot::draw_text("Difference",
                     size = 30, angle = -90,
                     family = "Helvetica", x= 0.795, y = 0.15)+
  cowplot::draw_text("\U2190 Higher", 
                     size = 22, angle = -90,
                     family = "Helvetica", fontface = "italic", x=0.925, y = 0.875)+
  cowplot::draw_text("Lower \U2192", 
                     size = 22, angle = -90,
                     family = "Helvetica", fontface = "italic", x= 0.925, y = 0.32)

ggsave(plot = Risk, "[LINK TO FIGURE 5]", width=21, height=13, dpi=300)

