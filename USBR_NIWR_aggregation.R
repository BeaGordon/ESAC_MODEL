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
options(scipen=999)
###############################################################################
# ++++++++++++++++++++ READ in REQ FUNCTIONS++++++++++++++++++++++++++++++ #
###############################################################################
hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}

###############################################################################
# ++++++++++++++++++++ READ in the ET files you need  +++++++++++++++++++++++++ #
###############################################################################
dir <- "[LINK TO DIRECTORY WHERE USBR ET DATA ARE STORED FROM USBR_EXTRACTION.R]"
setwd(dir)

list = list.files(path = dir, 
                  pattern=glob2rx("*ET*"), # Input GCM name here
                  recursive = TRUE)
lapply.dat = lapply(list, read.csv)
data.name = sub(".*/", "", list)
data.name = sub(".csv", "", data.name)
data.name = sub("ET_", "", data.name)
site.name = sub("/.*", "", list)

lapply.dat = Map(cbind, lapply.dat, Site.ID = site.name, Data = data.name)
data= list.rbind(lapply.dat) # create a dataframe of all lists
names(data) = c("X", "twoweeks", "ET_Var", "BoR.crop", "med.ETact", "ub.ETact", "lb.ETact",
                "Site.ID", "Period")

data$Period = ifelse(data$Period == c("2020"), "2020-2050", 
                     ifelse(data$Period == c("2050"), "2050-2080", 
                            ifelse(data$Period == c("2080"), "2080-2100",
                                   "baseline")))

Crop_LU_1 = read.csv("[READ IN LOOKUP TABLE CROP_LU_V2.CSV FROM LOOKUP TABLE FOLDER")[,-1]
data = left_join(data, Crop_LU_1, by = c("BoR.crop")) %>%
  dplyr::select(Site.ID, Period, Crop, twoweeks, ET_Var, med.ETact, ub.ETact, lb.ETact) 

###############################################################################
# ++++++++++++++++++++ READ in LU data from cropscape +++++++++++++++++++++++ #
###############################################################################
landuse.uncertainty = read.csv("[LINK TO PROCESSED CROPSCAPE DATA FROM THE CROPSCAPE FOLDER]",
                               stringsAsFactors = FALSE)[,-1]
landuse.uncertainty = landuse.uncertainty %>% 
  dplyr::select("Site.ID", "System", "Year", "BoR.crop", "Area.adj")
landuse.uncertainty = landuse.uncertainty %>%
  mutate(System =  plyr::revalue(System,
                                 c(# Northwest-blue
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
                                   "Yakima" = "Ki-WA"
                                 )))

names(landuse.uncertainty) = c("Site.ID", "System", "Year", "Crop", "Area.adj")
landuse.uncertainty.agregate = landuse.uncertainty %>%
  group_by(Year, System, Crop) %>%
  dplyr::summarize(Acreage= sum(Area.adj)) %>% na.omit()
large = c("Ka-CA" , "Ke-CA", "Ki-WA")
unique(landuse.uncertainty.agregate$System)
landuse.uncertainty.agregate.small = landuse.uncertainty.agregate %>% 
  dplyr::filter(System != "Ka-CA") %>%
  dplyr::filter(System != "Ke-CA")%>%
  dplyr::filter(System != "Ki-WA") %>% na.omit()

landuse.uncertainty.agregate.large = landuse.uncertainty.agregate %>% 
  dplyr::filter(System == "Ka-CA" |System == "Ke-CA" | System == "Ki-WA") %>% na.omit()


###############################################################################
# +++++++++++++++++++ GENERATE PLOTS OF LU FOR SI FIGURES +++++++++++++++++++ #
###############################################################################
library(randomcoloR)
n <- 33
palette <- distinctColorPalette(n)
unique(landuse.uncertainty.agregate$Crop)

## Constuct plot for SI
plot.si.1 = ggplot(data =landuse.uncertainty.agregate.small) +
  geom_col(aes(x= Year, y = Acreage,
               fill = Crop,
               group =Crop), size = 1)+
  theme_minimal(base_size = 24)+xlab("Year") +
  scale_fill_manual(values = palette)+
  ylab("Bias Corrected Acreage by Crop (Acres)")+
  theme(plot.title=element_text(size=30, face="bold"),
        axis.text.x=element_text(size=15, angle=90),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=25),
        axis.title.y=element_text(size=25))+theme(plot.title = element_text(hjust = 0.5))+
  scale_color_brewer("Crop", palette = "Spectral")+
  theme(legend.direction = "horizontal", legend.position = "bottom")+
  facet_wrap(~System, labeller = label_wrap_gen(width=21))

plot.1 = cowplot::ggdraw() +
  cowplot::draw_plot(plot.si.1,x=0,  y = 0,width = 1, height = 1) +
  cowplot::draw_text("A", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.1, y = 0.96)+
  cowplot::draw_text("B", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.38-0.035, y = 0.96)+
  cowplot::draw_text("C", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.59-0.02, y = 0.96)+
  cowplot::draw_text("D", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.85-0.05, y = 0.96)+
  cowplot::draw_text("E", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.1, y = 0.7)+
  cowplot::draw_text("F", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.38-0.035, y = 0.7)+
  cowplot::draw_text("G", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.59-0.02, y = 0.7)+
  cowplot::draw_text("H", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.85-0.05, y = 0.7)+
  cowplot::draw_text("I", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.1, y = 0.46)+
  cowplot::draw_text("J", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.38-0.035, y = 0.46)

ggsave(plot = plot.1, "[LINK TO FOLDER FOR LAND USE SUPPLEMENTAL FIGURE #1]", 
       width=21, height=13, dpi=300)

plot.si.2 = ggplot(data =landuse.uncertainty.agregate.large) +
  geom_col(aes(x= Year, y = Acreage,
               fill = Crop,
               group =Crop), size = 1)+
  theme_minimal(base_size = 24)+xlab("Year") +
  scale_fill_manual(values = palette)+
  ylab("Bias Corrected Acreage by Crop (Acres)")+
  theme(plot.title=element_text(size=30, face="bold"),
        axis.text.x=element_text(size=15, angle=90),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=25),
        axis.title.y=element_text(size=25))+theme(plot.title = element_text(hjust = 0.5))+
  scale_color_brewer("Crop", palette = "Spectral")+
  theme(legend.direction = "horizontal", legend.position = "bottom")+
  facet_wrap(~System, labeller = label_wrap_gen(width=21))

plot.2 = cowplot::ggdraw() +
  cowplot::draw_plot(plot.si.2,x=0,  y = 0,width = 1, height = 1) +
  cowplot::draw_text("A", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.1, y = 0.96)+
  cowplot::draw_text("B", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x=0.45-0.035, y = 0.96)+
  cowplot::draw_text("C", 
                     size = 26, angle = 0,
                     family = "Helvetica", fontface = "bold", x= 0.74-0.02, y = 0.96)
ggsave(plot = plot.2, "[LINK TO FOLDER FOR LAND USE SUPPLEMENTAL FIGURE #2]", width=21, height=13, dpi=300)

###############################################################################
# ++++++++++++++++++++ Analyze the et demand, combine +++++++++++++++++++++++ #
###############################################################################
## A. First filter out the odd Bridger years
data.bor.cs$Screen = ifelse(data.bor.cs$System == c("Bridger")&
                              data.bor.cs$Year %in% c("2011", "2012"), "Flag", "No Flag")
data.bor.cs$Screen = ifelse(data.bor.cs$System == c("Paonia")&
                              data.bor.cs$Year %in% c("2012"), "Flag", data.bor.cs$Screen)
data.bor.cs$Screen = ifelse(data.bor.cs$System == c("Costilla")&
                              data.bor.cs$Year %in% c("2017", "2018"), "Flag", data.bor.cs$Screen)

## 1. Add your land use estimates + irrigation demand
mm_2_m = 1000
data.bor.cs.ag = data.bor.cs %>% na.omit %>% dplyr::filter(Period == c("baseline"))%>%
  dplyr::filter(Screen != c("Flag")) %>%
  group_by(Site.ID, Year, Crop, twoweeks, System) %>%
  dplyr::summarize(Use = Area.adj*(med.ETact/mm_2_m))

## 3. Now obtain the cumulative sum
data.bor.cs.ag.sum = data.bor.cs.ag %>% ungroup() %>%
  group_by(Site.ID,Year, Crop, System) %>%
  dplyr::summarize(Volume = sum(Use))

# ## QAQC on crops
unique(data.bor.cs.ag.sum$System)
test = data.bor.cs.ag.sum %>% dplyr::filter(System == c("Price"))

ggplot(data = test) + geom_col(aes(x = Crop, y = Volume, fill = as.factor(Year)),
                               position = "dodge2") +facet_wrap(~Site.ID)

## 4. Use the annual volume to choose the crop mix
data.bor.cs.ag.sum = data.bor.cs.ag.sum %>% ungroup() %>%
  group_by(Site.ID,Year, System) %>%
  dplyr::mutate(Total.Volume = sum(Volume))%>% ungroup() %>%
  group_by(Site.ID, System) %>%
  dplyr::mutate(Volume.75 = quantile(Total.Volume, .75),
                Volume.50 = quantile(Total.Volume, .50),
                Volume.25 = quantile(Total.Volume,0.25)) 

data.bor.cs.ag.sum$difference.75 = abs(data.bor.cs.ag.sum$Total.Volume-data.bor.cs.ag.sum$Volume.75)
data.bor.cs.ag.sum$difference.25 = abs(data.bor.cs.ag.sum$Total.Volume-data.bor.cs.ag.sum$Volume.25)
data.bor.cs.ag.sum$difference.50 = abs(data.bor.cs.ag.sum$Total.Volume-data.bor.cs.ag.sum$Volume.50)

## 5. Scenario cropping selection
data.bor.scen.75 = data.bor.cs.ag.sum %>% 
  dplyr::mutate(min = abs(min(difference.75))) %>% 
  dplyr::filter(difference.75 == min) %>%
  dplyr::select("Site.ID", "Year", "Crop", "System") %>%
  dplyr::mutate(Scenario.High = "High.Land.Use")
check = data.bor.scen.75 %>% group_by(Site.ID,Scenario.High, Crop) %>% tally()
check = data.bor.scen.75 %>% group_by(Scenario.High, Site.ID) %>% tally()
data.bor.scen.50 = data.bor.cs.ag.sum %>% 
  dplyr::mutate(min = abs(min(difference.50))) %>% 
  dplyr::filter(difference.50 == min) %>%
  dplyr::select("Site.ID", "Year", "Crop", "System")%>%
  dplyr::mutate(Scenario.Med = "Median.Land.Use")
check = data.bor.scen.50 %>% group_by(Site.ID,Scenario.Med, Crop) %>% tally() 
check = data.bor.scen.50 %>% group_by(Scenario.Med, Site.ID) %>% tally()

data.bor.scen.25 = data.bor.cs.ag.sum %>% 
  dplyr::mutate(min = abs(min(difference.25))) %>% 
  dplyr::filter(difference.25 == min) %>%
  dplyr::select("Site.ID", "Year", "Crop", "System")%>%
  dplyr::mutate(Scenario.Low = "Low.Land.Use")
check = data.bor.scen.25 %>% group_by(Site.ID,Scenario.Low, Crop) %>% tally() 
check = data.bor.scen.25 %>% group_by(Scenario.Low, Site.ID) %>% tally()

data.bor.scen.25$Flag = ifelse(data.bor.scen.25$System == c("Price") & 
                                 data.bor.scen.25$Year == c("2011"),
                               "Flag", "No Flag")

data.bor.scen.25 = data.bor.scen.25 %>% 
  dplyr::filter(Flag != c("Flag"))

###############################################################################
# ++++++++++++++++++++ Set up and graph your scenarios +++++++++++++++++++++++ #
###############################################################################
## 1. Join your DF with your scenarios
data.scenario.bor = full_join(data.bor.cs, data.bor.scen.25, 
                              by = c("Site.ID", "Year", "Crop", "System"))
data.scenario.bor = full_join(data.scenario.bor, data.bor.scen.50, 
                              by = c("Site.ID", "Year", "Crop", "System"))
data.scenario.bor = full_join(data.scenario.bor, data.bor.scen.75, 
                              by = c("Site.ID", "Year", "Crop", "System"))

data.scenario.bor = data.scenario.bor %>% dplyr::select(-c("Flag")) %>%
  dplyr::mutate(Scenario = coalesce(Scenario.Low,Scenario.Med,Scenario.High)) %>%
  dplyr::select(-c("Scenario.Low", "Scenario.Med", "Scenario.High", "Year")) %>%
  na.omit()

check = data.scenario.bor %>% group_by(Site.ID,Scenario, Period) %>% tally() 

## 3. Now convert to m, obtain the cumulative sum
data.scenario = data.scenario.bor %>% 
  dplyr::select("System", "Site.ID", "Period", "twoweeks", "Crop", "Scenario", "Area.adj",
                "med.ETact", "ub.ETact", "lb.ETact")
data.scenario$med.ETact = data.scenario$med.ETact/mm_2_m
data.scenario$ub.ETact = data.scenario$ub.ETact/mm_2_m
data.scenario$lb.ETact = data.scenario$lb.ETact/mm_2_m

## Format as a pivot wider
data.scenario = data.scenario %>% 
  tidyr::pivot_wider(names_from = Scenario, values_from = Area.adj)

## Convert to NAs and replace
data.scenario$Median.Land.Use = as.numeric(as.character(data.scenario$Median.Land.Use))
data.scenario$Low.Land.Use = as.numeric(as.character(data.scenario$Low.Land.Use))
data.scenario$High.Land.Use = as.numeric(as.character(data.scenario$High.Land.Use))
data.scenario[is.na(data.scenario)] = 0

## 4. Now convert to m, obtain the cumulative sum
data.scenario.csum = data.scenario %>% 
  group_by(System, Site.ID, Period,Crop) %>%
  dplyr::mutate(HighLU.HighWU = cumsum(ub.ETact*High.Land.Use),
                MedLU.MedWU = cumsum(med.ETact*Median.Land.Use),
                LowLU.LowWU = cumsum(med.ETact*Low.Land.Use))

# el.lu = read.csv("D:/Projects/Project_Reservoir/Data/Spreadsheets/LU/Elastic_Inelastic.csv")
# 
data.scenario.csum.m = data.scenario.csum %>% dplyr::select(-c("med.ETact",
                                                               "ub.ETact",
                                                               "lb.ETact",
                                                               "High.Land.Use",
                                                               "Median.Land.Use",
                                                               "Low.Land.Use")) %>%
  reshape2::melt(id = c("System", "Site.ID", "Period", "twoweeks", "Crop"))

###############################################################################
# ++++++++++++++++++ Join Scenarios with ETO and plot +++++++++++++++++++++++ #
###############################################################################
## 1. Curate your BoR dataset
data.scenario.bor.j = data.scenario.csum.m %>% ungroup() %>%
  group_by(System, Period, twoweeks, Crop, variable) %>%
  dplyr::summarize(value = sum(value))

names(data.scenario.bor.j) = c("System","Period", "twoweeks", "Crop",
                               "Scenario", "ETa_m3")

write.csv(unique(data.scenario.bor.j$Crop), "[LINK TO CROP NAMES FILE]")

##2. Clean up the names 
data.scenario.bor.j$Scenario = str_replace(data.scenario.bor.j$Scenario, "HighLU.HighWU", "High Water Use")
data.scenario.bor.j$Scenario = str_replace(data.scenario.bor.j$Scenario, "MedLU.MedWU", "Business-as-Usual")
data.scenario.bor.j$Scenario = str_replace(data.scenario.bor.j$Scenario, "LowLU.LowWU", "Low Water Use")


data.scenario.bor.j$screen = ifelse(data.scenario.bor.j$Scenario == c("High Water Use") & 
                                      data.scenario.bor.j$Period == c("baseline"), 0,
                                    ifelse(data.scenario.bor.j$Scenario == c("Low Water Use") & 
                                             data.scenario.bor.j$Period == c("baseline"), 0, 1))

data.scenario.bor.j$Scenario = ifelse(data.scenario.bor.j$Scenario == c("Business-as-Usual") & 
                                        data.scenario.bor.j$Period == c("baseline"), "Historical",
                                      data.scenario.bor.j$Scenario)

data.scenario.bor.j$Period = ifelse(data.scenario.bor.j$Scenario == c("Historical") & 
                                      data.scenario.bor.j$Period == c("baseline"), "1950-1999",
                                    data.scenario.bor.j$Period)

## 3. Sub out for matching
data.scenario.bor.j = data.scenario.bor.j %>% dplyr::filter(screen !=0)

data.scenario.bor.j = data.scenario.bor.j %>% 
  dplyr::select(System, Period, twoweeks, Scenario, ETa_m3, Crop)

names(data.scenario.bor.j) = c("System", "Period", "twoweeks", "Scenario", "ET_m3", 
                               "Category")

## 4. Read in your ETo dataset
data.ETo.sys = read.csv("[READ IN YOUR ETO DATASET OBTAINED FROM THE RE_PROCESSING CODE]")[,-1]
data.ETo.sys$Category = "Reservoir Evaporation"
names(data.ETo.sys) = c("System", "Period", "twoweeks", "Scenario", "ET_m3", 
                        "Category")

## 5. Combine 
data.scenario.full = rbind(data.scenario.bor.j, data.ETo.sys)

write.csv(data.scenario.full, 
          "[LINK TO YOUR FULL DEMAND DATA]")

data.hist = data.scenario.full %>% dplyr::filter(Period == c("1950-1999"))

write.csv(data.hist, 
          "[LINK TO YOUR HISTORICAL ONLY DEMAND DATA]")


write.csv(data.bor.cs.ag, "[LINK TO YOUR FULL SYSTEM RESULTS DEMAND DATA]")

