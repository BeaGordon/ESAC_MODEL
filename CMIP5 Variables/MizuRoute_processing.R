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
library(exactextractr)
library(tidync)

###############################################################################
# ++++++++++++++++++++ READ in STREAMFLOW DATA ++++++++++++++++++++++++++++++ #
###############################################################################

dir = "[LINK TO TAR OR ZIP FILE OF ROUTED STREAMFLOW" # from the .tar
setwd(dir)

## Read in the time series information
files = list.files(path=setwd(dir),
                   pattern=glob2rx("*.nc"), recursive = TRUE)

## Setup the files.list function
files.list = function(files) {
  
  files = files # assigns files to files
  file_nonc = str_remove(files, ".nc") # removes .nc extension to use nc_open
  fname = nc_open(files)
  datetime = nc.get.time.series(fname) # stores the datetime according to the .nc
  sSeg = as.data.frame(fname$dim$sSeg$vals) # extract the sSeg variable names
  datetime = toString(datetime) # datetime to string
  datetime.spl <- unlist(strsplit(datetime, ",")) # formats the datetime string
  datetime.spl = as.data.frame(datetime.spl) # datetime string to df
  names(datetime.spl) = "datetime" # column name
  datetime.spl$datetime = as.Date(datetime.spl$datetime, format = c("%Y-%m-%d")) 
  ``
  ## Read in the values
  src <- tidync(files) # connect to data source and obtain metadata
  data = src %>% tidync::hyper_array() # select output as array
  df = data$KWTroutedRunoff # subset data based on desired routing option
  df = as.data.frame(df) # make df
  df.t = as.data.frame(t(df)) # transpose df
  names = sSeg # assign ordered numeric value based on sSeg variable
  names = paste0("Site_", names$`fname$dim$sSeg$vals`)# paste Site for clarity
  colnames(df.t) = names # assign correct column names
  
  ## Add in the time
  df.t$datetime = datetime.spl$datetime # replace df.t datetime with formatted raw datetime
  
  ## Melt to get into long format
  df.melt = df.t %>% reshape2::melt(id = c("datetime")) # melt based on df.t
  names(df.melt) = c("dates", "site", "streamflow") # assign the correct names
  
  write.csv(df.melt, file=paste("[DIRECTORY FOR PROCESSED STREAMFLOW", file_nonc, ".csv", sep="")) # write file
}
lapply(files, files.list)
