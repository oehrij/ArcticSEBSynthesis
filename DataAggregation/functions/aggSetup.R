#########################################################
#####
##### project:       ArcticSEBSynthesis
##### title:         aggSetup
##### description:   setup script for aggregation of data
##### publication:   Oehri, J. et al. (2022). Vegetation Type is an Important Predictor of the Arctic Summer Land Surface Energy Budget. Nature Communications.
##### author:        Jacqueline Oehri
##### date:          28.09.2022  
##### comments:      Run this script before aggregating data
#########################################################

##########################################################
## clean space
rm(list=ls(all=TRUE))

##########################################################
## libraries
library(pascal)         #library(devtools); install_github("pascal-niklaus/pascal/pascal")
library(pgeo)           #library(devtools); install_github("pascal-niklaus/pascal/pascal")
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(magick)
library(readxl)
library(raster)
library(ncdf4)
library(lattice)
library(chron)
library(lutz)           #tz_lookup_coords(lat, lon, method = "accurate", warn = TRUE)
library(doParallel)
library(zoo)            #for the rollmean function

##########################################################
## functions on which the below depend...
### colors
col   <- brewer.pal(7,"Accent")
crf   <- colorRampPalette(col)
qccol <- c("darkblue","darkorange")
crf1  <- colorRampPalette(qccol)

### make one dataframe
load_data <- function(path, pattern) { 
  files <- dir(path, pattern = pattern, full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

## add legend outside plot margins
add_legend <- function(...) {
  opar <- par(fig=c(0.02, 0.97, 0.02, 0.97), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

################################################################################
## Preparing the source data to be processed: timestamp variables: (including "TIMESTAMP")
## example: add timestamps FLUXNET format YYYYMMDDHHMM
# d2$YEAR      <-  as.numeric(gsub("^([0-9]{4})([0-9]{8})","\\1",d2[,grep("TIMESTAMP",colnames(d2))[1]]))
# d2$MONTH     <-  as.numeric(gsub("^([0-9]{4})([0-9]{2})([0-9]{6})$","\\2",d2[,grep("TIMESTAMP",colnames(d2))[1]]))
# d2$DAY       <-  as.numeric(gsub("^([0-9]{6})([0-9]{2})([0-9]{4})$","\\2",d2[,grep("TIMESTAMP",colnames(d2))[1]]))
# d2$HOUR      <-  as.numeric(gsub("^([0-9]{8})([0-9]{2})([0-9]{2})$","\\2",d2[,grep("TIMESTAMP",colnames(d2))[1]]))
# d2$HRMINUTE  <- d2$HOUR+(as.numeric(gsub("^([0-9]{10})([0-9]{2})$","\\2",d2[,grep("TIMESTAMP",colnames(d2))[1]]))/60)
# d2$DOY       <-  as.numeric(format(strptime(as.character( d2[,grep("TIMESTAMP",colnames(d2))[1]] ) ,"%Y%m%d%H%M"),'%j'))
# d2$DOYdec    <-  as.numeric(format(strptime(as.character( d2[,grep("TIMESTAMP",colnames(d2))[1]] ) ,"%Y%m%d%H%M"),'%j')) + (d2$HOUR)/24
#
################################################################################
## end
################################################################################
