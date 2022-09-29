#########################################################
#####
##### project:       ArcticSEBSynthesis
##### title:         aggFLUXNET_example
##### description:   Example of aggregation of FLUXNET data on surface energy budget (SEB) components from (half-)hourly to daily and monthly resolution and store them 
##### publication:   Oehri, J. et al. (2022). Vegetation Type is an Important Predictor of the Arctic Summer Land Surface Energy Budget. Nature Communications.
##### author:        Jacqueline Oehri
##### date:          28.09.2022  
##### comments:      last run: JO, 14.01.2021, files: FLUXNET_SITES_60N.txt, variable_codes_FULLSET_20200504.csv
##########################################################

##########################################################
## load setup and other functions
## directory where aggregation functions are stored
fundir <- "C:/Users/jacqu/P1_Project/functions"
source(sprintf("%s/aggSetup.R",fundir))
source(sprintf("%s/aggFunctions.R",fundir))

##########################################################
## set directories (these directories may need to be adjusted)
## directory where data is stored
datdir    <- "C:/Users/jacqu/P1_Project/data/FLUXNET"
## directory with large space for saving aggregated data
dirlarge  <- "D:/Project_B1_Appendix/SEB_data_integration/FLUXNET"
## directory for saving logfiles
for( di in c("logs_aggav")){
  if( ! dir.exists(sprintf("%s/%s",dirlarge,di)) ){ 
    dir.create(sprintf("%s/%s",dirlarge,di))
  }}

################################################################################
## DATA processing start                                                      
################################################################################
## read in data FLUXNET: preparation
## currently FLUXNET directory contains all station data >60°N latitude
## https://fluxnet.fluxdata.org/sites/site-list-and-pages/?view=table

## Site information
fs <- read.delim(sprintf("%s/FLUXNET_SITES_60N.txt",datdir))
fs$LOCATION_LAT<- as.numeric(gsub("[,]",".",fs$LOCATION_LAT))
fs$LOCATION_LONG  <- as.numeric(gsub("[,]",".",fs$LOCATION_LONG))
print(fs)

## variable information
fv <- read.csv(sprintf("%s/variable_codes_FULLSET_20200504.csv",datdir))
print(fv)

## explicit variables of interest:
vars <- c("TIMESTAMP","TIMESTAMP_START","TIMESTAMP_END",
          "TA_F_MDS","TA_F_MDS_QC",
          "SW_IN_POT", "SW_DIF","SW_DIF_QC",
          "SW_IN_F_MDS","SW_IN_F_MDS_QC","LW_IN_F_MDS","LW_IN_F_MDS_QC",
          "SW_OUT","SW_OUT_QC","LW_OUT","LW_OUT_QC",
          "NETRAD", "NETRAD_QC",
          "H_F_MDS","H_F_MDS_QC",
          "LE_F_MDS","LE_F_MDS_QC",
          "G_F_MDS","G_F_MDS_QC",
          "VPD_F_MDS","VPD_F_MDS_QC",
          "PA","P","WS","WD","RH",
          "USTAR","USTAR_QC",
          "TS_F_MDS_#","TS_F_MDS_#_QC","SWC_F_MDS_#","SWC_F_MDS_#_QC")

print(fv[fv$Variable %in% vars,c("Variable","Description")])
## nr variables
length(vars)

## variables of interest made for "grep" selection
varsg <- c("^TIMESTAMP",
          "TA_F_MDS",
          "SW_IN_POT", "SW_DIF", 
          "SW_IN_F_MDS","SW_IN","SW_OUT",
          "LW_IN_F_MDS","LW_IN","LW_OUT",
          "NETRAD",
          "H_F_MDS","LE_F_MDS","G_F_MDS",
          "^H$","^LE$","^G$",
          "VPD_F_MDS",
          "^PA$","^PA_F$","^P$","^P_F$","^WS$","^WS_F$","WD","RH",
          "TS_F_MDS_[0-9]","SWC_F_MDS_[0-9]")

################################################################################
## read in data FLUXNET:

## Sites
sites <- fs$SITE_ID

## Vars of interest
vs   <- varsg 

## Times of interest
# for the moment use only HH or HR (see notes below)

## list directories
## because files are zipped this does not work
dirs <- list.files(datdir,pattern="[.]zip")

## Metadata (e.g. author, variable height and instrumentation!!) to all sites in folder labelled "FLX_AA"
mdir  <- dirs[grep("FLX_AA",dirs)]
fdirs <- dirs[-which(dirs==mdir)]

simple.heading(sprintf("nr of fluxnet sites = %d",length(fdirs)))


##########################################
## go through different sites (in different fdirs)

for (d in fdirs) {
  ## set up correct names
  siteid   <-  gsub("^FLX[_](.+)[_]FLUXNET2015_FULLSET[_](.+)[_].+[.]zip$","\\1",d)
  yearsite <-  gsub("^FLX[_](.+)[_]FLUXNET2015_FULLSET[_](.+)[_].+[.]zip$","\\2",d) 
  simple.heading(sprintf("file %d site ID: %s, years %s",which(fdirs==d),siteid, yearsite ))
  coords <- fs[grep(siteid,fs$SITE_ID),c("LOCATION_LAT","LOCATION_LONG")]
  
  ## Timezone lookup
  #tz_lookup_coords(lat, lon, method = "accurate", warn = TRUE)
  TZ <- tz_lookup_coords(coords[1,1], coords[1,2], method = "accurate", warn = TRUE)
  
  files <- as.character(unzip(sprintf("%s/%s",datdir,d), list = TRUE)$Name)
  fs1    <- files[grep("FULLSET",files)]

  ## work with the hourly data!
  t <- "HH"
  if(length(grep(t,fs1))==0) {
  t<-"HR"
  }
  
  ## read in the data
  simple.heading(sprintf("temporal resolution = %s",t))
  simple.heading(print(fs1[grep(t,fs1)]))
  ## takes kind of long: think of using other function to read in data?  
  data <- read.csv(unz(sprintf("%s/%s",datdir,d),fs1[grep(t,fs1)]), header = TRUE, sep = ",") 
  ## which variables cannot be found in the dataset
  simple.heading("vars of interest NOT found in data")
  print(vars[which(!vars %in% colnames(data))])
  ## select all variables of interest
  idx <- c()
  for(vv in vs) {
    idx <- c(idx,grep(vv, colnames(data)))
  }
  idx <- unique(idx)
  simple.heading("vars of interest SELECTED from data")
  print(colnames(data)[idx])
  # some ERA downscaled measurements are still in there
  
  ## extract only variables of interest (& set naval NA)
  d1 <- data[,idx]
  d2 <- d1
  d2[d2 == -9999] <- NA
   
  ## add timestamps
  # FLUXNET format YYYYMMDDHHMM
  d2$YEAR     <- as.numeric(gsub("^([0-9]{4})([0-9]{8})","\\1",d2[,grep("TIMESTAMP",colnames(d2))[1]]))
  d2$MONTH    <- as.numeric(gsub("^([0-9]{4})([0-9]{2})([0-9]{6})$","\\2",d2[,grep("TIMESTAMP",colnames(d2))[1]]))
  d2$DAY      <- as.numeric(gsub("^([0-9]{6})([0-9]{2})([0-9]{4})$","\\2",d2[,grep("TIMESTAMP",colnames(d2))[1]]))
  d2$HOUR     <- as.numeric(gsub("^([0-9]{8})([0-9]{2})([0-9]{2})$","\\2",d2[,grep("TIMESTAMP",colnames(d2))[1]]))
  d2$HRMINUTE <- d2$HOUR+(as.numeric(gsub("^([0-9]{10})([0-9]{2})$","\\2",d2[,grep("TIMESTAMP",colnames(d2))[1]]))/60)
  d2$DOY      <- as.numeric(format(strptime(as.character( d2[,grep("TIMESTAMP",colnames(d2))[1]] ) ,"%Y%m%d%H%M",tz=TZ),'%j'))
  d2$DOYdec   <- as.numeric(format(strptime(as.character( d2[,grep("TIMESTAMP",colnames(d2))[1]] ) ,"%Y%m%d%H%M",tz=TZ),'%j')) + (d2$HOUR)/24
  ## The TZ is not necessary but also does not matter here...
  
  sink(sprintf("%s/logs_aggav/General_Log_site_%s.txt",dirlarge,siteid))
  
  simple.heading(sprintf("%s",siteid))
  simple.heading("variable inspection...")
  
  ## which variables are now entirely NA?
  NAS <- apply(d2, 2, function(x) all(is.na(x)) )
  simple.heading("variables that do exist but are NA...")
  print(names(NAS[NAS==TRUE]))
  
  simple.heading("Deleting entire NA columns...")
  d2[,NAS] <- NULL
   
  simple.heading("Filter SWin and SWout <0 set to NA >1400 set NA...")
  sw <- grep("^SW_",colnames(d2))                                        ## I can do this because also SW_(.+)_QC flags are either NA or not smaller than 0
  print(colnames(d2)[sw])
  print(head(d2[,sw][d2[,sw]<0]))
  d2[,sw][d2[,sw]<0] <- NA
  d2[,sw][d2[,sw]>1400] <- NA
  
  simple.heading("Filter LWin and LWout <0 set to NA >1400 set NA...")
  lw <- grep("^LW_",colnames(d2))
  print(colnames(d2)[lw])
  print(head(d2[,lw][d2[,lw]<0]))
  d2[,lw][d2[,lw]<0] <- NA
  d2[,lw][d2[,lw]>1400] <- NA
  
  simple.heading("Filter ALBEDO set between 0 and 1 (resp 0 and 100)...") ## Albedo: range: -9 to 6407
  print("there is no albedo in FLUXNET")
  
  simple.heading("Filter Rnet > 1400 NA ...")                     ## Rnet:     range: -774 to 1744
  rn<- grep("NETRAD",colnames(d2),ignore.case = TRUE)
  print(colnames(d2)[rn])
  print(summary(d2[,rn]))
  print(head(d2[,rn][d2[,rn]>1400]))
  d2[,rn][d2[,rn]>1400] <- NA
  
  simple.heading("Filter Tsurf/Tsoil < -100 NA...")               ## Tsurf:     range: -105 to 30
  ts <- grep("TS_F_MDS_1",colnames(d2),ignore.case = TRUE)
  print(colnames(d2)[ts])
  print(summary(d2[,ts]))
  print(head(d2[,ts][d2[,ts] < -100]))
  d2[,ts][d2[,ts] < -100] <- NA
  
  simple.heading("Filter Ws > 120 NA & Ws <0 NA...")              ## range: -3 to 352 (both min and max unreasonable..according to Beaufort  scale: >37.2m/s is a hurricane!!) PROMICE, GC-Net, Ameriflux, FLUXNET Ws=m s-1
  ws <- unique(c(grep("^WS$",colnames(d2),ignore.case = TRUE), 
                 grep("^WS_F$",colnames(d2),ignore.case = TRUE))) ## wikipedia: an automatic weather station on Barrow Island, Australia, registered a maximum wind gust of 113.3 m/s 
  print(colnames(d2)[ws])
  print(summary(d2[,ws]))
  print(head(d2[,ws][d2[,ws]<0]))
  d2[,ws][d2[,ws]<0] <- NA
  print(head(d2[,ws][d2[,ws]>120]))
  d2[,ws][d2[,ws]>120] <- NA
  
  simple.heading("summary...")
  print(summary(d2))
  
  sink()
  
  ###################################################
  ########## Aggregation Procedure ##################
  simple.heading(sprintf("%s",siteid))
  simple.heading("start aggregation procedure...")
  simple.heading("set up parameters...")
  
  #####################################    
  ####### Set-up Parameters
  
  ####################################
  ####### dataset minimum requirements
  
# 1) qualityfilter data 
#    only keep quality flag = 1; keep an identifier how many qc=1 vs qc=0 are available for which day!
#    qc1 <- 0     ## set accepted qualityflag (0 = I only take measured values, 1=good quality gapfill)
   
  simple.heading("Apply criteria...")
    
  ######### Apply aggDOY function
  ######### qc0 cannot be "NA" because then the criterion will be undecidable...
  ######### qc1 is set to 0 in FLUXNET aggregation....
 
  aggDOY(d2=d2,qc0=0,qc1=0,minmeasd=0.65,maxgapd=0.2,qcflag="^(.+)_QC$",siteid=siteid,yearsite=yearsite,latlon=coords,dirlarge=dirlarge)
    
    
          } # fdir i.e. fluxnet site
  
###############################################
## aggregate further and store the data
files  <- list.files(sprintf("%s/aggav",dirlarge),pattern="(.+)_ALLY_DOY_(.+)_(.+)_var(.+).csv")

## derive the sites from names stored in files
sites1 <- sort(unique(gsub("(.+)_ALLY_DOY_(.+)_([0-9]{4}to[0-9]{4})_(.+)_var(.+).csv","\\2",files)))
sites[which(!sites%in% sites1)]
length(sites[which(!sites%in% sites1)])

## store the daily resolution data
storeDagg(sites=sites1,dirlarge)
## aggregate to monthly resolution
aggMf(sites=sites1,minmeasM=20,maxgapM=10,dirlarge)
## store the monthly resolution data
storeMagg(sites=sites1,dirlarge)

######################################################################################################
### Notes
# ## general variable description: see fv or https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/fullset-data-product/
# ## HH = half hourly, DD= Daily, WW= weekly, YY= yearly
# ## QC HH : 0 = measured; 1 = good quality gapfill; 2 = medium; 3 = poor
# ## QC>HH : fraction between 0-1, indicating percentage of measured and good quality gapfill data
# 
# ## general file description: see  https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/
# #[PUBLISHER]_[SITEID]_[PROCESSING-PIPELINE]_[GROUPING]_[RESOLUTION]_[FIRST-LAST-YEARS]_[SITEVERSION-CODEVERSION].[#EXT]
# 
# ## temporal resolution on files
# # HH: Half-Hourly time steps
# # HR: Hourly time steps (NOTE: documentation for HH also applies to HR)
# # DD: Daily time steps
# # WW: Weekly time steps
# # MM: Monthly time steps
# # YY: Yearly time steps
# 
# ###############################################################
# ###
# ### ---> Timestamp convention
# # Data files in half-hourly, hourly, and weekly resolutions use start and end timestamps. 
# # Data files using daily, monthly, and yearly resolutions use a single timestamp. 
# # Below are examples of resolutions that will use a single TIMESTAMP variable for timekeeping, 
# # and resolutions requiring the use of both TIMESTAMP_START and TIMESTAMP_END (blank spaces added for legibility).
# #
# ### ---> Time zone convention 
# # Time is reported in local standard time (i.e.,without Daylight Saving Time). The timezone information (with respect to UTC time) is reported in the site metadata.
# #  
# ## Column ordering
# # For text file data representations (i.e., CSV formatted), the variable/column order is relevant. The order of columns will NOT be guaranteed to be the same for different files (e.g., different sites), even though they will be similar in many cases. This means that any data processing routines should rely on the variable label (which is always consistent) and not the order of occurrence of that variable in the file. Timestamps are the only exception and will always be the first variable(s)/column(s) of the data file.
# #
# ### ---> Missing data
# # Missing data values are indicated with -9999 (without decimal points) as a replacement value, 
# # independent of the cause for the missing value.
