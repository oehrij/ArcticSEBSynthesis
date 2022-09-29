#########################################################
#####
##### project:       ArcticSEBSynthesis
##### title:         aggFunctions
##### description:   Aggregate source data of surface energy budget (SEB) components from (half-)hourly to daily and monthly resolution and store them 
##### publication:   Oehri, J. et al. (2022). Vegetation Type is an Important Predictor of the Arctic Summer Land Surface Energy Budget. Nature Communications.
##### author:        Jacqueline Oehri
##### date:          28.09.2022  
##### comments:      These functions only work with data that is preprocessed in a specific way as described below
#########################################################

# ##########################################################
# ## clean space
# # rm(list=ls(all=TRUE))
# 
# ##########################################################
# ## libraries
# library(pascal)         #library(devtools); install_github("pascal-niklaus/pascal/pascal")
# library(pgeo)           #library(devtools); install_github("pascal-niklaus/pascal/pascal")
# library(RColorBrewer)
# library(tidyverse)
# library(magick)
# library(readxl)
# library(raster)
# library(ncdf4)
# library(lattice)
# library(chron)
# library(lutz)           #tz_lookup_coords(lat, lon, method = "accurate", warn = TRUE)
# library(doParallel)
# library(zoo)            #for the rollmean function
# 
# ##########################################################
# ## functions on which the below depend...
# ### colors
# col   <- brewer.pal(7,"Accent")
# crf   <- colorRampPalette(col)
# qccol <- c("darkblue","darkorange")
# crf1  <- colorRampPalette(qccol)
# 
# ### make one dataframe
# load_data <- function(path, pattern) { 
#   files <- dir(path, pattern = pattern, full.names = TRUE)
#   tables <- lapply(files, read.csv)
#   do.call(rbind, tables)
# }
# 
# ## add legend outside plot margins
# add_legend <- function(...) {
#   opar <- par(fig=c(0.02, 0.97, 0.02, 0.97), oma=c(0, 0, 0, 0), 
#               mar=c(0, 0, 0, 0), new=TRUE)
#   on.exit(par(opar))
#   plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
#   legend(...)
# }
# 
# ################################################################################
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
##
#' @name aggDOY
#' @title Aggregate data to daily resolution
#' @description Aggregate source data of surface energy budget (SEB) components from (half-)hourly to daily resolution 
#' @param dirlarge A character specifying the directory where aggregated data should be stored. The last folder name will be used for naming the files (e.g. FLUXNET in this example: "C:/Users/jacqu/data/FLUXNET")
#' @param siteid A character of the form "name" specifying siteid to be processed
#' @param yearsite A character of the form "YYYY-YYYY", specifying minimum and maximum year covered by the data at the site
#' @param latlon  A dataframe of the form dataframe(Lat=Y,Lon=X)
#' @param qc0 Numeric: The minimum quality flag value (i.e. highest quality) possible, default: 0.
#' @param qc1 Numeric: The maximum accepted quality flag value, default: 0. For FLUXNET and AmeriFlux: 0 = only take measured values, 1=good quality gapfill
#' @param minmeasd Numeric: define minimum percentage of measurements needed to be present per day, default: 0.65, i.e. 65%
#' @param maxgapd Numeric: define maximum percentage of measurements allowed to be missing per day, is calculated in amount of hours 0.2*24 = 4.8
#' @param qcflag Character: quality flag label structure, default = "^(.+)_QC$" (FLUXNET structure)
#' @return Daily aggregated data & controlplots are directly saved into the aggav folder in the specified "dirlarge" directory
#' @export
aggDOY <- function(d2,qc0=0,qc1=0,minmeasd=0.65,maxgapd=0.2,qcflag="^(.+)_QC$",siteid,yearsite,latlon,dirlarge) {
  
  ##########################################################
  ## create directories needed
  for( di in c("aggav")){
    if( ! dir.exists(sprintf("%s/%s",dirlarge,di)) ){ 
      dir.create(sprintf("%s/%s",dirlarge,di))
    }}
  
  ## structure of error flag (replace (.+) with the variable of interest)
  if (!exists("qcflag")) {qcflag="^(.+)_QC$"}
  ## flag for "measured" or "best quality data"
  if (!exists("qc0")) {qc0=0}
  ## flag above which no data will be accepted anymore...maybe needs to be changed to 0 for measured only??
  if (!exists("qc1")) {qc1=1}
  ## minimum measurements per day (fraction)
  if (!exists("minmeasd")) {minmeasd=0.65}
  ## maximum temporal gap per day (fraction of 24 hours)
  if (!exists("maxgapd")) {maxgapd=0.2}
  
  ## 2) make daily averages
  ## of all days where criteria are fullfilled: min nr of measurements available & no gap larger than maxgap
  ## Adjust t in case needed to proper temporal resolution
 
  ## which is the most common gap between HRMINUTE
  if( safen(names(sort(table(c((d2$HRMINUTE),NA)-c(NA,(d2$HRMINUTE))),decreasing=TRUE)[1])) == 1 ) {t <-"HR"; maxmeasD <- 24} else 
    if( safen(names(sort(table(c((d2$HRMINUTE),NA)-c(NA,(d2$HRMINUTE))),decreasing=TRUE)[1])) == 0.5 ) {t <-"HH"; maxmeasD <- 48} else 
      if(! safen(names(sort(table(c((d2$HRMINUTE),NA)-c(NA,(d2$HRMINUTE))),decreasing=TRUE)[1])) %in% c(0.5,1)) {t <- safen(names(sort(table(c((d2$HRMINUTE),NA)-c(NA,(d2$HRMINUTE))),decreasing=TRUE)[1])); maxmeasD <- 24/t }
  
  ## define minimum of measures: At moment 0.65*48 [1] 31.2 |  0.65*24 [1] 15.6
     minmeasD <- minmeasd*maxmeasD  
  ## define maxgap allowed to be present; is calculated in amount of hours 0.2*24 = 4.8
     maxgapD  <- maxgapd*24         
     
    simple.heading(sprintf("Aggregate  %s data to daily resolution...",t))
    simple.heading("Apply criteria...")
    simple.heading(sprintf("maxQC=%d, minmeasD=%0.02f, maxgapD=%0.02f...",qc1,minmeasD,maxgapD))
    

 #####################################    
 ####### parameters for overview plot
    
  timesY <- seq(range(d2$DOYdec)[1],range(d2$DOYdec)[2],length.out=5)
  timesD <- seq(range(d2$HRMINUTE)[1],range(d2$HRMINUTE)[2],length.out=5)
  
  Yrs <-  sort(unique(d2$YEAR))
  doys <- sort(unique(d2$DOY))
  
  simple.heading(sprintf("nr of years = %s", length(Yrs) ))
  
 ## exclude timestamp and Quality flags (structural pattern of quality flag, e.g. "(.+)_QC$")
  idx <- c()
  for(idxx in c("TIMESTAMP",sprintf("%s",qcflag),"YEAR","MONTH","HOUR","DAY","DOY","DOYdec","TimestampJ")) {
    idx <- c(idx,grep(idxx,colnames(d2)))
  }
  
  #####################################   
  ### set up cluster for paralleling work for each variable
  ### (time gains are not so strong if there are not soo many variables)  
  
  simple.heading("Setting up cluster...")
  
  clus <- makeCluster(3)
  registerDoParallel(clus)
  t0 <- Sys.time()
  
  simple.heading(sprintf("starting going through variables...."))
  print(sort(colnames(d2)[-idx]))
  
  foreach(x=iter(sort(colnames(d2)[-idx])),
            .packages=c("pgeo","pascal","raster"),
            .export=ls(.GlobalEnv)
             ) %dopar% {
     
  if ( nrow(d2[is.finite(d2[[x]]),]) > 0 ) { 
  
  sink(sprintf("%s/aggav/Log_%s_ALLY_DOY_%s_%s_%s_var%s.txt", dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid, yearsite,t,x))
  simple.heading(sprintf("QCcrit=%0.0f,maxmeasD=%0.02f,minmeasD=%0.02f,maxgapD=%0.02f",qc1,maxmeasD,minmeasD,maxgapD))
    
  simple.heading(sprintf("extract metadata %s...",x))
  
  ## assess if quality flags are there
  q <- gsub("^.|.$","",gsub("([(][.][+][])])",x,qcflag))
  QCcrit <- ifelse(q %in% colnames(d2),"hasQC","noQC")
  
  ## select only values with measurements QC flag = qc1
  if(QCcrit == "hasQC") {
     
    simple.heading(sprintf("var %s has QC info",x))
    d3 <- d2[,c("YEAR","MONTH","DAY", "HOUR","HRMINUTE","DOY","DOYdec", x, q)]
    
    pdf(sprintf("%s/aggav/%s_aggDOY_%s_%s_var%s_%s.pdf",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge),siteid,t,x,QCcrit))
    
    
    ##############################
    ## make a flag if there are x values with no q information
    XnaQC <- 0
    XnaQC <- nrow(d3[!is.finite(d3[[q]]) & is.finite(d3[[x]]),])
    
    par(mfrow=c(1,2))
    hist(d3[[q]],xlab=sprintf("%s",q),main=sprintf("before qc %s",yearsite))
    hist(d3[[x]],xlab=sprintf("%s",x),main=sprintf("before qc %s",yearsite),sub="mean+-sd")
    abline(v=mean(d3[[x]],na.rm=TRUE),col="red")
    abline(v=mean(d3[[x]],na.rm=TRUE)+sd(d3[[x]],na.rm=TRUE),col="red",lty=2)
    abline(v=mean(d3[[x]],na.rm=TRUE)-sd(d3[[x]],na.rm=TRUE),col="red",lty=2)
    legend("topright",legend=sprintf("%0.02f+-%0.03f",mean(d3[[x]],na.rm=TRUE),sd(d3[[x]],na.rm=TRUE)),bty="n",text.col="red")
    
    #############
    ## JO, 13.01.2021 consider displaying (and setting 0) points where QC flag is infinite/NA
    plot(d3[which(!is.na(d3[[q]]) & d3[[q]]==sort(unique(d3[[q]],na.rm=TRUE))[1] & is.finite(d3[[x]])),]$DOY, d3[which(!is.na(d3[[q]]) & d3[[q]]==sort(unique(d3[[q]],na.rm=TRUE))[1] & is.finite(d3[[x]])),][[x]],  col=adjustcolor(c("grey50"),alpha.f=0.2),cex=0.3,pch=20,xlab="DOY",ylab=x)
    
    if(XnaQC != 0 ) {
      points(d3[which(!is.finite(d3[[q]]) & is.finite(d3[[x]])),]$DOY, d3[which(!is.finite(d3[[q]]) & is.finite(d3[[x]])),][[x]], col=adjustcolor(c("black"),alpha.f=0.2),cex=0.2,pch=20,xlab="DOY",ylab=x)
      legend("bottomright",title="QC flag",legend="QC=NA",fill=c("black"),bty="n")
    }
    for(qc in sort(unique(d3[[q]],na.rm=TRUE))[-1]) {
    #plot(d3[d3[[q]]==qc,][[x]]~d3[d3[[q]]==qc,][["DOY"]], main=sprintf("QC=%0.0f",qc),pch=20,col=adjustcolor("grey20",alpha.f=0.5),xlab="DOY",ylab=x,cex=0.5)
    ## original data
    points(d3[which(!is.na(d3[[q]]) & d3[[q]]==qc & is.finite(d3[[x]])),]$DOY, d3[which(!is.na(d3[[q]]) & d3[[q]]==qc & is.finite(d3[[x]])),][[x]], col=adjustcolor(c("grey50","purple","green","yellow","red")[which(sort(unique(d3[[q]],na.rm=TRUE))==qc)],alpha.f=0.1),cex=0.2,pch=20,xlab="DOY",ylab=x)
    }
    legend("topright",title="QC flag",legend=sort(unique(d3[[q]],na.rm=TRUE)),fill=c("grey50","purple","green","yellow","red")[c(1:length(sort(unique(d3[[q]],na.rm=TRUE))))],bty="n")
    #############
     
    ##############################
    ## apply hard criteria: QC <=1 -> Actually, sometimes, it appears that QC=NA is like QC=0!! in case of FI-Let G_F_MDS, for example!
    ##############################
    ## there were cases where no QC0 but QCNA was present
    ## This should be carefully evaluated!! Im replacing here the QC>0 with 0, under the assumption it was labeled wrongly, but I should contact
    ## site scientists!!
    
    XnaQCmin <- 0
    
    if(min(d3[[q]],na.rm=TRUE)>qc0) {
      simple.heading(sprintf("ATTENTION: data seems to have no QC0 values...!"))
      simple.heading("...(only because no single QC labeled as 0)")
      if(nrow(d3[!is.finite(d3[[q]]),]) > 0) {
        XnaQCmin <-NA
        simple.heading(sprintf("ATTENTION: Replacing %s values with 0..!..",as.character(XnaQCmin)))
        d3[!is.finite(d3[[q]]),][[q]] <- 0
      } 
      if(nrow(d3[ d3[[q]] ==  min(d3[[q]],na.rm=TRUE),]) > 0) {
        XnaQCmin <-min(d3[[q]],na.rm=TRUE)
        simple.heading(sprintf("ATTENTION: Replacing %s values with 0..!..",as.character(XnaQCmin)))
        d3[ d3[[q]] ==  min(d3[[q]],na.rm=TRUE),][[q]] <- 0
      } 
    }
   
    ##############################
    ##############################
    ## exclude QC > qc1 !!! (Keep in mind to later carefully check XnaQCmin criterion!!!)
    ## really set x values NA that have a QC Flag higher than qc1
    ## JO 13.01.2021: in AON (and potentially FLUXNET) relevant: x values that have NA/INF 
    ##                q indication should not be affected by below code, correct? Correct..!
    
    simple.heading("apply QC criteria....")
    if(nrow(d3[d3[[q]]>qc1 & is.finite(d3[[q]]) & is.finite(d3[[x]]),]) >0 ) {
    simple.heading(sprintf("%d values with QC > %d set to NA", nrow(d3[d3[[q]]>qc1 & is.finite(d3[[q]]) & is.finite(d3[[x]]),]),qc1))
    d3[d3[[q]]>qc1 & is.finite(d3[[q]]) & is.finite(d3[[x]]) ,][[x]] <- NA
    } else { simple.heading(sprintf("%d values with QC > %d set to NA", nrow(d3[d3[[q]]>qc1 & is.finite(d3[[q]]) & is.finite(d3[[x]]),]),qc1))}
   
    ##############################
    ##############################
    
    ##############################
    ## make a logfile when XnaQCmin or XnaQC are strange!
    if (XnaQC != 0 | XnaQCmin!= 0) {
      sink(sprintf("%s/aggav/%s_DOY_%s_%s_%s_var%s_%s_XnaQC.txt", dirlarge, gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid,yearsite,t,x,QCcrit))
      simple.heading(sprintf("var %s has no 0 QC flags",x))
      simple.heading(sprintf("XnaQC = %s",as.character(XnaQC)))
      simple.heading(sprintf("XnaQCmin = %s",as.character(XnaQCmin)))
      sink()
    }
    
    ##############################
    ## check after qc application
    hist(d3[[x]],xlab=sprintf("%s",x),main=sprintf("after qc %s",yearsite),sub="mean+-sd")
    abline(v=mean(d3[[x]],na.rm=TRUE),col="red")
    abline(v=mean(d3[[x]],na.rm=TRUE)+sd(d3[[x]],na.rm=TRUE),col="red",lty=2)
    abline(v=mean(d3[[x]],na.rm=TRUE)-sd(d3[[x]],na.rm=TRUE),col="red",lty=2)
    legend("topright",legend=sprintf("%0.02f+-%0.03f",mean(d3[[x]],na.rm=TRUE),sd(d3[[x]],na.rm=TRUE)),bty="n",text.col="red")
    legend("bottomright",legend=sprintf("QCNA=%s QCmin=%s",as.character(XnaQC),as.character(XnaQCmin)),bty="n",text.col="red")
    
    ##############################
    ## apply for each year
    
  for(y in Yrs) {
    
    simple.heading(sprintf("Year=%s.....",y))
    dy <- d3[d3$YEAR==y,]
    
    if ( nrow(dy[is.finite(dy[[x]]),]) > 0 ) { 
      
    ## aggregate by DOY set NA all x values where qualityflag is higher than allowed
    d5 <- aggr(dy,c("DOY"),c("YEAR=unique(YEAR)","MONTH=unique(MONTH)","DAY=unique(DAY)",
                             sprintf("mn%s=mean((%s),na.rm=TRUE)",x,x),
                             sprintf("min%s=min((%s),na.rm=TRUE)",x,x),
                             sprintf("max%s=max((%s),na.rm=TRUE)",x,x), 
                             sprintf("sd%s=sd((%s),na.rm=TRUE)",x,x)
                             ))
    
    ## aggregate by DOY Quality Info
    d51 <- aggr(dy[is.finite(dy[[x]]),],c("DOY"),
                                        c(sprintf("nQC0=length(which((%s)==%s))",q,qc0),
                                          sprintf("nmeas=length(HRMINUTE)"),
                                          sprintf("maxgap=max(c((HRMINUTE),NA)-c(NA,(HRMINUTE)),na.rm=TRUE)")
                                          ))
    
    d5<- merge(d5,d51,by="DOY",all.x=TRUE)
     
    d5$DOY <- as.numeric(as.character(d5$DOY)) 
  
    ##############################
    ## apply soft criteria: set NA days where.... 
    ## set NA values with too few good quality values or with a gap that is too large
    simple.heading(sprintf("%0.0f values < %0.02f minmeas set to NA",nrow(d5[d5$nmeas<minmeasD,]),minmeasD))
    d5[d5$nmeas<minmeasD & !is.na(d5$nmeas),which(colnames(d5) %in% c(sprintf("mn%s",x),
                                                                      sprintf("min%s",x),
                                                                      sprintf("max%s",x),
                                                                      sprintf("sd%s",x))) ] <- NA
    # ## ## where gap of values is too large
    simple.heading(sprintf("%0.0f values > %0.02f h maxgap set to NA",nrow(d5[d5$maxgap>maxgapD,]),maxgapD))
    d5[d5$maxgap>maxgapD & !is.na(d5$nmeas), which(colnames(d5) %in% c(sprintf("mn%s",x),
                                                                      sprintf("min%s",x),
                                                                      sprintf("max%s",x),
                                                                      sprintf("sd%s",x))) ] <- NA
    #d5$QCcrit<- QCcrit 
    
    
    # change Inf values to NA where no values where there with min/max function
    d5[is.na(d5[[sprintf("mn%s",x)]]), which(colnames(d5) %in% c(sprintf("min%s",x),
                                                                 sprintf("max%s",x), 
                                                                 sprintf("sd%s",x))) ] <- NA
    
    # sort dataframe and apply fraction of measurements
    d5 <- d5[order(d5$DOY),]
    d5$fracQC0  <- d5$nQC0/maxmeasD
    d5$fracmeas <- d5$nmeas/maxmeasD
    
    ## in case QC1 = 0: 
    #all(d5$fracQC0 == d5$fracmeas, na.rm=TRUE) [1] TRUE
    #all(d5$nQC0 == d5$nmeas, na.rm=TRUE)[1] TRUE
   
    par(mfrow=c(1,1))
    plot(NA,NA,xlim=c(min(timesY),max(timesY)),xaxt="n",xlab="DOY",ylab=x,ylim=c(range(d2[[x]],na.rm=TRUE)),main=sprintf("%s, %s, %s",siteid, y, t),
         sub="lines: mean, min and max")
    axis(1,at=timesY, labels=timesY, cex.axis=0.7)
    
    if (nrow(d5[is.finite(d5[[sprintf("mn%s",x)]]),]) > 0) {
      simple.heading(sprintf("nr of observation mean for year %s = %d", y, nrow(d5[is.finite(d5[[sprintf("mn%s",x)]]),]) ))
      
      ## original data
      
      ## JO 13.01.2021 plot also data where QC=NA
      if(nrow(d2[!is.finite(d2[[q]]) & is.finite(d2[[x]]),])>0) {
      points(d2[which(d2$YEAR==y & !is.finite(d2[[q]]) & is.finite(d2[[x]])),]$DOY, d2[which(d2$YEAR==y & !is.finite(d2[[q]]) & is.finite(d2[[x]])),][[x]],  col=adjustcolor("black",alpha.f=0.2),cex=0.65,pch=20)
      legend("bottomright",title="QC flag",legend="QC=NA",fill=c("black"),bty="n")
      }
      
      points(d2[which(d2$YEAR==y & !is.na(d2[[q]]) & d2[[q]]==qc1 & is.finite(d2[[x]])),]$DOY, d2[which(d2$YEAR==y & !is.na(d2[[q]]) & d2[[q]]==qc1 & is.finite(d2[[x]])),][[x]],  col=adjustcolor("grey50",alpha.f=0.2),cex=0.65,pch=20)
      points(d2[which(d2$YEAR==y & !is.na(d2[[q]]) & d2[[q]]==1 & is.finite(d2[[x]])),]$DOY, d2[which(d2$YEAR==y & !is.na(d2[[q]]) & d2[[q]]==1 & is.finite(d2[[x]])),][[x]],  col=adjustcolor("purple",alpha.f=0.2),cex=0.65,pch=20)
      points(d2[which(d2$YEAR==y & !is.na(d2[[q]]) & d2[[q]]==2 & is.finite(d2[[x]])),]$DOY, d2[which(d2$YEAR==y & !is.na(d2[[q]]) & d2[[q]]==2 & is.finite(d2[[x]])),][[x]],  col=adjustcolor("green",alpha.f=0.2),cex=0.65,pch=20)
      
      ## aggregated data
      lines(d5$DOY, d5[[sprintf("mn%s",x)]],  col="grey30")
      lines(d5$DOY, d5[[sprintf("min%s",x)]], col="grey30",lty=2)
      lines(d5$DOY, d5[[sprintf("max%s",x)]], col="grey30",lty=2)
      points(d5$DOY, d5[[sprintf("mn%s",x)]], col=adjustcolor(crf1(101)[1+(d5$fracmeas*100)],alpha.f=1),pch=20,cex=1)
      abline(v=d5[d5$DAY==1,]$DOY,lty=3,col="grey30")
      text(d5[d5$DAY==5,]$DOY, range(d2[[x]],na.rm=TRUE)[1]+0.15*range(d2[[x]],na.rm=TRUE)[2],paste("M",d5[d5$DAY==5,]$MONTH,sep="="),cex=0.8)
      legend("topright",title="quality fraction",legend=c("fracmeas=0","fracmeas=100","origQC=0","origQC=1","origQC=2"),fill=c(qccol,"grey50","purple","green"),bty="n")
      
    } else {simple.heading(sprintf("DOY data from %s has <=0 observations",x))}
    
    
    d5$siteid <- siteid
    d5$LAT <- latlon[,1]
    d5$LONG <- latlon[,2]
    
    ## this is smarter: make sure it is always the same length of DOY
    ## JO 30.06.2020: I now still do it so that all years have the same amount of rows...
    d5 <- merge(data.frame(DOY=c(1:366),YEAR=y),d5,by=c("DOY","YEAR"),all=TRUE)
    d5 <- d5[order(d5$DOY),]
    write.csv(d5,sprintf("%s/aggav/%s_%s_DOY_%s_%s_%s_var%s_%s.csv", dirlarge, gsub("^(.+)[/](.+)$","\\2",dirlarge), y, siteid, yearsite,t,x,QCcrit),row.names = FALSE)
    
    } else {simple.heading(sprintf("year %d data from %s has <=0 observations",y,x))} # nrow
   
    }# Yrs
    
  dev.off()
    
    D <- load_data(sprintf("%s/aggav", dirlarge),pattern=sprintf("%s_(.+)_DOY_%s_%s_%s_var%s_%s.csv$",gsub("^(.+)[/](.+)$","\\2",dirlarge),siteid,yearsite,t,x,QCcrit))
    summary(D)
    D$XnaQC <- XnaQC 
    D$XnaQCmin <- XnaQCmin
    
    if(length(D)>0){
    minY <- min(D$YEAR,na.rm=TRUE)
    maxY <- max(D$YEAR,na.rm=TRUE)
    if(length(which(!is.finite(D[[sprintf("min%s",x)]]))) >0) { D[is.na(D$nmeas), c(sprintf("min%s",x), sprintf("max%s",x))] <- NA  }
    if(length(which(!is.finite(D[["maxgap"]]))) >0) {  D[!is.finite(D$maxgap),]$maxgap  <- NA  }
    } else if(length(D)<=0) {
    minY <- "NA"
    maxY <- "NA"
    }
    
    write.csv(D,sprintf("%s/aggav/%s_ALLY_DOY_%s_%s_%s_var%s_%s.csv", dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid, paste(minY,maxY,sep="to"),t,x,QCcrit),row.names = FALSE)
    
    simple.heading(sprintf("there were %0.03f x values where QC info was NA",XnaQC))
    
    print("are yearsite and years measured in data in accordance?")
    print(paste(yearsite,paste(minY,maxY,sep="-")))
    print((yearsite==paste(minY,maxY,sep="-")))
    
   
  } else if(QCcrit == "noQC") {
    
    simple.heading(sprintf("var %s NO QC info",x))
    d3 <- d2[,c("YEAR","MONTH","DAY", "HOUR","HRMINUTE","DOY","DOYdec", x)]
    
    ##############################
    ## apply hard criteria: QC <=0
    simple.heading(sprintf("NO values set to NA"))
    #d3[d3[,q]>qc1,][[x]] <- NA
    
    ##############################
    ## set up results data frame
    ## write each year as csv and rbind at end.
    
    pdf(sprintf("%s/aggav/%s_aggDOY_%s_%s_var%s_%s.pdf",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge),siteid,t,x,QCcrit))
    
    par(mfrow=c(1,2))
    hist(d3[[x]],xlab=sprintf("%s",x),main=sprintf("histogram %s",yearsite),sub="mean+-sd")
    abline(v=mean(d3[[x]],na.rm=TRUE),col="red")
    abline(v=mean(d3[[x]],na.rm=TRUE)+sd(d3[[x]],na.rm=TRUE),col="red",lty=2)
    abline(v=mean(d3[[x]],na.rm=TRUE)-sd(d3[[x]],na.rm=TRUE),col="red",lty=2)
    legend("topright",legend=sprintf("%0.02f+-%0.03f",mean(d3[[x]],na.rm=TRUE),sd(d3[[x]],na.rm=TRUE)),bty="n",text.col="red")
    
    
    for(y in Yrs) {
      
      simple.heading(sprintf("Year=%s.....",y))
      dy <- d3[d3$YEAR==y,]
      
      
      if ( nrow(dy[is.finite(dy[[x]]),]) > 0 ) { 
  
      ## aggregate by DOY set NA all x values where qualityflag is higher than allowed
      d5 <- aggr(dy,c("DOY"),c("YEAR=unique(YEAR)","MONTH=unique(MONTH)","DAY=unique(DAY)",
                               sprintf("mn%s=mean((%s),na.rm=TRUE)",x,x),
                               sprintf("min%s=min((%s),na.rm=TRUE)",x,x),
                               sprintf("max%s=max((%s),na.rm=TRUE)",x,x), 
                               sprintf("sd%s=sd((%s),na.rm=TRUE)",x,x)
                               ))
      
      ## aggregate by DOY Quality Info
      d51 <- aggr(dy[is.finite(dy[[x]]),],c("DOY"),
                  c(sprintf("nQC0=NA"),
                    sprintf("nmeas=length(HRMINUTE)"),
                    sprintf("maxgap=max(c((HRMINUTE),NA)-c(NA,(HRMINUTE)),na.rm=TRUE)")
                    ))
      
      d5<- merge(d5,d51,by="DOY",all.x=TRUE)
      
      d5$DOY <- as.numeric(as.character(d5$DOY)) 
      
      ##############################
      ## apply 2nd hard criteria: set NA days where.... 
      ## set NA values with too few good quality values or with a gap that is too large
      # simple.heading(sprintf("%0.0f values < %0.02f minmeas set to NA",nrow(d5[d5$nmeas<minmeasD,]),minmeasD))
      d5[d5$nmeas<minmeasD & !is.na(d5$nmeas), which(colnames(d5) %in% c(sprintf("mn%s",x),
                                                                         sprintf("min%s",x),
                                                                         sprintf("max%s",x),
                                                                         sprintf("sd%s",x))) ] <- NA
      # ## ## where gap of values is too large
      simple.heading(sprintf("%0.0f values > %0.02f h maxgap set to NA",nrow(d5[d5$maxgap>maxgapD,]),maxgapD))
      d5[d5$maxgap>maxgapD & !is.na(d5$nmeas),which(colnames(d5) %in% c(sprintf("mn%s",x),
                                                                        sprintf("min%s",x),
                                                                        sprintf("max%s",x),
                                                                        sprintf("sd%s",x))) ] <- NA
      #d5$QCcrit<- QCcrit 
      
      # change Inf values to NA where no values where there with min/max function
      d5[is.na(d5[[sprintf("mn%s",x)]]), which(colnames(d5) %in% c(sprintf("min%s",x),
                                                                   sprintf("max%s",x), 
                                                                   sprintf("sd%s",x))) ] <- NA
      
      # sort dataframe and apply fraction of measurements
      d5 <- d5[order(d5$DOY),]
      d5$fracQC0  <- d5$nQC0/maxmeasD
      d5$fracmeas <- d5$nmeas/maxmeasD
      
      par(mfrow=c(1,1))
      plot(NA,NA,xlim=c(min(timesY),max(timesY)),xaxt="n",xlab="DOY",ylab=x,ylim=c(range(d2[[x]],na.rm=TRUE)),main=sprintf("%s, %s, %s, %s",siteid, y, t, QCcrit),
           sub="lines: mean, min and max")
      axis(1,at=timesY, labels=timesY, cex.axis=0.7)
      
      if (nrow(d5[is.finite(d5[[sprintf("mn%s",x)]]),]) > 0) {
        simple.heading(sprintf("nr of observation mean for year %s = %d", y, nrow(d5[is.finite(d5[[sprintf("mn%s",x)]]),]) ))
        ## original data
        points(d2[d2$YEAR==y,]$DOY, d2[d2$YEAR==y,][[sprintf("%s",x)]],  col=adjustcolor("grey50",alpha.f=0.2),cex=0.65,pch=20)
        
        ## aggregated data
        lines(d5$DOY, d5[[sprintf("mn%s",x)]],  col="grey30")
        lines(d5$DOY, d5[[sprintf("min%s",x)]], col="grey30",lty=2)
        lines(d5$DOY, d5[[sprintf("max%s",x)]], col="grey30",lty=2)
        points(d5$DOY, d5[[sprintf("mn%s",x)]], col=adjustcolor(crf1(101)[1+(d5$fracmeas*100)],alpha.f=1),pch=20,cex=1)
        abline(v=d5[d5$DAY==1,]$DOY,lty=3,col="grey30")
        text(d5[d5$DAY==5,]$DOY, range(d2[[x]],na.rm=TRUE)[1]+0.15*range(d2[[x]],na.rm=TRUE)[2],paste("M",d5[d5$DAY==5,]$MONTH,sep="="),cex=0.8)
        legend("topright",title="quality fraction",legend=c("fracmeas=0","fracmeas=100"),fill=qccol,bty="n")
        
      } else {simple.heading(sprintf("DOY data from %s has <=0 observations",x))}
      
      
      d5$siteid <- siteid
      d5$LAT <- latlon[,1]
      d5$LONG <- latlon[,2]
      
      ## make sure it is always the same length of DOY
      ## JO 30.06.2020: I now still do it so that all years have the same amount of rows...
      d5 <- merge(data.frame(DOY=c(1:366),YEAR=y),d5,by=c("DOY","YEAR"),all=TRUE)
      d5 <- d5[order(d5$DOY),]
      
      write.csv(d5,sprintf("%s/aggav/%s_%s_DOY_%s_%s_%s_var%s_%s.csv", dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge), y, siteid, yearsite,t,x,QCcrit),row.names = FALSE)
   
      } else {simple.heading(sprintf("year %d data from %s has <=0 observations",y,x))} # nrow
      
       } # Yrs
    
    dev.off()
    
    
    D <- load_data(sprintf("%s/aggav", dirlarge),pattern=sprintf("%s_(.+)_DOY_%s_%s_%s_var%s_%s.csv$",gsub("^(.+)[/](.+)$","\\2",dirlarge),siteid,yearsite,t,x,QCcrit))
    summary(D)
    ## add the same variable as for quality info x variables
    D$XnaQC <- -999 
    D$XnaQCmin <- -999
    
    if(length(D)>0){
    minY <- min(D$YEAR,na.rm=TRUE)
    maxY <- max(D$YEAR,na.rm=TRUE)
    if(length(which(!is.finite(D[[sprintf("min%s",x)]]))) >0) { D[is.na(D$nmeas), c(sprintf("min%s",x), sprintf("max%s",x))] <- NA  }
    if(length(which(!is.finite(D[["maxgap"]]))) >0) {  D[!is.finite(D$maxgap),]$maxgap  <- NA  }
    } else if(length(D)<=0) {
      minY <- "NA"
      maxY <- "NA"
    }
    
    write.csv(D,sprintf("%s/aggav/%s_ALLY_DOY_%s_%s_%s_var%s_%s.csv", dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid, paste(minY,maxY,sep="to"),t,x,QCcrit),row.names = FALSE)
    print("are yearsite and years measured in data in accordance?")
    print(paste(yearsite,paste(minY,maxY,sep="-")))
    print((yearsite==paste(minY,maxY,sep="-")))
    
          } # if no QCcrit

  print(sprintf("Variable=%s, finished at %s",x,Sys.time()))
  
  sink()
  
  } else {
    
    sink(sprintf("%s/aggav/Log_%s_ALLY_DOY_%s_%s_%s_var%s.txt", dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid,yearsite,t,x))
    simple.heading(sprintf("QCcrit=%0.0f,maxmeasD=%0.02f,minmeasD=%0.02f,maxgapD=%0.02f",qc1,maxmeasD,minmeasD,maxgapD))
    print(sprintf("Variable=%s, finished at %s",x,Sys.time()))
    print("NO data for this variable in all years")
    sink()
    
  }  # nrow finite x values
  
              } # x variables in dataset
  
  simple.heading(sprintf("all variables %s processed going to next site...",siteid))
    
    t1 <- Sys.time()
    print(t1-t0)
    stopCluster(clus)
    
    
          } #aggDOY function

################################################################################
##
#' @name aggM
#' @title Aggregate data to monthly resolution
#' @description Aggregate daily measurement produced by the aggDOY function data to monthly resolution
#' @param dirlarge A character specifying the directory that was used for the aggDOY function. The last folder name will be used for naming the files (e.g. FLUXNET in this example: "C:/Users/jacqu/data/FLUXNET")
#' @param sites A string of characters of the form c("name1","name2",...) specifying siteid's to be processed
#' @param minmeasM Numeric: minimum number of days of measurements per month
#' @param maxgapM  Numeric: maximum number of days without measurement allowed to be present
#' @return Monthly aggregated data & controlplots are directly saved into the aggav2 folder in the specified "dirlarge" directory
#' @export
#' 
aggM <- function(sites,minmeasM=20,maxgapM=10,dirlarge) { 
  
  ##########################################################
  ## create directories needed
  for(di in c("aggav2")){
    if( ! dir.exists(sprintf("%s/%s",dirlarge,di)) ){ 
      dir.create(sprintf("%s/%s",dirlarge,di))
    }}
  
  if (!exists("minmeasM")) {minmeasM=20}
  if (!exists("maxgapM")) {maxgapM=10}
  
  simple.heading("Apply criteria for monthly aggregation...")
  simple.heading(sprintf("minmeasM=%0.02f, maxgapD=%0.02f days...",minmeasM,maxgapM))
  
  ##########################################
  ## go through different sites (in different fdirs)
  
  for (siteid in sites) {
    simple.heading(sprintf("site = %s",siteid))
    files <- list.files(sprintf("%s/aggav", dirlarge),pattern=sprintf("^%s_ALLY_DOY_%s_(.+)", gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid))
    minYY <- min(as.numeric(gsub("^(.+)[_]([0-9]{4})to([0-9]{4})[_](.+)$","\\2", files)),na.rm=TRUE)
    maxYY <- min(as.numeric(gsub("^(.+)[_]([0-9]{4})to([0-9]{4})[_](.+)$","\\3", files)),na.rm=TRUE)
    
    pdf(sprintf("%s/aggav2/%s_aggM_%s_%s_varALL.pdf",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid,paste(minYY,maxYY,sep="-")))
    sink(sprintf("%s/aggav2/Log_%s_aggM_%s_%s_varALL.txt",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid,paste(minYY,maxYY,sep="-")))
    
    simple.heading(sprintf("       "))
    simple.heading(sprintf("site = %s, years =%d-%d",siteid,minYY,maxYY))
    simple.heading(sprintf("       "))
    
    ##########################################  
    ## aggregate each variable  (each file is a variable)
    
    for(f in files) {
      d <- read.csv(sprintf("%s/aggav/%s", dirlarge,f))
      d$Time <- d$YEAR + d$DOY/366
      
      x <- gsub("^(.+)var(.+)[_](.+)$","\\2",f)
      simple.heading(sprintf("var = %s",x))  
      
      minY <- min(d$YEAR,na.rm=TRUE)
      maxY <- max(d$YEAR,na.rm=TRUE)
      Yrs <- c(minY:maxY)
      print(Yrs)
      
      QCcrit <- gsub("^(.+)var(.+)[_](.+)[.]csv$","\\3",f)
      XnaQC <- mean(d$XnaQC,na.rm=TRUE)
      XnaQCmin <- mean(d$XnaQCmin,na.rm=TRUE)
      
      if ( nrow(d[is.finite(d[[sprintf("mn%s",x)]]),]) > 0 ) { 
        
        for( y in Yrs ) {
          simple.heading(sprintf("year = %d",y))  
          
          dy <- d[d$YEAR==y & is.finite(d$MONTH),]  
          
          ## I have to add the 366th doy again because I select only finite Year values above 
          ## (this is not necessary for months!)
          
          if(nrow(dy) > 0 ) {
            simple.heading(sprintf("nr rows dataframe = %d...",nrow(dy)))
            
            if (nrow(dy[is.finite(dy[[sprintf("mn%s",x)]]),]) > 0) {
              
              ##########################################
              ####### Timeseries plot
              par(mfrow=c(1,1))
              plot(NA,NA,xlim=c(min(dy$Time,na.rm=TRUE),max(dy$Time,na.rm=TRUE)),xaxt="n",xlab="YEARDOY",ylab=x,ylim=c(min(dy[,grep("^min",colnames(d))],na.rm=TRUE),max(dy[,grep("^max",colnames(d))],na.rm=TRUE)),main=sprintf("%s, %d",siteid,y),
                   sub="lines: mean, min and max")
              axis(1,at=seq(min(dy$Time,na.rm=TRUE),max(dy$Time,na.rm=TRUE),length.out=5), labels=as.character(round(seq(min(dy$Time,na.rm=TRUE),max(dy$Time,na.rm=TRUE),length.out=5),3)), cex.axis=0.7)
              
              simple.heading(sprintf("nr of observations = %d", nrow(dy[is.finite(dy[[sprintf("mn%s",x)]]),]) ))
              lines(dy$Time, dy[[sprintf("mn%s",x)]],  col="grey30")
              lines(dy$Time, dy[[sprintf("min%s",x)]], col="grey30",lty=2)
              lines(dy$Time, dy[[sprintf("max%s",x)]], col="grey30",lty=2)
              points(dy$Time, dy[[sprintf("mn%s",x)]], col=adjustcolor("black",alpha.f=0.6),pch=20,cex=0.8)
              points(dy$Time, dy[[sprintf("mn%s",x)]], col=adjustcolor(crf1(101)[1+(dy$fracmeas*100)],alpha.f=1),pch=20,cex=1)
              abline(v=dy[dy$DAY==1,]$Time,lty=3,col="grey30")
              
              text(dy[dy$DAY==5,]$Time, 
                   c(min(dy[,grep("min",colnames(d))],na.rm=TRUE)+0.15*max(dy[,grep("max",colnames(d))],na.rm=TRUE)),paste("M",dy[dy$DAY==5,]$MONTH,sep="="),
                   cex=0.8, col=adjustcolor("grey20",alpha.f=0.8))
              legend("topright",title="quality fraction",legend=c("fracmeas=0","fracmeas=100"),fill=qccol,bty="n")
              
            } else if (nrow(dy[is.finite(dy[[sprintf("mn%s",x)]]),]) <= 0) {
              
              par(mfrow=c(1,1))
              plot(NA,NA,xlim=c(min(dy$Time),max(dy$Time)),xaxt="n",xlab="YEARDOY",ylab=x,ylim=c(0,1),main=sprintf("%s, %d",siteid,y),
                   sub="lines: mean, min and max")
              axis(1,at=seq(min(dy$Time),max(dy$Time),length.out=5), labels=as.character(round(seq(min(dy$Time),max(dy$Time),length.out=5),3)), cex.axis=0.7)
              
            }
            
            # plot
            
            ##########################################
            ## aggregate by MONTH set NA all x values where qualityflag is higher than allowed
            d5 <- aggr(dy,c("MONTH"),c("YEAR=unique(YEAR)", "DOYstart=min((DOY),na.rm=TRUE)",
                                       sprintf("mnmn%s=mean((mn%s),na.rm=TRUE)",x,x),
                                       sprintf("mnmin%s=mean((min%s),na.rm=TRUE)",x,x),
                                       sprintf("mnmax%s=mean((max%s),na.rm=TRUE)",x,x), 
                                       sprintf("sdmn%s=sd((mn%s),na.rm=TRUE)",x,x)
            ))
            
            ## aggregate by DOY Quality Info (minmeas before was 0.65, therefore wherever more than 0.65 QC=0 available its like "all QC=1")
            ## build in "-1" as number in maxgap to deal with cases where no non missing values are there
            d51 <- aggr(dy[is.finite(dy[[sprintf("mn%s",x)]]),],c("MONTH"),
                        c(sprintf("nQC0=length(which((fracmeas)>0.65))"),
                          sprintf("nmeas=length(Time)"),
                          sprintf("maxgap=max(c(c((Time),NA)-c(NA,(Time)),-1),na.rm=TRUE)")))
            
            d5<- merge(d5,d51,by="MONTH",all.x=TRUE)
            
            d5$MONTH <- as.numeric(as.character(d5$MONTH)) 
            d5 <- d5[order(d5$MONTH),]
            
            
            if(any(!is.na(d5$maxgap))){ 
              d5$maxgap <- d5$maxgap*366
            } else if(all(is.na(d5$maxgap))){ 
              d5$nQC0 <- unlist(d5$nQC0)
              d5$nmeas <- unlist(d5$nmeas)
              d5$maxgap <- unlist(d5$maxgap)
            }
            
            ##############################
            ## apply soft criteria: set NA days where.... 
            ## set NA values with too few good quality values or with a gap that is too large
            simple.heading(sprintf("%0.0f values < %0.02f minmeas set to NA",nrow(d5[d5$nmeas<minmeasM,]),minmeasM))
            d5[d5$nmeas<minmeasM & !is.na(d5$nmeas),which(colnames(d5) %in% c(sprintf("mnmn%s",x),
                                                                              sprintf("mnmin%s",x),
                                                                              sprintf("mnmax%s",x), 
                                                                              sprintf("sdmn%s",x))) ] <- NA
            ## ## where gap of values is too large
            simple.heading(sprintf("%0.0f values > %0.02f h maxgap set to NA",nrow(d5[d5$maxgap>maxgapM,]),maxgapM))
            d5[(abs(d5$maxgap)>maxgapM | d5$maxgap<0) & !is.na(d5$nmeas), which(colnames(d5) %in% c(sprintf("mnmn%s",x),
                                                                                                    sprintf("mnmin%s",x),
                                                                                                    sprintf("mnmax%s",x), 
                                                                                                    sprintf("sdmn%s",x))) ] <- NA 
            
            
            ## set also NA the values where nmeas = NA
            d5[is.na(d5$nmeas), which(colnames(d5) %in% c(sprintf("mnmn%s",x),
                                                          sprintf("mnmin%s",x),
                                                          sprintf("mnmax%s",x), 
                                                          sprintf("sdmn%s",x)))] <- NA 
            
            
            
            ############
            ############ Add monthly averages to plot!
            d5$Time <- d5$YEAR+(d5$DOYstart+15)/366
            
            if (nrow(d5[is.finite(d5[[sprintf("mnmn%s",x)]]),]) > 0) {
              simple.heading(sprintf("nr of observations = %d", nrow(d5[is.finite(d5[[sprintf("mnmn%s",x)]]),]) ))
              lines(d5$Time, d5[[sprintf("mnmn%s",x)]],  col="red",)
              lines(d5$Time, d5[[sprintf("mnmin%s",x)]], col="red",lty=2)
              lines(d5$Time, d5[[sprintf("mnmax%s",x)]], col="red",lty=2)
              points(d5$Time, d5[[sprintf("mnmn%s",x)]], col="darkred",pch=20,cex=1.5)
              points(d5$Time, d5[[sprintf("mnmn%s",x)]], col="red",pch=20,cex=1)
              arrows(d5$Time, d5[[sprintf("mnmin%s",x)]], d5$Time, d5[[sprintf("mnmax%s",x)]],col="red",length=0,lwd=1.5,cex=1.5)
              legend("topleft",title="monthly averages",legend=c(sprintf("minmeas=%0.00f, maxgap=%0.0f days",minmeasM,maxgapM)),fill="red",bty="n",cex=0.75)
              
            } # plot
            
            
            ############
            ############ write to csv
            d5 <- merge(data.frame(MONTH=c(1:12),YEAR=y),d5,by=c("MONTH","YEAR"),all=TRUE)
            d5 <- d5[order(d5$MONTH),]
            write.csv(d5, sprintf("%s/aggav2/%s_aggM_%s_%d_var%s_%s.csv",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge),siteid,y,x,QCcrit),row.names = FALSE)
            
            
          } else if( nrow(dy) <= 0) {
            simple.heading(sprintf("strange DATA for Variable %s Year = %d",x,y))
            simple.heading(sprintf("nr rows dataframe = %d...",nrow(dy)))
            
          } 
          # NO data for year
          
          
        } ## for years    
        
        
        D <- load_data(sprintf("%s/aggav2", dirlarge),pattern=sprintf("^%s_aggM_%s_([0-9]{4})_var%s_%s.csv$", gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid,x,QCcrit))
        summary(D)
        D$mnXnaQC <- XnaQC
        D$mnXnaQCmin <- XnaQCmin
        
        write.csv(D, sprintf("%s/aggav2/%s_aggM_%s_AllY_%s_var%s_%s.csv",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge),siteid,paste(minY,maxY,sep="-"),x,QCcrit), row.names = FALSE)
        
        
        ##########################################
        ####### Histogram plot comparing daily and monthly data
        par(mfrow=c(1,3))
        hist(d[[sprintf("mn%s",x)]],xlab=sprintf("%s",x),main=sprintf("mean %s",x),sub=sprintf("%s mean+-sd",siteid))
        abline(v=mean(d[[sprintf("mn%s",x)]],na.rm=TRUE),col="red")
        abline(v=mean(d[[sprintf("mn%s",x)]],na.rm=TRUE)+sd(d[[sprintf("mn%s",x)]],na.rm=TRUE),col="red",lty=2)
        abline(v=mean(d[[sprintf("mn%s",x)]],na.rm=TRUE)-sd(d[[sprintf("mn%s",x)]],na.rm=TRUE),col="red",lty=2)
        
        abline(v=mean(D[[sprintf("mnmn%s",x)]],na.rm=TRUE),col=adjustcolor("purple",alpha.f=0.6))
        abline(v=mean(D[[sprintf("mnmn%s",x)]],na.rm=TRUE)+sd(D[[sprintf("mnmn%s",x)]],na.rm=TRUE),col=adjustcolor("purple",alpha.f=0.6),lty=2)
        abline(v=mean(D[[sprintf("mnmn%s",x)]],na.rm=TRUE)-sd(D[[sprintf("mnmn%s",x)]],na.rm=TRUE),col=adjustcolor("purple",alpha.f=0.6),lty=2)
        
        legend("topright",legend=c(sprintf("d: %0.02f+-%0.03f",mean(d[[sprintf("mn%s",x)]],na.rm=TRUE),sd(d[[sprintf("mn%s",x)]],na.rm=TRUE)),
                                   sprintf("m: %0.02f+-%0.03f",mean(D[[sprintf("mnmn%s",x)]],na.rm=TRUE),sd(D[[sprintf("mnmn%s",x)]],na.rm=TRUE))),bty="n",
               fill=c("red","purple"),title="mean (daily, monthly)")
        
        legend("bottomleft",legend=c(sprintf("XnaQC=%0.02f",XnaQC),sprintf("XnaQCmin=%0.02f",XnaQCmin)),bty="n",fill=c("green","green"),title="Variable QCflags")
        
        
        hist(d[[sprintf("min%s",x)]],xlab=sprintf("%s",x),main=sprintf("min %s",x),sub=sprintf("%s mean+-sd",siteid))
        abline(v=mean(d[[sprintf("min%s",x)]],na.rm=TRUE),col="red")
        abline(v=mean(d[[sprintf("min%s",x)]],na.rm=TRUE)+sd(d[[sprintf("min%s",x)]],na.rm=TRUE),col="red",lty=2)
        abline(v=mean(d[[sprintf("min%s",x)]],na.rm=TRUE)-sd(d[[sprintf("min%s",x)]],na.rm=TRUE),col="red",lty=2)
        
        abline(v=mean(D[[sprintf("mnmin%s",x)]],na.rm=TRUE),col=adjustcolor("purple",alpha.f=0.6))
        abline(v=mean(D[[sprintf("mnmin%s",x)]],na.rm=TRUE)+sd(D[[sprintf("mnmin%s",x)]],na.rm=TRUE),col=adjustcolor("purple",alpha.f=0.6),lty=2)
        abline(v=mean(D[[sprintf("mnmin%s",x)]],na.rm=TRUE)-sd(D[[sprintf("mnmin%s",x)]],na.rm=TRUE),col=adjustcolor("purple",alpha.f=0.6),lty=2)
        
        legend("topright",legend=c(sprintf("d: %0.02f+-%0.03f",mean(d[[sprintf("min%s",x)]],na.rm=TRUE),sd(d[[sprintf("min%s",x)]],na.rm=TRUE)),
                                   sprintf("m: %0.02f+-%0.03f",mean(D[[sprintf("mnmin%s",x)]],na.rm=TRUE),sd(D[[sprintf("mnmin%s",x)]],na.rm=TRUE))),bty="n",
               fill=c("red","purple"),title="min (daily, monthly)")
        
        
        hist(d[[sprintf("max%s",x)]],xlab=sprintf("%s",x),main=sprintf("max %s",x),sub=sprintf("%s mean+-sd",siteid))
        abline(v=mean(d[[sprintf("max%s",x)]],na.rm=TRUE),col="red")
        abline(v=mean(d[[sprintf("max%s",x)]],na.rm=TRUE)+sd(d[[sprintf("max%s",x)]],na.rm=TRUE),col="red",lty=2)
        abline(v=mean(d[[sprintf("max%s",x)]],na.rm=TRUE)-sd(d[[sprintf("max%s",x)]],na.rm=TRUE),col="red",lty=2)
        
        
        abline(v=mean(D[[sprintf("mnmax%s",x)]],na.rm=TRUE),col=adjustcolor("purple",alpha.f=0.6))
        abline(v=mean(D[[sprintf("mnmax%s",x)]],na.rm=TRUE)+sd(D[[sprintf("mnmax%s",x)]],na.rm=TRUE),col=adjustcolor("purple",alpha.f=0.6),lty=2)
        abline(v=mean(D[[sprintf("mnmax%s",x)]],na.rm=TRUE)-sd(D[[sprintf("mnmax%s",x)]],na.rm=TRUE),col=adjustcolor("purple",alpha.f=0.6),lty=2)
        
        legend("topright",legend=c(sprintf("d: %0.02f+-%0.03f",mean(d[[sprintf("max%s",x)]],na.rm=TRUE),sd(d[[sprintf("max%s",x)]],na.rm=TRUE)),
                                   sprintf("m: %0.02f+-%0.03f",mean(D[[sprintf("mnmax%s",x)]],na.rm=TRUE),sd(D[[sprintf("mnmax%s",x)]],na.rm=TRUE))),bty="n",
               fill=c("red","purple"),title="max (daily, monthly)")
        
        
        ##########################################
        ####### Timeseries plot
        par(mfrow=c(1,1))
        plot(NA,NA,xlim=c(min(d$Time,na.rm=TRUE),max(d$Time,na.rm=TRUE)),xaxt="n",xlab="YEARDOY",ylab=x,ylim=c(min(d[,grep("^min",colnames(d))],na.rm=TRUE),max(d[,grep("^max",colnames(d))],na.rm=TRUE)),main=sprintf("%s, %s",siteid,x),
             sub="lines: mean, min and max")
        axis(1,at=seq(min(d$Time,na.rm=TRUE),max(d$Time,na.rm=TRUE),length.out=5), labels=as.character(round(seq(min(d$Time,na.rm=TRUE),max(d$Time,na.rm=TRUE),length.out=5),3)), cex.axis=0.7)
        
        if (nrow(d[is.finite(d[[sprintf("mn%s",x)]]),]) > 0) {
          simple.heading(sprintf("nr of observations = %d", nrow(d[is.finite(d[[sprintf("mn%s",x)]]),]) ))
          lines(d$Time, d[[sprintf("mn%s",x)]],  col="grey70")
          lines(d$Time, d[[sprintf("min%s",x)]], col="grey70",lty=2)
          lines(d$Time, d[[sprintf("max%s",x)]], col="grey70",lty=2)
          points(d$Time, d[[sprintf("mn%s",x)]], col=adjustcolor(crf1(101)[1+(d$fracmeas*100)],alpha.f=0.6),pch=20,cex=1)
          abline(v=d[d$DAY==1,]$Time,lty=3,col="grey70")
          
          legend("topright",title="quality fraction",legend=c("fracmeas=0","fracmeas=100"),fill=qccol,bty="n")
        }
        
        lines(D$Time, D[[sprintf("mnmn%s",x)]],  col="purple",lwd=1.2)
        lines(D$Time, D[[sprintf("mnmin%s",x)]], col="purple",lty=2,lwd=1.2)
        lines(D$Time, D[[sprintf("mnmax%s",x)]], col="purple",lty=2,lwd=1.2)
        points(D$Time,D[[sprintf("mnmn%s",x)]], col="purple",pch=20,cex=1.2)
        points(D$Time, D[[sprintf("mnmn%s",x)]], col="pink",pch=20,cex=1)
        arrows(D$Time, D[[sprintf("mnmin%s",x)]], D$Time, D[[sprintf("mnmax%s",x)]],col="purple",length=0,cex=1.5,lwd=1.2)
        
        legend("topleft",legend=c(sprintf("d: %0.02f+-%0.03f",mean(d[[sprintf("mn%s",x)]],na.rm=TRUE),sd(d[[sprintf("mn%s",x)]],na.rm=TRUE)),
                                  sprintf("m: %0.02f+-%0.03f",mean(D[[sprintf("mnmn%s",x)]],na.rm=TRUE),sd(D[[sprintf("mnmn%s",x)]],na.rm=TRUE))),bty="n",
               fill=c("darkorange","purple"),title="mean (daily, monthly)")
        
        
      } else {simple.heading(sprintf("aggDOY data from %s has <=0 observations",x))}
      
    } # for files i.e. variables
    
    sink()
    dev.off()
    
    simple.heading("site finished going to next...")
    
  }# for sites 
  
  
} # aggM function end

################################################################################
##
#' @name storeDagg
#' @title Aggregation Procedure DAILY (depends on aggDOY function outputs in aggav folder)
#' @description Stores daily aggregated data from aggDOY function in the aggav4 folder
#' @param dirlarge A character specifying the directory that was used for the aggDOY function. The last folder name will be used for naming the files (e.g. FLUXNET in this example: "C:/Users/jacqu/data/FLUXNET")
#' @param sites A string of characters of the form c("name1","name2",...) specifying siteid's to be processed
#' @return Monthly aggregated data & controlplots are directly saved into the aggav2 folder in the specified "dirlarge" directory
#' @export
#' 
storeDagg <- function(sites,dirlarge) {
  
  ##########################################################
  ## create directories needed
  for( di in c("aggav4")){
    if( ! dir.exists(sprintf("%s/%s",dirlarge,di)) ){ 
      dir.create(sprintf("%s/%s",dirlarge,di))
    }}
  
  simple.heading("Daily aggregation...")
  
  for (siteid in sites) {
    simple.heading(sprintf("nr. %d site = %s",which(sites==siteid),siteid))
    files <- list.files(sprintf("%s/aggav", dirlarge),pattern=sprintf("^%s_ALLY_DOY_%s_(.+)", gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid))
    
    DTA    <-data.frame(Time=NA,DOY=NA,YEAR=NA,MONTH=NA,DAY=NA,
                        mean=NA,min=NA,max=NA,sd=NA,
                        maxgap=NA, fracQC0=NA, fracmeas=NA, XnaQC=NA, XnaQCmin=NA,
                        siteid=NA, LAT=NA, LONG=NA,
                        var=NA, QCcrit=NA) 
    
    for(f in files) {
      d <- read.csv(sprintf("%s/aggav/%s", dirlarge,f))
      d$Time <- d$YEAR + d$DOY/366
      
      x <- gsub("^(.+)var(.+)[_](.+)$","\\2",f)
      
      simple.heading(sprintf("var = %s",x))  
      simple.heading(sprintf("nrow data = %d",nrow(d)))  
      simple.heading(sprintf("Years = %d - %d",min(d$YEAR,na.rm=TRUE),max(d$YEAR,na.rm=TRUE)))  
      simple.heading(sprintf("avg nmeas if ~all DOYs 366 = %0.02f",length(unique(d$YEAR))*366))  
      
      minY <- min(d$YEAR,na.rm=TRUE)
      maxY <- max(d$YEAR,na.rm=TRUE)
      Yrs <- c(minY:maxY)
      print(Yrs)
      
      QCcrit <- gsub("^(.+)var(.+)[_](.+)[.]csv$","\\3",f)
      
      dta <- data.frame(d[,c("Time", "DOY", "YEAR", "MONTH", "DAY",
                             sprintf("mn%s",x),sprintf("min%s",x),sprintf("max%s",x),sprintf("sd%s",x),
                             "maxgap", "fracQC0", "fracmeas", "XnaQC", "XnaQCmin",
                             "siteid", "LAT","LONG")])
      dta$var <- x
      dta$QCcrit <- QCcrit
      
      idxx <- c()
      for(vs in  c(sprintf("mn%s",x),sprintf("min%s",x),sprintf("max%s",x),sprintf("sd%s",x))) {
        idxx <- c(idxx, which(colnames(dta)==vs))
      }
      colnames(dta)[idxx] <- c("mean","min","max","sd")
      
      DTA <- rbind(DTA,dta)
      
    } # files/i.e. variables
    
    DTA <- DTA[-1,]
    
    write.csv(DTA, sprintf("%s/aggav4/%s_aggD_%s_Allyv.csv",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge),siteid),row.names = FALSE)
    
  } # sites 
  
  
  D1 <- load_data(sprintf("%s/aggav4/",dirlarge),pattern=sprintf("^.+_aggD_.+Allyv.csv"))
  
  pdf(sprintf("%s/aggav4/%s_aggD_Allsyv_controlplot.pdf",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge)))
  
  for(vv1 in unique(D1$var)) {
    cat(vv1,"\n")
    d1 <- D1[D1$var==sprintf("%s",vv1),]
    
    for(yy in c("mean","min","max","sd")) {
      if (TRUE %in% is.finite(d1[[sprintf("%s",yy)]]) ) { 
        summary(d1[[sprintf("%s",yy)]])
        plot(d1[[sprintf("%s",yy)]]~d1$siteid,main=sprintf("%s",vv1),cex.axis=0.8,las=2,xlab="",ylab=sprintf("%s",yy))
        points(d1[[sprintf("%s",yy)]]~jitter(as.numeric(d1$siteid)),col=adjustcolor(crf(366)[d1$DOY],alpha.f=0.3),pch=20)
        plot(d1[[sprintf("%s",yy)]]~d1$siteid,main=sprintf("%s",vv1),cex.axis=0.8,las=2,xlab="",ylab=sprintf("%s",yy),add=TRUE,col=NA)
        legend("bottomright",bty="n",cex=0.8,title="months",legend=as.character(c(1:12)),fill=crf(12)[c(1:12)])
      } else { 
        plot(NA,NA,main=sprintf("%s",vv1),xlim=c(0,1),ylim=c(0,1),cex.axis=0.8,las=2,xlab="",ylab=sprintf("%s",yy))
      }
    } 
  }
  
  write.csv(D1, sprintf("%s/aggav4/%s_aggD_Allsyv.csv",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge)), row.names = FALSE)
  
  dev.off()
  
  print(summary(D1))
  print(nrow(D1))
  
} # storeDagg finish

################################################################################
##
#' @name storeMagg
#' @title Aggregation Procedure MONTHLY (depends on aggM function outputs in aggav2 folder)
#' @description Stores monthly aggregated data from aggM function in the aggav4 folder
#' @param dirlarge A character specifying the directory that was used for the aggDOY function. The last folder name will be used for naming the files (e.g. FLUXNET in this example: "C:/Users/jacqu/data/FLUXNET")
#' @param sites A string of characters of the form c("name1","name2",...) specifying siteid's to be processed
#' @return Monthly aggregated data & controlplots are directly saved into the aggav2 folder in the specified "dirlarge" directory
#' @export
#' 
storeMagg <- function(sites,dirlarge) { 
  
  simple.heading("Monthly aggregation...")
  
  #####################################    
  ####### load data and aggregate to one

  for (siteid in sites) {
    simple.heading(sprintf("nr. %d site = %s",which(sites==siteid),siteid))
    
    files <- list.files(sprintf("%s/aggav2", dirlarge),pattern=sprintf("^%s_aggM_%s_AllY(.+).csv$", gsub("^(.+)[/](.+)$","\\2",dirlarge), siteid))
    
    DTA    <-data.frame(Time=NA,DOYstart=NA,YEAR=NA,MONTH=NA,
                        mnmean=NA,mnmin=NA,mnmax=NA,sdmn=NA,
                        maxgap=NA,nQC0=NA, nmeas=NA,  mnXnaQC=NA, mnXnaQCmin=NA,
                        siteid=NA, LAT=NA, LONG=NA,
                        var=NA, QCcrit=NA)

    for(f in files) {
      d <- read.csv(sprintf("%s/aggav2/%s", dirlarge,f))
   
      x <- gsub("^(.+)var(.+)[_](.+)$","\\2",f)
      
      simple.heading(sprintf("var = %s",x))  
      simple.heading(sprintf("nrow data = %d",nrow(d)))  
      simple.heading(sprintf("Years = %d - %d",min(d$YEAR,na.rm=TRUE),max(d$YEAR,na.rm=TRUE)))  
      simple.heading(sprintf("avg nmeas if ~all Months present = %0.02f",length(unique(d$YEAR))*12))  
      
      minY <- min(d$YEAR,na.rm=TRUE)
      maxY <- max(d$YEAR,na.rm=TRUE)
      Yrs <- c(minY:maxY)
      print(Yrs)
      
      QCcrit <- gsub("^(.+)var(.+)[_](.+)[.]csv$","\\3",f)
      
      
      dta <- data.frame(d[,c("Time", "DOYstart", "YEAR", "MONTH",
                             sprintf("mnmn%s",x),sprintf("mnmin%s",x),sprintf("mnmax%s",x),sprintf("sdmn%s",x),
                             "maxgap","nQC0", "nmeas", "mnXnaQC", "mnXnaQCmin")])
      dta$siteid <- siteid
      dta$LAT <- NA
      dta$LONG <- NA
      dta$var <- x
      dta$QCcrit <- QCcrit
      
      idxx <- c()
      for(vs in  c(sprintf("mnmn%s",x),sprintf("mnmin%s",x),sprintf("mnmax%s",x),sprintf("sdmn%s",x))) {
        idxx <- c(idxx, which(colnames(dta)==vs))
      }
      colnames(dta)[idxx] <- c("mnmean","mnmin","mnmax","sdmn")
      
      DTA <- rbind(DTA,dta)
      
    } # files/i.e. variables
    
    DTA <- DTA[-1,]
    
    write.csv(DTA, sprintf("%s/aggav4/%s_aggM_%s_Allyv.csv",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge),siteid),row.names = FALSE)
    
  } # sites 
  
  
  D1 <- load_data(sprintf("%s/aggav4/",dirlarge),pattern=sprintf("^.+_aggM_.+Allyv.csv"))
  
  pdf(sprintf("%s/aggav4/%s_aggM_Allsyv_controlplot.pdf",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge)))
  
  for(vv1 in unique(D1$var)) {
    cat(vv1,"\n")
    d1 <- D1[D1$var==sprintf("%s",vv1),]
    
    for(yy in c("mnmean","mnmin","mnmax","sdmn")) {
      if (TRUE %in% is.finite(d1[[sprintf("%s",yy)]]) ) { 
        summary(d1[[sprintf("%s",yy)]])
        plot(d1[[sprintf("%s",yy)]]~d1$siteid,main=sprintf("%s",vv1),cex.axis=0.8,las=2,xlab="",ylab=sprintf("%s",yy))
        points(d1[[sprintf("%s",yy)]]~jitter(as.numeric(d1$siteid)),col=adjustcolor(crf(12)[d1$MONTH],alpha.f=0.3),pch=20)
        plot(d1[[sprintf("%s",yy)]]~d1$siteid,main=sprintf("%s",vv1),cex.axis=0.8,las=2,xlab="",ylab=sprintf("%s",yy),add=TRUE,col=NA)
        legend("bottomright",bty="n",cex=0.8,title="months",legend=as.character(c(1:12)),fill=crf(12)[c(1:12)])
      } else { 
        plot(NA,NA,main=sprintf("%s",vv1),xlim=c(0,1),ylim=c(0,1),cex.axis=0.8,las=2,xlab="",ylab=sprintf("%s",yy))
      }
    } 
  }
  
  write.csv(D1, sprintf("%s/aggav4/%s_aggM_Allsyv.csv",dirlarge,gsub("^(.+)[/](.+)$","\\2",dirlarge)), row.names = FALSE)
  
  dev.off()
  
  
  print(summary(D1))
  print(nrow(D1))
  
} # storeMagg finish

##################################### 

  
  

