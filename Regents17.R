library(data.table)
library(reshape2)
setwd("/Users/katemarfleet/Documents/ZSL/Analysis/BasicCamtrapProcessing") #sets working directory
source("/Users/katemarfleet/Documents/ZSL/Analysis/BasicCamtrapProcessing/basic_camtrap_processing.R") 

infolder <- "/Users/katemarfleet/Documents/ZSL/Analysis/RegentsTagged"
outfile <- "/Users/katemarfleet/Documents/ZSL/Analysis/metadata.csv"
exifolder <- "/Users/katemarfleet/Documents/ZSL/Analysis/exiftool.exe"


##################################################################################################################################
#MAKE SITE AND EVENT DATAFRAMES
##################################################################################################################################
setwd("/Users/katemarfleet/Documents/ZSL/Analysis")
imgdat <- read.csv("metadata.csv")
#Read site data
sitedat <- read.csv("sitedat.csv")
sitedat$deploy_start <- strptime(sitedat$deploy_start, "%d/%m/%Y")
sitedat$deploy_end <- strptime(sitedat$deploy_end, "%d/%m/%Y")
sitedat$nights <- as.numeric(sitedat$deploy_end - sitedat$deploy_start)/(24*60^2)
View(sitedat)
##write.csv(imgdat, file = "imgdat.csv")

 
#Read tag data
tagdat <- get.tagdat(imgdat, sitedat)
View(tagdat)

#Make row-per-event data
eventdat <- get.eventdat(1, tagdat, sitedat)
View(eventdat)


#Tabulate site-by-species events
events <- event.table(eventdat)
trate <- events/sitedat$nights

#Overall species events and trap rates (per 100 nights)
ev_sp <- apply(events, 2, sum)
tr_sp <- 100 * ev_sp / sum(sitedat$nights)
data.frame(Events=ev_sp, Traprate=round(tr_sp,1))
range(sitedat$nights)
sum(sitedat$nights)
min(sitedat$deploy_start)
max(sitedat$deploy_end)

##MAPPING

##install.packages("XML")
##install.packages("jpeg")
##install.packages("geosphere")
##install.packages("SDMTools")

library(XML)
library(jpeg)
library(geosphere)
library(SDMTools)


source("/Users/katemarfleet/Documents/ZSL/Analysis/BasicCamtrapProcessing/mapping.r")
wd <- getwd() 
setwd("/Users/katemarfleet/Documents/ZSL/Analysis/Mapping")

regentsbase <- loadmap("regentsmap.jpg", "regentspath.kml")
regentsboundary <- getXMLcoords("studyarea.kml")
setwd(wd)


par(mar=c(2,2,2,2))
sp <- "hedgehog"
for(sp in colnames(events)){
  tr <- trate[,sp]
  plotmap(regentsbase)
  mtext(sp)
  addshape(regentsbase, regentsboundary, "poly", col=2)
  latlong <- sitedat[,c("lat","long")]
  names(latlong) <- c("long", "lat")
  addshape(regentsbase, latlong, "point", col=ifelse(tr==0, 4, 2), cex=0.5*tr+0.5, pch=16)
}
trrange <- range(trate[trate>0])
tr <- round(c(0,trrange[1], mean(trrange), trrange[2]), 2) 



plot.new()
legend(0.5,0.5, tr, col=c(4,2,2,2), pt.cex=0.5*tr+0.5, pch=16, 
       title="Records per night", y.intersp=1.5)    ##HAVE THESE 3 LINES WORKED?

install.packages("unmarked")
library(unmarked)
dmat <- get.dmatrix("hedgehog", eventdat, sitedat)
View(dmat)
dmat3 <- condense.matrix(dmat, 3)
View(dmat3$detection)
View(dmat3$effort)

##DONE UP TO HERE


weatherdat <- read.csv("/Users/katemarfleet/Documents/ZSL/Analysis/regentsweather.csv")
obscovs <- list(temp=matrix(weatherdat$temp.avg, nrow=nrow(dmat), ncol=nrow(weatherdat), byrow=T),
                precip=matrix(weatherdat$precip.mm.sum, nrow=nrow(dmat), ncol=nrow(weatherdat), byrow=T))

sitecovs <- sitedat[match(rownames(dmat), sitedat$site), c("group","habitat")]