############## Activity pattern analysis weasel Dieren ################

rm(list=ls())

library(overlap);library(activity);library(sp)

## read data
load("weasel-activity-data.RData")

activity.list <- list()

for(i in 1:3){
  temp.clock.time <- as.numeric(2*pi*weasel.list[[i]]$time)
  temp.dates <- weasel.list[[i]]$Datum.tijd
  temp.xy <- cbind(weasel.list[[i]]$x,weasel.list[[i]]$y)
  temp.xy.sp <- SpatialPoints(temp.xy,proj4string=CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs"))
  temp.xy.sp <- spTransform(temp.xy.sp,"+proj=longlat +datum=WGS84")
  temp.weasel <- sunTime(temp.clock.time,temp.dates,temp.xy.sp)
  activity.list[[i]] <- fitact(temp.weasel, reps=1000, sample="model", adj = 1.5)
}

## compare activity patterns
compare.17.18 <- compareCkern(activity.list[[1]],activity.list[[3]],reps=1000)
compare.sp.su <- compareCkern(activity.list[[2]],activity.list[[3]],reps=1000)
compare.17.18
compare.sp.su
