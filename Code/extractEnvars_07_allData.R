#Extract all environmental data
#and save this as a raster stack for further use

#Mark Hindell, Ryan Reisinger & Ben Raymond

#last updated 2018-10-01

#Runs per species and stage
#------------------------------------------------

rm(list=ls(all=TRUE))

library(rgdal)
library(ncdf4)
library(raster)
library(raadtools)
library(dplyr)
library(sp)
library(maptools)
library(viridis)
library(orsifronts)

#memory.limit(32000)
#-----------------------------------------------

## Set working directory
setwd(paste0("~/RAATD_01/RAATD/envarsAll"))

#-----------------------------------------------

## Generate a list of weekly dates for 10 years (2004-2014)

dates <- data.frame(dates=seq(as.POSIXct("1991-11-25"), as.POSIXct("2016-06-20"), 86400*7))
stg_date <- dates$dates

#-----------------------------------------------

## Make an exmpty raster
lims <- c(-180, 180, -80, -40)
grd <- raster(extent(lims), resolution = c(0.1, 0.1))

## The same empty raster on 0-360, for some products
grd360 <- raster(extent(0, 360, -80, -40), resolution = c(0.1, 0.1))
#----------------------------------------------- 

#----------------------------------------------- 
## Create masks

## land mask
# bath <- readderivaadc("bathymetry", xylim = lims)
# bath <- resample(bath, grd, method = "bilinear")
# bath <- crop(bath, extent(grd))
# landmask <- bath
# landmask[landmask > -1] <- NA
# rm(bath)

## load mask already created (see folder 'antarcticMask')
## this mask includes the ice-shelves
landmask <- raster("../antarcticMask/mask.grd")

## ice mask
ice <- mean(readice(stg_date, setNA = T, latest = F), na.rm = T)
rd3dummy <- grd
projection(grd) <- "+proj=longlat +datum=WGS84"
rd3dummy <- grd
rd3dummy[] <- 0
rd3points <-  rasterToPoints(rd3dummy, spatial = TRUE)
mice <- extract(ice, rd3points, method = "bilinear")
mice <- setValues(grd, mice)
mice[is.na(mice)] <- 0
mice[mice > 99] <- NA
rm(rd3dummy, rd3points)

## Combine the masks
msk <- mask(mice, landmask)
msk[!is.na(msk)] <- 0

## Clean up
rm(landmask, ice, mice)

#----------------------------------------------- 
## Function to check NAs against expectation
naProp <- function(rst, msk, threshold = 10) {
  expectedNAs <- length(which(is.na(values(msk))))
  isNAs <- length(which(is.na(values(rst))))
  prp <- (isNAs - expectedNAs) / (length(values(rst)) - expectedNAs)*100
  if (prp < threshold) {
    giveback <- rst
    print(paste0("Okay - NAs below threshold: ", round(prp, 2), "%"))
  } else {
    giveback <- NULL
    print(paste0("Set to NULL! NAs exceed threshold: ", round(prp, 2), "%"))
  }
  return(giveback)
}
  
  
#----------------------------------------------- 

## Get the Envars

#######################
## SST
sst <- readsst(stg_date, xylim=extend(extent(grd), 5), lon180 = TRUE, time.resolution = "monthly", latest = FALSE) #uses the extent of the original, un-projected grid, defined above
## BR extend SST manually in longitude, because the above does not work
set_extent <- function(z, ext) {extent(z) <- ext; z}
bigsst <- merge(set_extent(sst, extent(sst) + c(-360, -360, 0, 0)), sst, set_extent(sst, extent(sst) + c(360, 360, 0, 0)))
m.sst <- mean(bigsst, na.rm = T)
msst <- raster::resample(m.sst, grd, method = "bilinear")
rm(m.sst)

msst <- mask(msst, msk) # mask
SST <- naProp(msst, msk) # check threshold
names(SST) <- "SST"

#######################
## SST gradient
sstg <- stack(bigsst)
for(i in 1:nlayers(sstg)){
  sstg[[i]] <- terrain(bigsst[[i]], opt = "slope", neighbors = 8)
}
msstg <- mean(sstg, na.rm = T)
msstg <- raster::resample(msstg, grd, method = "bilinear")
projection(msstg) <- projection(bigsst)

msstg <- mask(msstg, msk) # mask
SSTg <- naProp(msstg, msk) # check threshold
names(SSTg) <- "SSTg"

rm(sst, sstg, msst, msstg)

#######################
## SSHa
ssha <-readssh(stg_date, xylim=extent(grd360), ssha=T, lon180 = FALSE, latest = FALSE)
ssha <- rotate(ssha)
mssha <- mean(ssha, na.rm = T)
mssha <- resample(mssha, grd, method = "bilinear")

mssha <- mask(mssha, msk) # mask
SSHa <- naProp(mssha, msk) # check threshold
names(SSHa) <- "SSHa"

rm(mssha)

#######################
## SSHa sd
v.ssha <- calc(ssha, sd, na.rm = T)
msshv <- resample(v.ssha, grd, method = "bilinear")
rm(ssha, v.ssha)

msshv <- mask(msshv, msk) # mask
SSHsd <- naProp(msshv, msk) # check threshold
names(SSHsd) <- "SSHsd"

rm(msshv)

#######################
## Wind
wind <- mean(readwind(stg_date, magonly = TRUE, lon180 = TRUE, latest = FALSE), na.rm = T) # don't use xylim, crops edge
mwind <- resample(wind, grd, method = "bilinear")
rm(wind)

mwind <- mask(mwind, msk) # mask
WIND <- naProp(mwind, msk) # check threshold
names(WIND) <- "WIND"

rm(mwind)


#######################
## Eddy Kinetic Energy
ekev <- readcurr(stg_date, xylim=extent(grd), vonly = T, lon180 = TRUE, latest = FALSE)
ekeu <- readcurr(stg_date, xylim=extent(grd), uonly = T, lon180 = TRUE, latest = FALSE)
sekev <- stack(ekev)
sekeu <- stack(ekeu)
for(i in 1:nlayers(sekev)){
  sekev[[i]] <- 1/2*(sekev[[i]]^2 + sekeu[[i]]^2)
}
meke <- mean(sekeu, na.rm = T)
meke <- resample(meke, grd, method = "bilinear")

meke <- mask(meke, msk) # mask
EKE <- naProp(meke, msk) # check threshold
names(EKE) <- "EKE"

rm(ekev, ekeu, sekev, sekeu, meke)

#######################
## Current
m.curr <- mean(readcurr(stg_date, xylim=extent(grd), magonly = TRUE, lon180 = TRUE, latest = F), na.rm = T)
mcurr <- resample(m.curr, grd, method = "bilinear")
rm(m.curr)

mcurr <- mask(mcurr, msk) # mask
CURR <- naProp(mcurr, msk) # check threshold
names(CURR) <- "CURR"

rm(mcurr)

#######################
## Ice
#ice <- mean(readice(stg_date, latest = FALSE), na.rm = T)
ice <- mean(readice(stg_date, setNA = T, latest = FALSE), na.rm = T)
rd3dummy <- grd
projection(grd) <- "+proj=longlat +datum=WGS84"
rd3dummy <- grd
rd3dummy[] <- 0
rd3points <-  rasterToPoints(rd3dummy, spatial = TRUE)
mice <- extract(ice, rd3points, method = "bilinear")
mice <- setValues(grd, mice)
mice[is.na(mice)] <- 0
rm(ice)

mice <- mask(mice, msk) # mask
ICE <- naProp(mice, msk) # check threshold
names(ICE) <- "ICE"

rm(rd3dummy, rd3points, mice)

#######################
## Ice sd
ice <- calc(readice(stg_date, latest = FALSE), sd, na.rm = T)
rd3dummy <- grd
projection(grd) <- "+proj=longlat +datum=WGS84"
rd3dummy <- grd
rd3dummy[] <- 0
rd3points <-  rasterToPoints(rd3dummy, spatial = TRUE)
vice <- extract(ice, rd3points, method = "bilinear")
vice <- setValues(grd, vice)
vice[is.na(vice)] <- 0
rm(ice)

vice <- mask(vice, msk) # mask
ICEsd <- naProp(vice, msk) # check threshold
names(ICEsd) <- "ICEsd"

rm(rd3dummy, rd3points, vice)

#######################
## Ice CV
# ice <- calc(readice(stg_date, latest = F), cv, na.rm = T)
# rd3dummy <- grd
# projection(grd) <- "+proj=longlat +datum=WGS84"
# rd3dummy <- grd
# rd3dummy[] <- 0
# rd3points <-  rasterToPoints(rd3dummy, spatial = TRUE)
# cice <- extract(ice, rd3points, method = "bilinear")
# cice <- setValues(grd, cice)
# cice[is.na(cice)] <- 0
# rm(ice)
# 
# cice <- mask(cice, msk) # mask
# ICEcv <- naProp(cice, msk) # check threshold
# names(ICEcv) <- "ICEcv"
# 
# rm(rd3dummy, rd3points, cice)

#######################
## Chlorophyll-a

chla <- mean(readchla(stg_date, xylim=extent(grd), algorithm="johnson", latest = F), na.rm = T) #uses the extent of the original, un-projected grid, defined above
mchla <- resample(chla, grd, method = "bilinear")
# cgrid <- grd
# cgrid[] <- 0
# chla <- readoc(stg_date, platform = "SeaWiFS", grid = cgrid, latest = FALSE)
# mchla <- chla
rm(chla)

mchla <- mask(mchla, msk) # mask
CHLA <- naProp(mchla, msk) # check threshold
names(CHLA) <- "CHLA"

rm(mchla)

#######################
## Ice accessibility

icefolder <- "/perm_storage/home/shared/data_extra/accessibility/"
icefiles <- list.files(icefolder, full.names = T)
icetoget <- unique(format(stg_date, "%Y_%m"))
filestoget <- unique (grep(paste(icetoget, collapse = "|"), icefiles, value=TRUE))
rm(icefolder, icefiles, icetoget)

ice.access <- do.call(stack, lapply(filestoget, readRDS))
rm(filestoget)
ice.access <- mean(ice.access)

## this does not cope with dateline
##ice.access <- projectRaster(from = ice.access, to = grd, method = "bilinear")
dummy <- grd
dummy[] <- 0
dummy_points <-  rasterToPoints(dummy, spatial = TRUE)
ice.access <- extract(ice.access, dummy_points, method = "bilinear")
ice.access <- setValues(grd, ice.access)
## check
#prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"
#temp <- projectRaster(mice, crs=prj)
#plot(temp, col = viridis(101))
ice.access[is.na(ice.access)] <- maxValue(ice.access) ## so that non-ice-zone areas get filled with values

ICEA <- mask(ice.access, msk) # mask
ICEA <- naProp(ICEA, msk) # check threshold
names(ICEA) <- "ICEA"

rm(ice.access)

#######################
## Salinity difference between 200 m and 600 m

# From World Ocean Atlas 2013 v2
# https://www.nodc.noaa.gov/OC5/woa13/

sal.files <- list.files("/rdsi/PUBLIC/raad/data/data.nodc.noaa.gov/woa/WOA13/DATAv2/salinity/netcdf/decav/0.25/", full.names = T)
saltoget <- unique(format(stg_date[!is.na(stg_date)], "%m"))
saltoget <- paste0("/rdsi/PUBLIC/raad/data/data.nodc.noaa.gov/woa/WOA13/DATAv2/salinity/netcdf/decav/0.25//woa13_decav_s",
                   saltoget,
                   "_04v2.nc")

sal200 <- do.call(stack, lapply(saltoget, raster, varname="s_an", level=25))
sal600 <- do.call(stack, lapply(saltoget, raster, varname="s_an", level=39))

sal <- sal600 - sal200

sal <- mean(sal)

crs(sal) <- "+proj=longlat +datum=WGS84 +no_defs"

sal <- raster::resample(sal, grd, method = "bilinear")

sal <- mask(sal, msk) # mask
SAL <- naProp(sal, msk) # check threshold
names(SAL) <- "SAL"

rm(sal.files, saltoget, sal200, sal600, sal)


#######################
## Read in static layers

# Bathymetry
bath <- readderivaadc("bathymetry", xylim = lims)
bath <- resample(bath, grd, method = "bilinear")
bath <- crop(bath, extent(grd))

bath <- mask(bath, msk) # mask
DEPTH <- naProp(bath, msk) # check threshold
names(DEPTH) <- "DEPTH"

rm(bath)

# Slope
bathg <- readderivaadc("bathymetry_slope", xylim = lims)
bathg <- resample(bathg, grd, method = "bilinear")
bathg <- crop(bathg, extent(grd))

bathg <- mask(bathg, msk) # mask
DEPTHg <- naProp(bathg, msk) # check threshold
names(DEPTHg) <- "DEPTHg"

rm(bathg)

# Distance to shelf
dstslope <- readderivaadc("distance_shelf", xylim = lims)
dstslope <- resample(dstslope, grd, method = "bilinear")
dstslope <- crop(dstslope, extent(grd))

dstslope <- mask(dstslope, msk) # mask
dSHELF <- naProp(dstslope, msk) # check threshold
names(dSHELF) <- "dSHELF"

rm(dstslope)

# Vertical velocity at 500 m
# vertvel <- readderivaadc("vertical_velocity_500", xylim = lims)
# vertvel <- resample(vertvel, grd, method = "bilinear")
# vertvel <- crop(vertvel, extent(grd))
# 
# vertvel <- mask(vertvel, msk) # mask
# VERTVEL <- naProp(vertvel, msk) # check threshold
# names(VERTVEL) <- "VERTVEL"
# 
# rm(vertvel)

# Distance to ice-edge
# icedist <- readderivaadc("distance_max_ice_edge", xylim = lims)
# icedist <- resample(icedist, grd, method = "bilinear")
# icedist <- crop(icedist, extent(grd))
# 
# icedist <- mask(icedist, msk) # mask
# dICE <- naProp(icedist, msk) # check threshold
# names(dICE) <- "dICE"
# 
# rm(icedist)

# SH Flux
library(R.utils)
link <- paste0("http://webdav.data.aad.gov.au/data/environmental/derived_v2/common_grid/",
               "ecaisom_shflux_mean.asc.gz")
download.file(link, "temp.asc.gz")
gunzip("temp.asc.gz", remove = F, overwrite = T)
myras <- readAll(raster("temp.asc"))
#plot(clamp(myras, -200, 50))
file.remove("temp.asc.gz")
file.remove("temp.asc")

flux <- raster::resample(myras, grd, method = "bilinear")
flux <- crop(flux, extent(grd))

flux <- mask(flux, msk) # mask
SHFLUX <- naProp(flux, msk) # check threshold
names(SHFLUX) <- "SHFLUX"

rm(flux, myras)

# SH Flux SD
#library(R.utils)
link <- paste0("http://webdav.data.aad.gov.au/data/environmental/derived_v2/common_grid/",
               "ecaisom_shflux_sd.asc.gz")
download.file(link, "temp.asc.gz")
gunzip("temp.asc.gz", remove = F, overwrite = T)
myras <- readAll(raster("temp.asc"))
#plot(clamp(myras, -200, 50))
file.remove("temp.asc.gz")
file.remove("temp.asc")

flux <- raster::resample(myras, grd, method = "bilinear")
flux <- crop(flux, extent(grd))

flux <- mask(flux, msk) # mask
SHFLUXsd <- naProp(flux, msk) # check threshold
names(SHFLUXsd) <- "SHFLUXsd"

rm(flux, myras)

# VMIX
#library(R.utils)
link <- paste0("http://webdav.data.aad.gov.au/data/environmental/derived_v2/common_grid/",
               "ecaisom_vmix_250_mean.asc.gz")
download.file(link, "temp.asc.gz")
gunzip("temp.asc.gz", remove = F, overwrite = T)
myras <- readAll(raster("temp.asc"))
#plot(clamp(myras, -200, 50))
file.remove("temp.asc.gz")
file.remove("temp.asc")

mix <- raster::resample(myras, grd, method = "bilinear")
mix <- crop(mix, extent(grd))

mix <- mask(mix, msk) # mask
VMIX <- naProp(mix, msk) # check threshold
names(VMIX) <- "VMIX"

rm(mix, myras)

# VMIX SD
#library(R.utils)
link <- paste0("http://webdav.data.aad.gov.au/data/environmental/derived_v2/common_grid/",
               "ecaisom_vmix_250_sd.asc.gz")
download.file(link, "temp.asc.gz")
gunzip("temp.asc.gz", remove = F, overwrite = T)
myras <- readAll(raster("temp.asc"))
#plot(clamp(myras, -200, 50))
file.remove("temp.asc.gz")
file.remove("temp.asc")

mix <- raster::resample(myras, grd, method = "bilinear")
mix <- crop(mix, extent(grd))

mix <- mask(mix, msk) # mask
VMIXsd <- naProp(mix, msk) # check threshold
names(VMIXsd) <- "VMIXsd"

rm(mix, myras)
#-----------------------------------------------

# Function to fix plotting
fix_land_blip <- function(rst) {
  temp <- matrix(values(rst), nrow = dim(rst)[1], byrow = TRUE)
  temp[380:400, 1:13] <- NA_real_
  values(rst) <- temp
  rst
}

## Combine into a single stack

# If some rasters were set to NULL
vars <- c("SST",
          "SSTg",
          "SSHa",
          "SSHsd",
          "WIND",
          "EKE",
          "CURR",
          "ICE",
          "ICEsd",
          "ICEA",
          "CHLA",
          "DEPTH",
          "DEPTHg",
          "dSHELF",
          "SHFLUX",
          "SHFLUXsd",
          "VMIX",
          "VMIXsd",
          "SAL")

vars <- vars[order(vars)]

envdata <- stack()

for (i in 1:length(vars)) {
  if(length(get0(vars[i])) > 0) {
    envdata <- stack(envdata, get0(vars[i]))
  }
}

save(envdata, file = "./envdataAll.Rdata")
#load("./envdataAll.Rdata")

# rm(SST,
#   SSTg,
#   SSHa,
#   SSHsd,
#   WIND,
#   EKE,
#   CURR,
#   ICE,
#   ICEsd,
#   ICEA,
#   CHLA,
#   DEPTH,
#   DEPTHg,
#   dSHELF,
#   SHFLUX,
#   SHFLUXsd,
#   VMIX,
#   VMIXsd,
#   SAL)
#-----------------------------------------------

## Map
## Re-project to polar projection
projection(envdata) <- "+proj=longlat +datum=WGS84"
prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"

## Set up the world map 
# Create polar circle
library(rgeos)
library(geosphere)
data(countriesLow, package = "rworldmap")
polar.crop <- SpatialPoints(cbind(0, 0), proj4string = CRS(prj)) %>%
  gBuffer(., quadsegs = 1000, width = distGeo(cbind(0, -90), cbind(0, -41.5)))

## use a zero width buffer to ensure there's no topology problems (long story)
wrld <- spTransform(countriesLow, prj) %>%
  gBuffer(., width = 0, quadsegs = 1000) %>%
  crop(., polar.crop)

##set up fronts
of <- crop(orsifronts, extent(envdata))
ofp <- spTransform(of, CRS(projection(prj)))
ofp <- spTransform(ofp, CRS(prj))

## fix the land blip

pdiv <- projectRaster(envdata, crs=prj)


pdf("envdataAll.pdf", width=8.27, height=9.5)
par(mfrow=c(5,4), mar=c(0,1,0.8,6), cex.main = 1, font.main = 1)
for(i in 1:dim(pdiv)[3]){
  plot(pdiv[[i]], col=viridis(125), axes = FALSE, main=names(pdiv)[i], bty = "n", box = FALSE,
       legend = FALSE)
  plot(pdiv[[i]], col=viridis(125), axes = FALSE, main=names(pdiv)[i], bty = "n", box = FALSE,
       legend.only = T, legend.width = 2)
  #box(col = "white", lwd = 5)
  plot(wrld, add=T, col="darkgrey", border = FALSE)
  plot(ofp, add=T, col="#E6E6E6")
}
dev.off()
#-----------------------------------------------