## Maps of threat index

setwd("~/RAATD_01/RAATD/additionalAnalysis/fastRun/")

require(tidyverse)
require(raster)
require(sp)
require(rgeos)
require(geosphere)
require(viridis)
require(orsifronts)
require(rworldmap)
library(raadtools)
library(colorspace)
library(sf)

## Map stuff

prj <- "+proj=laea +lat_0=-90 +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0"

# Create polar circle
data(countriesLow, package = "rworldmap")
polar.crop <- SpatialPoints(cbind(0, 0), proj4string = CRS(prj)) %>%
  gBuffer(., quadsegs = 1000, width = distGeo(cbind(0, -90), cbind(0, -41.5)))

## use a zero width buffer to ensure there's no topology problems (long story)
wrld <- spTransform(countriesLow, prj) %>%
  gBuffer(., width = 0, quadsegs = 1000) %>%
  crop(., polar.crop)

# Get bathymetry
grd <- raster(extent(-180,180,-90,-40), resolution = c(0.1, 0.1))
# bathy <- readbathy()
# bathy <- resample(bathy, grd, method = "bilinear")
# bathy[bathy > 0 ] <- 0
# crs(bathy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# bathy <- projectRaster(bathy, crs = prj)
# writeRaster(bathy, "bath.grd", format = "raster")

# Read in pre-processed bathymetry (preceding chunk)
bathy <- raster("~/RAATD_01/RAATD/additionalAnalysis/threats/threatIndexMap/bath.grd")

# Get Orsi fronts & reproject to Lambert equal area
orsi.polar <- spTransform(orsifronts, CRS(prj))

# Get the threat layers
fish <- raster("~/RAATD_01/RAATD/additionalAnalysis/threats/threatIndexMap/fish.grd")
ice <- raster("~/RAATD_01/RAATD/additionalAnalysis/threats/threatIndexMap/ice.grd")
sst <- raster("~/RAATD_01/RAATD/additionalAnalysis/threats/threatIndexMap/sst.grd")
wind <- raster("~/RAATD_01/RAATD/additionalAnalysis/threats/threatIndexMap/wind.grd")

# Mask ice
msk <- raster("/perm_storage/home/shared/data_extra/env_change/ice_mask.grd")
msk[msk <= 5] <- NA
msk <- projectRaster(msk, sst)
msk[msk <= 5] <- NA
ice <- mask(ice, msk)

# Cube root of fishing
# fish <- fish^(1/3)

# Mean threat
meanThreat <- sum(sqrt(ice * ice), sqrt(sst * sst), fish, na.rm = T)

# Project
fishP <- projectRaster(fish, crs = prj)
iceP <- projectRaster(ice, crs = prj)
sstP <- projectRaster(sst, crs = prj)
windP <- projectRaster(wind, crs = prj)

meanThreatP <- projectRaster(meanThreat, crs = prj)

# -------------------------------------------

## Function to plot
plotR <- function(thisRaster, thisFilename, thisTitle, colorStyle) {
  
  # Plot stuff
  
  par(cex.main = 1, font.main = 1)
  
  this.alpha = 255/255
  grey.trans <- rgb(125, 125, 125, alpha = this.alpha, maxColorValue = 255)

  if (colorStyle == "sst") {
    colrs <- diverge_hsv(125)
  } else if (colorStyle == "ice") {
    colrs <- rev(diverge_hsv(125))
  } else {
    colrs <- diverge_hsv(250)[125:250]
  }
  
  tiff(
    file = thisFilename,
    width = 7,
    height = 7,
    units = "in",
    res = 600,
    bg = "white"
  )
  
  plot(
    thisRaster,
    #main = thisTitle,
    col = colrs,
    axes = FALSE,
    bty = "n",
    box = FALSE,
    legend = FALSE,
    plt = c(1,1,1,1)
  )
  
  # Add the legend
  plot(
    thisRaster,
    #main = thisTitle,
    col = colrs,
    legend.only=TRUE,
    horizontal = TRUE,
    legend.shrink=0.7,
    legend.args = list(text = thisTitle)
  )
 
  plot(r75p,
       col = grey.trans, border = "grey30",
       add = TRUE)

  plot(r90p,
       col = grey.trans, border = "grey0",
       add = TRUE)
  
  # add land
  plot(wrld,
       col = "darkgrey",
       border = FALSE,
       add = TRUE)
  
  # border
  plot(polar.crop, add = T)
  
  dev.off()
  
}

# -------------------------------------------------
# 1. Weighted

# Get AES polygons
r75 <- readRDS("polyW75.RDS")
r90 <- readRDS("polyW90.RDS")

r75p <- st_transform(r75, prj)
r90p <- st_transform(r90, prj)

r75p <- as(r75p, "Spatial")
r90p <- as(r90p, "Spatial")

plotR(thisRaster = iceP, thisFilename = "iceWeighted.tiff", thisTitle = "Sea ice change index", colorStyle = "ice")
plotR(thisRaster = sstP, thisFilename = "sstWeighted.tiff", thisTitle = "SST change index", colorStyle = "sst")
plotR(thisRaster = windP, thisFilename = "windWeighted.tiff", thisTitle = "Wind change index", colorStyle = "sst")
plotR(thisRaster = fishP, thisFilename = "fishWeighted.tiff", thisTitle = "Fishing effort index", colorStyle = "fish")
plotR(thisRaster = meanThreatP, thisFilename = "mean.tiff", thisTitle = "Mean threat index", colorStyle = "fish")

# Quick composite
library(grid)
library(gridExtra)
library(tiff)

p1 <-  grid::rasterGrob(tiff::readTIFF("iceWeighted.tiff"))
p2 <-  grid::rasterGrob(tiff::readTIFF("sstWeighted.tiff")) 
p3 <-  grid::rasterGrob(tiff::readTIFF("windWeighted.tiff")) 
p4 <-  grid::rasterGrob(tiff::readTIFF("fishWeighted.tiff"))
p5 <-  grid::rasterGrob(tiff::readTIFF("meanWeighted.tiff"))

tiff("quickCompositeWeighted.tiff",
     width = 14,
     height = 14,
     units = "in",
     res = 150,
     bg = "white")
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

# -------------------------------------------------
# 2. Unweighted

# Get AES polygons
r75 <- readRDS("polyNW75.RDS")
r90 <- readRDS("polyNW90.RDS")

r75p <- st_transform(r75, prj)
r90p <- st_transform(r90, prj)

r75p <- as(r75p, "Spatial")
r90p <- as(r90p, "Spatial")

plotR(thisRaster = iceP, thisFilename = "iceUnweighted.tiff", thisTitle = "Sea ice change index", colorStyle = "ice")
plotR(thisRaster = sstP, thisFilename = "sstUnweighted.tiff", thisTitle = "SST change index", colorStyle = "sst")
plotR(thisRaster = windP, thisFilename = "windUnweighted.tiff", thisTitle = "Wind change index", colorStyle = "sst")
plotR(thisRaster = fishP, thisFilename = "fishUnweighted.tiff", thisTitle = "Fishing effort index", colorStyle = "fish")
plotR(thisRaster = meanThreatP, thisFilename = "mean.tiff", thisTitle = "Mean threat index", colorStyle = "fish")

# Quick composite
library(grid)
library(gridExtra)
library(tiff)

p1 <-  grid::rasterGrob(tiff::readTIFF("iceUnweighted.tiff"))
p2 <-  grid::rasterGrob(tiff::readTIFF("sstUnweighted.tiff")) 
p3 <-  grid::rasterGrob(tiff::readTIFF("windUnweighted.tiff")) 
p4 <-  grid::rasterGrob(tiff::readTIFF("fishUnweighted.tiff"))
p5 <-  grid::rasterGrob(tiff::readTIFF("meanUnweighted.tiff"))

tiff("quickCompositeUnweighted.tiff",
     width = 14,
     height = 14,
     units = "in",
     res = 150,
     bg = "white")
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()
