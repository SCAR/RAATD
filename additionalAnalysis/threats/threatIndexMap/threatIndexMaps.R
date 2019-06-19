## Maps of threat index

setwd("~/RAATD_01/RAATD/additionalAnalysis/threats/threatIndexMap")

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
bathy <- raster("bath.grd")

# Get Orsi fronts & reproject to Lambert equal area
orsi.polar <- spTransform(orsifronts, CRS(prj))

# Get AES polygons
library(sf)
r1 <- readRDS("../../aesPoly/region1_poly_smooth.RDS")
r2 <- readRDS("../../aesPoly/region2_poly_smooth.RDS")
r3 <- readRDS("../../aesPoly/region3_poly_smooth.RDS")

r1p <- st_transform(r1, prj)
r2p <- st_transform(r2, prj)
r3p <- st_transform(r3, prj)

r1p <- as(r1p, "Spatial")
r2p <- as(r2p, "Spatial")
r3p <- as(r3p, "Spatial")

# Get the threat layers
fish <- raster("fish.grd")
ice <- raster("ice.grd")
sst <- raster("sst.grd")
wind <- raster("wind.grd")

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
  
  this.alpha = 255*0.3
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
 
  plot(r2p,
       col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30",
       add = TRUE)

  plot(r1p,
       col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30",
       add = TRUE)

  plot(r3p,
       col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30",
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

plotR(thisRaster = iceP, thisFilename = "ice.tiff", thisTitle = "Sea ice change index", colorStyle = "ice")
plotR(thisRaster = sstP, thisFilename = "sst.tiff", thisTitle = "SST change index", colorStyle = "sst")
plotR(thisRaster = windP, thisFilename = "wind.tiff", thisTitle = "Wind change index", colorStyle = "sst")
plotR(thisRaster = fishP, thisFilename = "fish.tiff", thisTitle = "Fishing effort index", colorStyle = "fish")
plotR(thisRaster = meanThreatP, thisFilename = "mean.tiff", thisTitle = "Mean threat index", colorStyle = "fish")

# Quick composite
library(grid)
library(gridExtra)
library(tiff)

p1 <-  grid::rasterGrob(tiff::readTIFF("ice.tiff"))
p2 <-  grid::rasterGrob(tiff::readTIFF("sst.tiff")) 
p3 <-  grid::rasterGrob(tiff::readTIFF("wind.tiff")) 
p4 <-  grid::rasterGrob(tiff::readTIFF("fish.tiff"))
p5 <-  grid::rasterGrob(tiff::readTIFF("mean.tiff"))

tiff("quickComposite.tiff",
     width = 14,
     height = 14,
     units = "in",
     res = 300,
     bg = "white")
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()