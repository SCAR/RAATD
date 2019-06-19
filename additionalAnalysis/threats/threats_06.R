## Scatterplots and maps of threats and context against habitat importance
## Threat index

setwd("~/RAATD_01/RAATD/additionalAnalysis/threats/")

library(raster)
require(rgeos)
require(geosphere)
library(ggplot2)
# library(colorplaner)
library(colorspace)
library(viridis)
library(sf)
library(pals)

library(paletteer)



# ---------------------------
# Specific stuff

weighted <- TRUE

group.1 <- c("ADPE", "ANPE", "CRAS", "EMPE", "HUWH", "SOES", "WESE")
group.2 <- c("ANFS", "BBAL", "DMSA", "GHAL", "HUWH", "KIPE", "LMSA", "MAPE.ROPE", "SOES", "WAAL", "WHCP")

species.groups <- list(group.1, group.2)

# ---------------------------
## Source functions

# Function to calculate AES
source("~/RAATD_01/RAATD/Code/function_aesPoly.R")

# Function to combine different species groups
source("~/RAATD_01/RAATD/Code/function_meanR.R")

# ---------------------------
# Calculate habitat importance

if (weighted) {
  dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")
  dat <- meanR(dat, species.groups)
} else {
  dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformed.RDS")
  dat <- meanR(dat, species.groups)
}


# ---------------------------
# Mask
landmask <- raster("../../antarcticMask/mask.grd")

# Get threat data
fldr <- paste0("/perm_storage/home/shared/data_extra/env_change/")

# For plotting, calculate deciles of mean habitat importance
qnt <- quantile(dat$MEAN, na.rm = T)

#-------------------------------
#-------------------------------
## Threat data

# Fishing
load(paste0(fldr, "fishing_effort_.5deg.Rdata"))
r[is.na(r)] <- 0
fish <- r
# fish <- resample(fish, landmask)
# fish <- mask(fish, landmask)
rm(r)

# Write as raster
# writeRaster(fish, "fishing.grd", format = "raster")

dat$fish <- raster::extract(fish, dat[ , c("x", "y")])

#-------------------------------
# SST change
sst <- raster(paste0(fldr, "sst_change.grd"))
sst <- resample(sst, landmask)
sst <- mask(sst, landmask)
dat$sst <- raster::extract(sst, dat[ , c("x", "y")])

#-------------------------------
# Wind change
wind <- raster(paste0(fldr, "wind.grd"))
dat$wind <- raster::extract(wind, dat[ , c("x", "y")])

#-------------------------------
# Ice change

msk <- raster("/perm_storage/home/shared/data_extra/env_change/ice_mask.grd")
msk[msk <= 5] <- NA
msk <- projectRaster(msk, sst)
msk[msk <= 5] <- NA

ice <- raster(paste0(fldr, "ice_duration_change.grd"))
ice <- projectRaster(ice, to = sst)
ice <- mask(ice, msk)
ice <- mask(ice, landmask)
dat$ice <- raster::extract(ice, dat[ , c("x", "y")])

#-------------------------------
# Halpern cumulative impact
impact <- raster(paste0(fldr, "HalpernCumulativeThreats.grd"))

# Projected version for visualization
impactProj <- projectRaster(impact, crs = "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs")
impactProj[impactProj <= 0] <- NA
png("cumulative_impact_polar_map.png", width = 800, height = 800, units = "px")
plot(impactProj, col = viridis(125))
dev.off()

dat$impact <- raster::extract(impact, dat[ , c("x", "y")])

#-------------------------------
#-------------------------------
## Maps of threats, for paper

# Project the layers

prj <- "+proj=laea +lat_0=-90 +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0"

landmaskP <- projectRaster(landmask, crs = prj)

projection(fish) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
projection(ice) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
projection(sst) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
projection(wind) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

fishP <- projectRaster(fish, landmaskP, method = "ngb")
iceP <- projectRaster(ice, landmaskP)
sstP <- projectRaster(sst, landmaskP)
windP <- projectRaster(wind, landmaskP)

fishP <- mask(fishP, landmaskP)
iceP <- mask(iceP, landmaskP)
sstP <- mask(sstP, landmaskP)
windP <- mask(windP, landmaskP)

# ------------------------
# Plots

# Create polar circle
data(countriesLow, package = "rworldmap")
polar.crop <- SpatialPoints(cbind(0, 0), proj4string = CRS(prj))
polar.crop <- gBuffer(polar.crop, quadsegs = 1000, width = distGeo(cbind(0, -90), cbind(0, -41.5)))

## use a zero width buffer to ensure there's no topology problems (long story)
wrld <- spTransform(countriesLow, prj)
wrld <- gBuffer(wrld, width = 0, quadsegs = 1000)
wrld <- crop(wrld, polar.crop)

# Get Orsi fronts & reproject
library(orsifronts)
orsi.polar <- spTransform(orsifronts, CRS(prj))

# Get AES

# r1 <- readRDS("../aesPoly/region1_poly_smooth.RDS")
# r2 <- readRDS("../aesPoly/region2_poly_smooth.RDS")
# r3 <- readRDS("../aesPoly/region3_poly_smooth.RDS")
# 
# r1p <- st_transform(r1, prj)
# r2p <- st_transform(r2, prj)
# r3p <- st_transform(r3, prj)
# 
# r1p <- as(r1p, "Spatial")
# r2p <- as(r2p, "Spatial")
# r3p <- as(r3p, "Spatial")

if (weighted) {
  poly75 <- readRDS("~/RAATD_01/RAATD/additionalAnalysis/fastRun/polyW75.RDS")
  poly90 <- readRDS("~/RAATD_01/RAATD/additionalAnalysis/fastRun/polyW90.RDS")
  poly75 <- st_transform(poly75, prj)
  poly90 <- st_transform(poly90, prj)
  poly75 <- as(poly75, "Spatial")
  poly90 <- as(poly90, "Spatial")
} else {
  poly75 <- readRDS("~/RAATD_01/RAATD/additionalAnalysis/fastRun/polyNW75.RDS")
  poly90 <- readRDS("~/RAATD_01/RAATD/additionalAnalysis/fastRun/polyNW90.RDS")
  poly75 <- st_transform(poly75, prj)
  poly90 <- st_transform(poly90, prj)
  poly75 <- as(poly75, "Spatial")
  poly90 <- as(poly90, "Spatial")
}

trans.orange <- rgb(238, 119, 51, alpha = 255/3, maxColorValue = 255)
this.orange <- "#EE7733"
trans.none <- rgb(0, 0, 0, alpha = 255/255, maxColorValue = 255)

# ------------------------
# Fishing
# Inverse sine hyperbolic transform
# https://robjhyndman.com/hyndsight/transformations/
# http://wresch.github.io/2013/03/08/asinh-scales-in-ggplot2.html

fishP <- asinh(fishP)

tiff(
  file = "polar_fishing.tiff",
  width = 7,
  height = 7,
  units = "in",
  res = 600,
  bg = "white",
  type = "cairo"
)

# prediction
plot(
  fishP,
  col = viridis(125),
  axes = FALSE,
  bty = "n",
  box = FALSE,
  legend = FALSE,
  plt = c(1,1,1,1)
)

# Add the legend
plot(
  fishP,
     col = viridis(125),
     legend.only=TRUE,
     horizontal = TRUE,
     legend.shrink=0.7,
  axis.args=list(at = c(0, asinh(10), asinh(100), asinh(1000), asinh(10000), asinh(100000)),
                 labels = c(0, 10, 100, 1000, 10000, "100000")),
  legend.args = list(text = 'Fishing effort (hours)')
)

# Add AES
# plot(r1p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(r2p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(r3p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")

# plot(poly75, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(poly90, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")

# plot(poly90, add = T, col = trans.orange, border = this.orange)

plot(poly90, add = T, col = trans.none, border = "white")

# add land
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)

# add border
plot(polar.crop,
     add = TRUE)

# add Orsi fronts
# plot(subset(orsi.polar, front != "stf"),
#      col = grey(0.4),
#      lwd = 0.75,
#      add = TRUE)

# add colony locations
# points(
#   x = cloc$lon,
#   y = cloc$lat,
#   pch = 19,
#   cex = 0.5,
#   col = "gold"
# )


dev.off()


# ------------------------
# Ice

iceP[iceP > 50] <- 50
iceP[iceP < -50] <- -50

tiff(
  file = "polar_ice.tiff",
  width = 7,
  height = 7,
  units = "in",
  res = 600,
  bg = "white",
  type = "cairo"
)

l <- length(seq(-50, 50, by = 1))-1

# prediction
plot(
  iceP,
  col = paletteer_c(scico, cork, l),
  #col=diverge_hsv(l),
  #col = brewer.prgn(l),
  axes = FALSE,
  bty = "n",
  breaks = seq(-50, 50, by = 1),
  box = FALSE,
  legend = FALSE,
  plt = c(1,1,1,1)
)

# Add the legend
plot(
  iceP,
  col = paletteer_c(scico, cork, l),
  #col=diverge_hsv(l),
  #col = brewer.prgn(l),
  breaks = seq(-50, 50, by = 1),
  legend.only=TRUE,
  horizontal = TRUE,
  legend.shrink=0.7,
  axis.args=list(at = c(-50, -25, 0, +25, 50),
                 labels = c("\u2264 -50", -25, 0, +25, "50 \u2264")),
  legend.args = list(text = 'Change in mean ice duration (days)')
)

# Add AES
# plot(r1p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(r2p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(r3p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")

# plot(poly75, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(poly90, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")

# plot(poly90, add = T, col = trans.orange, border = this.orange)

plot(poly90, add = T, col = trans.none, border = "black")

# add land
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)

# add border
plot(polar.crop,
     add = TRUE)

# add Orsi fronts
# plot(subset(orsi.polar, front != "stf"),
#      col = grey(0.4),
#      lwd = 0.75,
#      add = TRUE)

# add colony locations
# points(
#   x = cloc$lon,
#   y = cloc$lat,
#   pch = 19,
#   cex = 0.5,
#   col = "gold"
# )


dev.off()


# ------------------------
# SST

sstP[sstP > 2] <- 2
sstP[sstP < -2] <- -2

tiff(
  file = "polar_sst.tiff",
  width = 7,
  height = 7,
  units = "in",
  res = 600,
  bg = "white",
  type = "cairo"
)

l <- length(seq(-2, 2, by = 0.05))-1

# prediction
plot(
  sstP,
  col = paletteer_c(scico, cork, l),
  #col=diverge_hsv(l),
  #col = brewer.prgn(l),
  axes = FALSE,
  bty = "n",
  breaks = seq(-2, 2, by = 0.05),
  box = FALSE,
  legend = FALSE,
  plt = c(1,1,1,1)
)

# Add the legend
plot(
  sstP,
  col = paletteer_c(scico, cork, l),
  #col=diverge_hsv(l),
  #col = brewer.prgn(l),
  breaks = seq(-2, 2, by = 0.05),
  legend.only=TRUE,
  horizontal = TRUE,
  legend.shrink=0.7,
  axis.args=list(at = c(-2, -1, 0, 1, 2),
                 labels = c("\u2264 -2", -1, 0, 1, "2 \u2264")),
  legend.args = list(text = 'Change in mean SST (\u00B0)')
)

# Add AES
# plot(r1p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(r2p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(r3p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")

# plot(poly75, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(poly90, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")

plot(poly90, add = T, col = trans.none, border = "black")

# add land
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)

# add border
plot(polar.crop,
     add = TRUE)

# add Orsi fronts
# plot(subset(orsi.polar, front != "stf"),
#      col = grey(0.4),
#      lwd = 0.75,
#      add = TRUE)

# add colony locations
# points(
#   x = cloc$lon,
#   y = cloc$lat,
#   pch = 19,
#   cex = 0.5,
#   col = "gold"
# )


dev.off()

# ------------------------
# Wind

tiff(
  file = "polar_wind.tiff",
  width = 7,
  height = 7,
  units = "in",
  res = 600,
  bg = "white",
  type = "cairo"
)

l <- length(seq(-2, 2, by = 0.05))-1

# prediction
plot(
  windP,
  col = paletteer_c(scico, cork, l),
  #col=diverge_hsv(l),
  #col = brewer.prgn(l),
  axes = FALSE,
  bty = "n",
  breaks = seq(-2, 2, by = 0.05),
  box = FALSE,
  legend = FALSE,
  plt = c(1,1,1,1)
)

# Add the legend
plot(
  windP,
  col = paletteer_c(scico, cork, l),
  #col=diverge_hsv(l),
  #col = brewer.prgn(l),
  breaks = seq(-2, 2, by = 0.05),
  legend.only=TRUE,
  horizontal = TRUE,
  legend.shrink=0.7,
  axis.args=list(at = pretty(-2:2),
                 labels = pretty(-2:2)),
  legend.args = list(text = 'Change in mean wind speed (m/s)')
)

# Add AES
# plot(r1p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(r2p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(r3p, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")

# plot(poly75, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")
# plot(poly90, add = T, col = adjustcolor( "grey", alpha.f = 0.4), border = "grey30")

plot(poly90, add = T, col = trans.none, border = "black")

# add land
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)

# add border
plot(polar.crop,
     add = TRUE)

# add Orsi fronts
# plot(subset(orsi.polar, front != "stf"),
#      col = grey(0.4),
#      lwd = 0.75,
#      add = TRUE)

# add colony locations
# points(
#   x = cloc$lon,
#   y = cloc$lat,
#   pch = 19,
#   cex = 0.5,
#   col = "gold"
# )


dev.off()