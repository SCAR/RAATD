## Fast run-through to produce quick discussion results,
## based on species groups

library(raster)
library(ggplot2)
library(viridis)
library(sf)

setwd("~/RAATD_01/RAATD/additionalAnalysis/fastRun/")

# ---------------------------
# Specific stuff
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
## Mapping stuff

# Get Antarctic mask
landmask <- raster("../../antarcticMask/mask.grd")

# Projections
projW <- "+proj=longlat +datum=WGS84 +no_defs"
projP <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"

# ---------------------------
## Generate maps with AES

# a. Weighted
datW <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")
datW <- meanR(datW, species.groups)
datWR <- rasterFromXYZ(datW[ , c("x", "y", "MEAN")])

datWR <- extend(datWR, landmask, NA) # raster has one less rows
datWR <- mask(datWR, landmask) # mask

polyW75 <- aesPoly(datWR, 0.75)
polyW90 <- aesPoly(datWR, 0.90)

saveRDS(polyW75, "polyW75.RDS")
saveRDS(polyW90, "polyW90.RDS")

st_write(polyW75, "polyW75.shp")
st_write(polyW90, "polyW90.shp")

trans.none <- rgb(0, 0, 0, alpha = 255/255, maxColorValue = 255)

# Project
crs(datWR) <- projW
datWR <- projectRaster(datWR, res = 2500, crs = projP, over = F)

polyW75 <- st_transform(polyW75, projP)
polyW90 <- st_transform(polyW90, projP)

png("MapWeighted.png", width = 1000, height = 1000, units = "px", res = 150)
plot(datWR, col = viridis(125), main = "Weighted mean habitat importance",
     axes = FALSE,
     bty = "n",
     box = FALSE)
plot(polyW75, add = T, col = trans.none, border = "grey75")
plot(polyW90, add = T, col = trans.none, border = "grey100")
dev.off()

# b. Non-weighted
datNW <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformed.RDS")
datNW <- meanR(datNW, species.groups)
datNWR <- rasterFromXYZ(datNW[ , c("x", "y", "MEAN")])

datNWR <- extend(datNWR, landmask, NA) # raster has one less rows
datNWR <- mask(datNWR, landmask) # mask

polyNW75 <- aesPoly(datNWR, 0.75)
polyNW90 <- aesPoly(datNWR, 0.90)

saveRDS(polyNW75, "polyNW75.RDS")
saveRDS(polyNW90, "polyNW90.RDS")

st_write(polyW75, "polyNW75.shp")
st_write(polyW90, "polyNW90.shp")

trans.none <- rgb(0, 0, 0, alpha = 255/255, maxColorValue = 255)

# Project
crs(datNWR) <- projW
datNWR <- projectRaster(datNWR, res = 2500, crs = projP, over = F)

polyNW75 <- st_transform(polyNW75, projP)
polyNW90 <- st_transform(polyNW90, projP)

png("MapUnweighted.png", width = 1000, height = 1000, units = "px", res = 150)
plot(datNWR, col = viridis(125), main = "Unweighted mean habitat importance",
     axes = FALSE,
     bty = "n",
     box = FALSE)
plot(polyNW75, add = T, col = trans.none, border = "grey75")
plot(polyNW90, add = T, col = trans.none, border = "grey100")
dev.off()