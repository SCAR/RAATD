## Compare species groups as is, and
## when leaving out HUWH and SOES

library(raster)
library(ggplot2)
library(viridis)
library(sf)
library(pals)

library(rgeos)
library(geosphere)

setwd("~/RAATD_01/RAATD/additionalAnalysis/compareSpeciesGroupings/")

# ---------------------------
# Specific stuff

## Groups as they are now:
group.1 <- c("ADPE", "ANPE", "CRAS", "EMPE", "HUWH", "SOES", "WESE")
group.2 <- c("ANFS", "BBAL", "DMSA", "GHAL", "HUWH", "KIPE", "LMSA", "MAPE.ROPE", "SOES", "WAAL", "WHCP")

species.groups.all <- list(group.1, group.2)

## Groups without HUWH and SOES:
group.1x <- c("ADPE", "ANPE", "CRAS", "EMPE", "WESE")
group.2x <- c("ANFS", "BBAL", "DMSA", "GHAL", "KIPE", "LMSA", "MAPE.ROPE", "WAAL", "WHCP")

species.groups.x <- list(group.1x, group.2x)

## Groups with HUWH and SOES in their own group:
group.1y <- c("ADPE", "ANPE", "CRAS", "EMPE", "WESE")
group.2y <- c("ANFS", "BBAL", "DMSA", "GHAL", "KIPE", "LMSA", "MAPE.ROPE", "WAAL", "WHCP")
group.3y <- c("HUWH", "SOES")

species.groups.y <- list(group.1y, group.2y, group.3y)
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

# Create polar circle
data(countriesLow, package = "rworldmap")
polar.crop <- SpatialPoints(cbind(0, 0), proj4string = CRS(projP))
polar.crop <- gBuffer(polar.crop, quadsegs = 1000, width = distGeo(cbind(0, -90), cbind(0, -41.5)))

## use a zero width buffer to ensure there's no topology problems (long story)
wrld <- spTransform(countriesLow, projP)
wrld <- gBuffer(wrld, width = 0, quadsegs = 1000)
wrld <- crop(wrld, polar.crop)


# Plots
# ---------------------------
## Generate maps with AES

# --------------------------------------
# Current grouping
datW <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")
datW <- meanR(datW, species.groups.all)
datWR <- rasterFromXYZ(datW[ , c("x", "y", "MEAN")])

datWR <- extend(datWR, landmask, NA) # raster has one less rows
datWR <- mask(datWR, landmask) # mask

polyW90 <- aesPoly(datWR, 0.90)

# saveRDS(polyW90, "polyW90.RDS")
# st_write(polyW90, "polyW90.shp")

trans.none <- rgb(0, 0, 0, alpha = 255/255, maxColorValue = 255)

# Project
crs(datWR) <- projW
datWR <- projectRaster(datWR, res = 2500, crs = projP, over = F)

polyW75 <- st_transform(polyW75, projP)
polyW90 <- st_transform(polyW90, projP)

png("MapWeighted.png", width = 1000, height = 1000, units = "px", res = 150)
plot(datWR, col = viridis(125), main = "Weighted mean habitat importance \n Original grouping",
     axes = FALSE,
     bty = "n",
     box = FALSE)
plot(polyW90, add = T, col = trans.none, border = "grey100")
dev.off()


# --------------------------------------
# Without SOES and HUWH
datWx <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")
datWx <- meanR(datWx, species.groups.x)
datWRx <- rasterFromXYZ(datWx[ , c("x", "y", "MEAN")])

datWRx <- extend(datWRx, landmask, NA) # raster has one less rows
datWRx <- mask(datWRx, landmask) # mask

polyW90x <- aesPoly(datWRx, 0.90)

# saveRDS(polyW90x, "polyW90x.RDS")
# st_write(polyW90x, "polyW90x.shp")

trans.none <- rgb(0, 0, 0, alpha = 255/255, maxColorValue = 255)

# Project
crs(datWRx) <- projW
datWRx <- projectRaster(datWRx, res = 2500, crs = projP, over = F)

polyW90x <- st_transform(polyW90x, projP)

png("MapWeightedReduced.png", width = 1000, height = 1000, units = "px", res = 150)
plot(datWRx, col = viridis(125), main = "Weighted mean habitat importance \n Without HUWH & SOES",
     axes = FALSE,
     bty = "n",
     box = FALSE)
plot(polyW90x, add = T, col = trans.none, border = "grey100")
dev.off()


# --------------------------------------
# Difference

trans.red <- rgb(204, 51, 17, alpha = 255*0.3, maxColorValue = 255)
trans.blue <- rgb(0, 119, 187, alpha = 255*0.3, maxColorValue = 255)

tiff("MapWeightedDifference.tiff",
    width = 7,
    height = 7,
    units = "in",
    res = 600,
    bg = "white")

l <- length(seq(-25, 25, by = 1))-1

plot(datWR - datWRx, col = brewer.prgn(l), main = "Difference between species grouping results\n(Current grouping - revised grouping)",
     axes = FALSE,
     bty = "n",
     breaks = seq(-25, 25, by = 1),
     legend = FALSE,
     box = FALSE)

plot(polyW90, add = T, col = trans.red, border = "#CC3311")
plot(polyW90x, add = T, col = trans.blue, border = "#0077BB")

plot(datWR - datWRx, col = brewer.prgn(l), main = "Difference between species grouping results",
     axes = FALSE,
     bty = "n",
     breaks = seq(-25, 25, by = 1),
     box = FALSE,
     legend.only=TRUE,
     horizontal = TRUE,
     legend.shrink=0.7,
     axis.args=list(at = c(-25, 0, +25),
                    labels = c(-25, 0, +25)),
     legend.args = list(text = 'Difference in habitat importance')
)

# add land
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)

# add border
plot(polar.crop,
     add = TRUE)

dev.off()

# --------------------------------------
# With SOES and HUWH in their own groups
datWy <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")
datWy <- meanR(datWy, species.groups.y)
datWRy <- rasterFromXYZ(datWy[ , c("x", "y", "MEAN")])

datWRy <- extend(datWRy, landmask, NA) # raster has one less rows
datWRy <- mask(datWRy, landmask) # mask

polyW90y <- aesPoly(datWRy, 0.90)

# saveRDS(polyW90y, "polyW90y.RDS")
# st_write(polyW90y, "polyW90y.shp")

trans.none <- rgb(0, 0, 0, alpha = 255/255, maxColorValue = 255)

# Project
crs(datWRy) <- projW
datWRy <- projectRaster(datWRy, res = 2500, crs = projP, over = F)

polyW90y <- st_transform(polyW90y, projP)

png("MapWeightedReduced3groups.png", width = 1000, height = 1000, units = "px", res = 150)
plot(datWRy, col = viridis(125), main = "Weighted mean habitat importance \n HUWH & SOES own group",
     axes = FALSE,
     bty = "n",
     box = FALSE)
plot(polyW90y, add = T, col = trans.none, border = "grey100")
dev.off()


# --------------------------------------
# Difference with three groups

trans.red <- rgb(204, 51, 17, alpha = 255*0.3, maxColorValue = 255)
trans.blue <- rgb(0, 119, 187, alpha = 255*0.3, maxColorValue = 255)

tiff("MapWeightedDifferenceThreeGroups.tiff",
     width = 7,
     height = 7,
     units = "in",
     res = 600,
     bg = "white")

l <- 125

plot(datWR - datWRy, col = rev(brewer.purples(l)), main = "Difference between species grouping results\n(Current grouping - revised grouping)",
     axes = FALSE,
     bty = "n",
     legend = FALSE,
     box = FALSE)

plot(polyW90, add = T, col = trans.red, border = "#CC3311")
plot(polyW90y, add = T, col = trans.blue, border = "#0077BB")

plot(datWR - datWRy, col = rev(brewer.purples(l)), main = "Difference between species grouping results",
     axes = FALSE,
     bty = "n",
     box = FALSE,
     legend.only=TRUE,
     horizontal = TRUE,
     legend.shrink=0.7,
     legend.args = list(text = 'Difference in habitat importance')
)

# add land
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)

# add border
plot(polar.crop,
     add = TRUE)

dev.off()