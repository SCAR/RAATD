## Create polygons of AESs, based on quantiles

setwd("~/RAATD_01/RAATD/additionalAnalysis/aesPoly/")

library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(ggplot2)
library(viridis)

library(sf)
library(liblwgeom)
library(smoothr)

# ------------------------------------
# 1. Non colony-weighted

# Get the mean habitat importance
dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedGroup.RDS")
dat$MEAN <- rowMeans(dat[ , 3:18])

# Create a raster
scores <- rasterFromXYZ(dat[ , c("x", "y", "MEAN")])
regions <- rasterFromXYZ(dat[ , c("x", "y", "group")])

myFunc <- function(which.region) {
  
  require(sf)
  require(liblwgeom)
  require(smoothr)
  
# Create copies
foo <- scores
bar <- regions

# Select the region of interest
bar[bar != which.region] <- NA

# Mask the scores
foo <- mask(foo, bar)

# Calculate quantiles
q <- quantile(foo, p=c(0.75))

# Threshold raster based on quantiles
foo[foo < q] <- NA
foo[!is.na(foo)] <- 1

# Convert to polygons
c <- rasterToPolygons(foo, dissolve = T)

# Set CRS
proj4string(c) <- "+proj=longlat +datum=WGS84 +no_defs"

# Convert to SF object 
r_poly <- st_as_sf(c)

# Drop small, isolated cells
area_thresh <- units::set_units(100*100, km^2)
r_poly_dropped <- drop_crumbs(r_poly, area_thresh)

# Fill small holes
r_poly_filled <- fill_holes(r_poly_dropped, area_thresh)

# And smooth the edges
r_poly_smooth <- smoothr::smooth(r_poly_filled, method = "ksmooth", smoothness = 20)

# Plot to check
plot(c, col = NA, border = NA) # set up plot extent
plot(r_poly_smooth, col = "#4DAF4A", border = "grey20", lwd = 1.5, add = TRUE)

return(r_poly_smooth)

}

region.1 <- myFunc(1)
region.2 <- myFunc(2)
region.3 <- myFunc(3)

saveRDS(region.1, "region1_poly_smooth.RDS")
saveRDS(region.2, "region2_poly_smooth.RDS")
saveRDS(region.3, "region3_poly_smooth.RDS")

colrs <- c("#0077BB", "#009988", "#33BBEE")
this.alpha = 255/3

col1 <- rgb(col2rgb(colrs[1])[1], col2rgb(colrs[1])[2], col2rgb(colrs[1])[3], alpha = this.alpha, maxColorValue = 255)
col2 <- rgb(col2rgb(colrs[2])[1], col2rgb(colrs[2])[2], col2rgb(colrs[2])[3], alpha = this.alpha, maxColorValue = 255)
col3 <- rgb(col2rgb(colrs[3])[1], col2rgb(colrs[3])[2], col2rgb(colrs[3])[3], alpha = this.alpha, maxColorValue = 255)

plot(scores, col = grey.colors(125))
plot(region.1, add = T, col = col1, border = F)
plot(region.2, add = T, col = col2, border = F)
plot(region.3, add = T, col = col3, border = F)

# ------------------------------------
# 2. Colony weighted

rm(list = ls())

# Get the mean habitat importance
dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeightedGroup.RDS")
dat$MEAN <- rowMeans(dat[ , 3:18])

# Create a raster
scores <- rasterFromXYZ(dat[ , c("x", "y", "MEAN")])
regions <- rasterFromXYZ(dat[ , c("x", "y", "group")])

myFunc <- function(which.region) {
  
  require(sf)
  require(liblwgeom)
  require(smoothr)
  
  # Create copies
  foo <- scores
  bar <- regions
  
  # Select the region of interest
  bar[bar != which.region] <- NA
  
  # Mask the scores
  foo <- mask(foo, bar)
  
  # Calculate quantiles
  q <- quantile(foo, p=c(0.75))
  
  # Threshold raster based on quantiles
  foo[foo < q] <- NA
  foo[!is.na(foo)] <- 1
  
  # Convert to polygons
  c <- rasterToPolygons(foo, dissolve = T)
  
  # Set CRS
  proj4string(c) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  # Convert to SF object 
  r_poly <- st_as_sf(c)
  
  # Drop small, isolated cells
  area_thresh <- units::set_units(100*100, km^2)
  r_poly_dropped <- drop_crumbs(r_poly, area_thresh)
  
  # Fill small holes
  r_poly_filled <- fill_holes(r_poly_dropped, area_thresh)
  
  # And smooth the edges
  r_poly_smooth <- smoothr::smooth(r_poly_filled, method = "ksmooth", smoothness = 20)
  
  # Plot to check
  plot(c, col = NA, border = NA) # set up plot extent
  plot(r_poly_smooth, col = "#4DAF4A", border = "grey20", lwd = 1.5, add = TRUE)
  
  return(r_poly_smooth)
  
}

region.1 <- myFunc(1)
region.2 <- myFunc(2)
region.3 <- myFunc(3)

saveRDS(region.1, "region1_poly_smooth_colonyweighted.RDS")
saveRDS(region.2, "region2_poly_smooth_colonyweighted.RDS")
saveRDS(region.3, "region3_poly_smooth_colonyweighted.RDS")

colrs <- c("#0077BB", "#009988", "#33BBEE")
this.alpha = 255/3

col1 <- rgb(col2rgb(colrs[1])[1], col2rgb(colrs[1])[2], col2rgb(colrs[1])[3], alpha = this.alpha, maxColorValue = 255)
col2 <- rgb(col2rgb(colrs[2])[1], col2rgb(colrs[2])[2], col2rgb(colrs[2])[3], alpha = this.alpha, maxColorValue = 255)
col3 <- rgb(col2rgb(colrs[3])[1], col2rgb(colrs[3])[2], col2rgb(colrs[3])[3], alpha = this.alpha, maxColorValue = 255)

plot(scores, col = grey.colors(125))
plot(region.1, add = T, col = col1, border = F)
plot(region.2, add = T, col = col2, border = F)
plot(region.3, add = T, col = col3, border = F)

# ------------------------------------
# 3. Plot overlap between colony weighted and non colony weighted AES

r1 <- readRDS("region1_poly_smooth.RDS")
r2 <- readRDS("region2_poly_smooth.RDS")
r3 <- readRDS("region3_poly_smooth.RDS")

r1cw <- readRDS("region1_poly_smooth_colonyweighted.RDS")
r2cw <- readRDS("region2_poly_smooth_colonyweighted.RDS")
r3cw <- readRDS("region3_poly_smooth_colonyweighted.RDS")


trans.blue <- rgb(0, 119, 187, alpha = this.alpha, maxColorValue = 255)
trans.red <- rgb(204, 51, 17, alpha = this.alpha, maxColorValue = 255)

# Region 1
tiff("compare_openocean.tiff", width = 12, height = 6, units = "in", res = 150)
plot(scores, col = grey.colors(125))
plot(r1, add = T, col = trans.red, border = F)
plot(r1cw, add = T, col = trans.blue, border = F)
dev.off()

# Region 2
tiff("compare_antarctic.tiff", width = 12, height = 6, units = "in", res = 150)
plot(scores, col = grey.colors(125))
plot(r2, add = T, col = trans.red, border = F)
plot(r2cw, add = T, col = trans.blue, border = F)
dev.off()

# Region 3
tiff("compare_subantarctic.tiff", width = 12, height = 6, units = "in", res = 150)
plot(scores, col = grey.colors(125))
plot(r3, add = T, col = trans.red, border = F)
plot(r3cw, add = T, col = trans.blue, border = F)
dev.off()