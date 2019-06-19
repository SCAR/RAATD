library(raster)
library(viridis)
library(colorspace)

# -------------------------------
this.species <- "ADPE"
this.stage <- "chick-rearing"

this.species.low <- tolower(this.species)

setwd(paste0("~/RAATD_01/RAATD/", this.species))

plot.sims <- FALSE

# ------------------------------------------------------------------------
# Read in the track data
datadir <- "/perm_storage/home/shared/github/raatd_data/"
ssm <- readRDS(paste0(datadir,
                      "data_filtered/filtered_by_stage/", this.species.low, "_ssm_by_stage.RDS"))

# Get the stages
stages <- unique(ssm$stage)

# Get the metadata
meta <- read.csv(paste0(datadir,
                        "metadata/SCAR_Metadata_2017_forWEBDAV.csv"),
                 stringsAsFactors = F)
meta <- meta[meta$abbreviated_name == this.species, ]

# Get colony locations
cloc <- "/perm_storage/home/shared/github/raatd_modelling/dataLayers/distanceColony/colonies.csv"
cloc <- read_csv(cloc) %>%
  filter(abbreviation == this.species) %>%
  filter(lat <= -40) %>%
  filter(count > 0) %>%
  dplyr::select("lon", "lat") %>%
  SpatialPoints(., proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Get simulated tracks, if neccessary
if(plot.sims) {
siml <- readRDS(paste0(this.species, "_simTracks.RDS"))
siml <- siml[siml$sim == 1, ]
siml <- siml[siml$stage == this.stage, ]
}

# --------------------------------
# Transform from list to dataframe
ssm$howlong <- NA

for (i in 1:nrow(ssm)) {
  ssm$howlong[i] <- length(ssm$ssm[[i]])
}

# Remove the failures
ssm <- ssm[which(ssm$howlong != 1), ]

# Get the IDs and stages
ids <- ssm$id
stage <- ssm$stage

# Bind them up
hold <- list()

for(i in 1:length(ids)){
  dat <-ssm$ssm[[i]]$predicted
  dat$id <- rep(ids[i], nrow(dat))
  dat$stage<- rep(stage[i], nrow(dat))
  hold[i] <- list(dat)
  rm(dat)
}

# Keep only data for this stage
data <- do.call(rbind, hold)
data <- data[data$stage == this.stage, ]

# Crop
data <- data[data$lat < -40, ]

if(plot.sims) {
  siml <- siml[siml$lat < -40, ]
}

# Keep only relevant metadata
meta <- meta[meta$individual_id %in% unique(data$id), ]

rm(hold)
rm(ssm)
rm(ids)

# ------------------------------------------------------------------------
# Get environmental data

load(paste0("envdata/", this.species.low, "_envdata_", this.stage, ".Rdata"))
stk <- stack(envdata)
bath <- stk[["DEPTH"]]
rm(envdata)

# ------------------------------------------------------------------------
# Get predictions

pred <- raster(paste0("modPreds/", this.species, "_", this.stage, "_gbm_rast.gri"))
av <- raster(paste0("modPreds/", this.species, "_", this.stage, "_av.grd"))
pred.gbm.av <- raster(paste0("modPreds/", this.species, "_", this.stage, "_gbm-av_rast.gri"))
pred.trans <- raster(paste0("modPreds/", this.species, "_", this.stage, "_gbm_rast_trans.gri"))

# Colony weighted layers
col <- raster(paste0("/perm_storage/home/shared/colonyWeighting/colonyWeights_", this.species, "_", this.stage, ".grd"))
pred.cw <- raster(paste0("modPreds/", this.species, "_", this.stage, "_gbmcw_rast_trans.gri"))

# And swimming distance
swim <- raster(paste0("/perm_storage/home/shared/distanceColony/distanceColony_", this.species, ".grd"))
projection(swim) <- "+proj=longlat +datum=WGS84"

# ------------------------------------------------------------------------
# Plot

pdf(file = paste0(this.species, "_", this.stage, "_map.pdf"), paper = "a4")

par(mfrow = c(3, 1))
plot(bath, main = paste0("Tracks - ", this.species, " - ", this.stage), col = gray.colors(125))
points(data$lon, data$lat, cex = 0.01)
points(cloc, col = "gold", pch = 16, cex = 0.7)
points(x = meta$deployment_decimal_longitude, y = meta$deployment_decimal_latitude, col = "red", pch = 16, cex = 0.7)

plot(pred, main = paste0("Habitat preference - ", this.species, " - ", this.stage), col = viridis(125))
points(cloc, col = "red", pch = 16, cex = 0.7)
points(x = meta$deployment_decimal_longitude, y = meta$deployment_decimal_latitude, col = "red", pch = 16, cex = 0.7)

plot(av, main = paste0("Availability - ", this.species, " - ", this.stage), col = sequential_hcl(125))
points(cloc, col = "red", pch = 16, cex = 0.7)
points(x = meta$deployment_decimal_longitude, y = meta$deployment_decimal_latitude, col = "red", pch = 16, cex = 0.7)

plot(col, main = paste0("Weighted availability - ", this.species, " - ", this.stage), col = sequential_hcl(125))
points(cloc, col = "red", pch = 16, cex = 0.7)
points(x = meta$deployment_decimal_longitude, y = meta$deployment_decimal_latitude, col = "red", pch = 16, cex = 0.7)

plot(pred.gbm.av, main = paste0("Habitat preference given availability - ", this.species, " - ", this.stage), col = viridis(125))
points(cloc, col = "red", pch = 16, cex = 0.7)
points(x = meta$deployment_decimal_longitude, y = meta$deployment_decimal_latitude, col = "red", pch = 16, cex = 0.7)

plot(pred.trans, main = paste0("Habitat importance - ", this.species, " - ", this.stage), col = viridis(125))
points(cloc, col = "red", pch = 16, cex = 0.7)
points(x = meta$deployment_decimal_longitude, y = meta$deployment_decimal_latitude, col = "red", pch = 16, cex = 0.7)

plot(pred.cw, main = paste0("Weighted habitat importance - ", this.species, " - ", this.stage), col = viridis(125))
points(cloc, col = "red", pch = 16, cex = 0.7)
points(x = meta$deployment_decimal_longitude, y = meta$deployment_decimal_latitude, col = "red", pch = 16, cex = 0.7)

dev.off()

# ------------------------------------------------------------------------
# Polar plot

projection(bath) <- "+proj=longlat +datum=WGS84"
prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"

## Set up the world map 
# Create polar circle
library(rgeos)
library(geosphere)
library(tidyverse)
data(countriesLow, package = "rworldmap")
polar.crop <- SpatialPoints(cbind(0, 0), proj4string = CRS(prj)) %>%
  gBuffer(., quadsegs = 1000, width = distGeo(cbind(0, -90), cbind(0, -41.5)))

## use a zero width buffer to ensure there's no topology problems (long story)
wrld <- spTransform(countriesLow, prj) %>%
  gBuffer(., width = 0, quadsegs = 1000) %>%
  crop(., polar.crop)

## Project rasters
bathP <- projectRaster(bath, crs=prj)
predP <- projectRaster(pred, crs=prj)
avP <- projectRaster(av, crs=prj)
pred.gbm.avP <- projectRaster(pred.gbm.av, crs=prj)
pred.transP <- projectRaster(pred.trans, crs=prj)
swimP <- projectRaster(swim, crs = prj)
colP <- projectRaster(col, crs = prj)
pred.cwP <- projectRaster(pred.cw, crs = prj)

##set up fronts
library(orsifronts)
of <- crop(orsifronts, extent(bath))
ofp <- spTransform(of, CRS(projection(prj)))
ofp <- spTransform(ofp, CRS(prj))

## Tracks and deployment locations
library(sp)
tracks <- SpatialPoints(coords = data[ , c("lon", "lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))
tracks <- spTransform(tracks, prj)

deploy <- SpatialPoints(coords = meta[ , c("deployment_decimal_longitude", "deployment_decimal_latitude")], proj4string = CRS("+proj=longlat +datum=WGS84"))
deploy <- spTransform(deploy, prj)

## Simulated tracks
if (plot.sims) {
sims <- SpatialPoints(coords = siml[ , c("lon", "lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))
sims <- spTransform(sims, prj)
}

## Colony locations
cloc <- "/perm_storage/home/shared/github/raatd_modelling/dataLayers/distanceColony/colonies.csv"
cloc <- read_csv(cloc) %>%
  filter(abbreviation == this.species) %>%
  filter(lat <= -40) %>%
  filter(count > 0) %>%
  dplyr::select("lon", "lat") %>%
  SpatialPoints(., proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")) 
cloc <- sp::spTransform(cloc, CRS(prj))

# -------------------------------
# -------------------------------
# Plot to pdf
pdf(file = paste0(this.species, "_", this.stage, "_map_polar.pdf"), paper = "a4")
#par(mfrow=c(3,2), mar=c(2.5,2,1.5,0.5), cex.main = 1, font.main = 1)
par(cex.main = 1, font.main = 1)

# -------------------------------
# 1. Tracks
plot(bathP,
     col = colorRampPalette(c("dodgerblue", "white"), alpha = TRUE)(125),
     axes = FALSE,
     bty = "n",
     box = FALSE,
     main = paste0("Tracks - ", this.species, " - ", this.stage))

# Orsi fronts
plot(ofp,
     add=T,
     col="#E6E6E6")

# Add sims
if(plot.sims) {
points(
  x = sims$lon,
  y = sims$lat,
  pch = 16,
  cex = 0.3,
  col = "grey"
)
}

# Add tracks
points(
  x = tracks$lon,
  y = tracks$lat,
  pch = 16,
  cex = 0.3,
  col = rgb(50, 50, 50, 51, maxColorValue = 255)
)

# World map
plot(wrld,
     add=T,
     col="darkgrey",
     border = FALSE)

# Colony locations
points(
  x = cloc$lon,
  y = cloc$lat,
  pch = 19,
  cex = 0.6,
  col = "gold"
)

# Deployment locations
points(
  x = deploy$deployment_decimal_longitude,
  y = deploy$deployment_decimal_latitude,
  pch = 19,
  cex = 0.6,
  col = "red"
)

# -------------------------------  
# 2a. Availability
plot(avP,
     col = viridis(125),
     axes = FALSE,
     bty = "n",
     box = FALSE,
     main = paste0("Habitat availability - ", this.species, " - ", this.stage))

# Orsi fronts
plot(ofp,
     add=T,
     col="#E6E6E6")

# World map
plot(wrld,
     add=T,
     col="darkgrey",
     border = FALSE)

# Colony locations
points(
  x = cloc$lon,
  y = cloc$lat,
  pch = 19,
  cex = 0.6,
  col = "red"
)

# -------------------------------  
# 2b. Weighted Availability
plot(colP,
     col = viridis(125),
     axes = FALSE,
     bty = "n",
     box = FALSE,
     main = paste0("Weighted habitat availability - ", this.species, " - ", this.stage))

# Orsi fronts
plot(ofp,
     add=T,
     col="#E6E6E6")

# World map
plot(wrld,
     add=T,
     col="darkgrey",
     border = FALSE)

# Colony locations
points(
  x = cloc$lon,
  y = cloc$lat,
  pch = 19,
  cex = 0.6,
  col = "red"
)

# -------------------------------  
# 3. Habitat selectivity
plot(predP,
     col = viridis(125),
     axes = FALSE,
     bty = "n",
     box = FALSE,
     main = paste0("Habitat selectivity - ", this.species, " - ", this.stage))

# Orsi fronts
plot(ofp,
     add=T,
     col="#E6E6E6")

# World map
plot(wrld,
     add=T,
     col="darkgrey",
     border = FALSE)

# Colony locations
points(
  x = cloc$lon,
  y = cloc$lat,
  pch = 19,
  cex = 0.6,
  col = "red"
)

# -------------------------------  
# 4. Habitat selectivity given availability
plot(pred.gbm.avP,
     col = viridis(125),
     axes = FALSE,
     bty = "n",
     box = FALSE,
     main = paste0("Habitat selectivity given availability - ", this.species, " - ", this.stage))

# Orsi fronts
plot(ofp,
     add=T,
     col="#E6E6E6")

# World map
plot(wrld,
     add=T,
     col="darkgrey",
     border = FALSE)

# Colony locations
points(
  x = cloc$lon,
  y = cloc$lat,
  pch = 19,
  cex = 0.6,
  col = "red"
)

# -------------------------------  
# 5. Habitat importance
plot(pred.transP,
     col = viridis(125),
     axes = FALSE,
     bty = "n",
     box = FALSE,
     main = paste0("Habitat importance - ", this.species, " - ", this.stage))

# Orsi fronts
plot(ofp,
     add=T,
     col="#E6E6E6")

# World map
plot(wrld,
     add=T,
     col="darkgrey",
     border = FALSE)

# Colony locations
points(
  x = cloc$lon,
  y = cloc$lat,
  pch = 19,
  cex = 0.6,
  col = "red"
)

# -------------------------------  
# 6. Habitat importance (weighted)
plot(pred.cwP,
     col = viridis(125),
     axes = FALSE,
     bty = "n",
     box = FALSE,
     main = paste0("Weighted habitat importance - ", this.species, " - ", this.stage))

# Orsi fronts
plot(ofp,
     add=T,
     col="#E6E6E6")

# World map
plot(wrld,
     add=T,
     col="darkgrey",
     border = FALSE)

# Colony locations
points(
  x = cloc$lon,
  y = cloc$lat,
  pch = 19,
  cex = 0.6,
  col = "red"
)

dev.off()
#-----------------------------------------------


# # For overview figure...
# 
# # Mask
# landmask <- raster("../antarcticMask/mask.grd")
# landmaskP <- projectRaster(from = landmask, swimP)
# swimP <- mask(swimP, landmaskP)
# 
# # Plot to pdf
# pdf(file = "hold.pdf", paper = "a4")
# par(mfrow=c(3,2), mar=c(2.5,2,1.5,0.5), cex.main = 1, font.main = 1)
# 
# # -------------------------------
# # 1. Tracks
# plot(swimP,
#      col = viridis(125),
#      axes = FALSE,
#      bty = "n",
#      box = FALSE,
#      main = paste0("Tracks - ", this.species, " - ", this.stage))
# 
# # Orsi fronts
# plot(ofp,
#      add=T,
#      col="#E6E6E6")
# 
# # World map
# plot(wrld,
#      add=T,
#      col="darkgrey",
#      border = FALSE)
# 
# # Colony locations
# points(
#   x = cloc$lon,
#   y = cloc$lat,
#   pch = 19,
#   cex = 0.6,
#   col = "red"
# )
# 
# dev.off()
# # -------------------------------  