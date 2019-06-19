## Mean layers

## Set working directory
setwd("~/RAATD_01/RAATD/")

library(raster)
library(viridis)
library(orsifronts)

# --------------------------------------------
# Function to load rasters for a given species and calculate the mean
stkR <- function(sp, colony.weight = TRUE) {
  fls <- list.files(paste0("./", sp, "/modPreds/"))
  if (colony.weight) {
    these.files <- fls[grepl("_gbmcw_rast_trans.grd", fls)]
  } else {
    these.files <- fls[grepl("_gbm_rast_trans.grd", fls)]
  }
  stk <- stack()
  for (j in seq_along(these.files)) {
    r <- raster(paste0("./", sp, "/modPreds/", these.files[j]))
    stk <- stack(stk, r)
  }
  stk <- mean(stk)
  names(stk) <- sp
  
  # The means for MAPE and ROPE are combined
  if (sp == "MAPE") {
    fls2 <- list.files(paste0("./", "ROPE", "/modPreds/"))
    if (colony.weight) {
      these.files2 <- fls2[grepl("_gbmcw_rast_trans.grd", fls2)]
    } else {
      these.files2 <- fls2[grepl("_gbm_rast_trans.grd", fls2)]
    }
    stk2 <- stack()
    for (k in seq_along(these.files2)) {
      r <- raster(paste0("./", "ROPE", "/modPreds/", these.files2[k]))
      stk2 <- stack(stk2, r)
    }
    stk2 <- mean(stk2)
    
    stk <- overlay(stk, stk2, fun = "max") # Calculate the max value for each pixel
    names(stk) <- "MAPE&ROPE"
  }
  
  return(stk)
}

# --------------------------------------------
# Apply to all species and stack
sp.list <- list("ADPE",
                "ANFS",
                "ANPE",
                "BBAL",
                "CRAS",
                "DMSA",
                "EMPE",
                "GHAL",
                "HUWH",
                "KIPE",
                "LMSA",
                "MAPE",
                #"ROPE", # ROPE is combined with MAPE internally in the stkR function
                "SOES",
                "WAAL",
                "WESE",
                "WHCP")

# -----------------
# 1. Run non-colony weighted

allSp <- do.call(stack, lapply(sp.list, stkR, colony.weight = FALSE))

# allSpNotImp <- do.call(stack, lapply(sp.list, stkR, importance = FALSE))

plot(allSp, col = viridis(125))

# Calculate mean of all
all <- mean(allSp)

# Mask
landmask <- raster("antarcticMask/mask.grd")
all <- mask(all, landmask)

# Plot in polar stereographic

projection(allSp) <- "+proj=longlat +datum=WGS84"
prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"
allSpProj <- projectRaster(allSp, crs = prj)

pdf("./meanPredictions/allSpecies.pdf", paper = "a4r")
#par(mfrow=c(3,6), mar=c(0.5,0.5,0.5,3))
for(i in 1:dim(allSpProj)[3]){
  plot(allSpProj[[i]], col=viridis(125), axes = FALSE, main=names(allSpProj)[i], bty = "n", box = FALSE)
  #box(col = "white", lwd = 5)
  #plot(mp, add=T, col="grey")
  #plot(ofp, add=T, col="white")
}
dev.off()

# Write the raster
writeRaster(all, "./meanPredictions/meanAll.grd", format = "raster", overwrite = T)

# Write output grid
outSp <- as.data.frame(rasterToPoints(allSp))
outAll <- as.data.frame(rasterToPoints(all))

# Load environmental data
load("./envarsAll/envdataAll.Rdata")

# Extract mean environmental conditions
env <- raster::extract(envdata, outSp[ , c("x", "y")])
outSp <- cbind(outSp, env)
saveRDS(outSp, "./meanPredictions/griddedScoresEnvarsTransformed.RDS")


# --------------------------------------------
# Plot the mean layer

## Map
## Re-project to polar projection
projection(all) <- "+proj=longlat +datum=WGS84"
prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"

## Set up the world map 
#data(wrld_simpl)
#mp <- crop(wrld_simpl, extent(envdata))
#w <- spTransform(mp, CRS(projection(prj)))
#pw <- spTransform(w, CRS(prj))

##set up fronts
of <- crop(orsifronts, extent(all))
ofp <- spTransform(of, CRS(projection(prj)))
ofp <- spTransform(ofp, CRS(prj))

pdiv <- projectRaster(all, crs=prj)

pdf("./meanPredictions/meanMap.pdf", width=10, height=10)
plot(pdiv, col=viridis(125), axes = FALSE, main="Mean habitat importance", bty = "n", box = FALSE)
#plot(ofp, add=T, col="white")
dev.off()

jpeg("./meanPredictions/meanMap.jpeg", width=10, height=10, units = "in", res = 150)
plot(pdiv, col=viridis(125), axes = FALSE, main="Mean habitat importance", bty = "n", box = FALSE)
#plot(ofp, add=T, col="white")
dev.off()

# -----------------
# 2. Run colony-weighted

allSp <- do.call(stack, lapply(sp.list, stkR, colony.weight = TRUE))

# allSpNotImp <- do.call(stack, lapply(sp.list, stkR, importance = FALSE))

plot(allSp, col = viridis(125))

# Calculate mean of all
all <- mean(allSp)

# Mask
landmask <- raster("antarcticMask/mask.grd")
all <- mask(all, landmask)

# Plot in polar stereographic

projection(allSp) <- "+proj=longlat +datum=WGS84"
prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"
allSpProj <- projectRaster(allSp, crs = prj)

pdf("./meanPredictions/allSpeciesColonyWeighted.pdf", paper = "a4r")
#par(mfrow=c(3,6), mar=c(0.5,0.5,0.5,3))
for(i in 1:dim(allSpProj)[3]){
  plot(allSpProj[[i]], col=viridis(125), axes = FALSE, main=names(allSpProj)[i], bty = "n", box = FALSE)
  #box(col = "white", lwd = 5)
  #plot(mp, add=T, col="grey")
  #plot(ofp, add=T, col="white")
}
dev.off()

# Write the raster
writeRaster(all, "./meanPredictions/meanAllColonyWeighted.grd", format = "raster", overwrite = T)

# Write output grid
outSp <- as.data.frame(rasterToPoints(allSp))
outAll <- as.data.frame(rasterToPoints(all))

# Load environmental data
load("./envarsAll/envdataAll.Rdata")

# Extract mean environmental conditions
env <- raster::extract(envdata, outSp[ , c("x", "y")])
outSp <- cbind(outSp, env)
saveRDS(outSp, "./meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")


# --------------------------------------------
# Plot the mean layer

## Map
## Re-project to polar projection
projection(all) <- "+proj=longlat +datum=WGS84"
prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"

## Set up the world map 
#data(wrld_simpl)
#mp <- crop(wrld_simpl, extent(envdata))
#w <- spTransform(mp, CRS(projection(prj)))
#pw <- spTransform(w, CRS(prj))

##set up fronts
of <- crop(orsifronts, extent(all))
ofp <- spTransform(of, CRS(projection(prj)))
ofp <- spTransform(ofp, CRS(prj))

pdiv <- projectRaster(all, crs=prj)

pdf("./meanPredictions/meanMapColonyWeighted.pdf", width=10, height=10)
plot(pdiv, col=viridis(125), axes = FALSE, main="Mean weighted habitat importance", bty = "n", box = FALSE)
#plot(ofp, add=T, col="white")
dev.off()

jpeg("./meanPredictions/meanMapColonyWeighted.jpeg", width=10, height=10, units = "in", res = 150)
plot(pdiv, col=viridis(125), axes = FALSE, main="Mean weighted habitat importance", bty = "n", box = FALSE)
#plot(ofp, add=T, col="white")
dev.off()

# --------------------------------------
# Comparsion plot
library(grid)
library(jpeg)

p1 <-  grid::rasterGrob(jpeg::readJPEG("./meanPredictions/meanMap.jpeg"), interpolate = TRUE)
p2 <-  grid::rasterGrob(jpeg::readJPEG("./meanPredictions/meanMapColonyWeighted.jpeg"), interpolate = TRUE) 

tiff("./meanPredictions/compareMeanHabitat.tiff", width=10, height=5, units = "in", res = 150)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()