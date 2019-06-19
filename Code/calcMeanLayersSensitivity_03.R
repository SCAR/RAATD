## Mean layers

## Set working directory
setwd("~/RAATD_01/RAATD/additionalAnalysis/stageSensitivity")

library(raster)
library(viridis)
library(orsifronts)

# --------------------------------------------
# Function to load rasters for a given species and calculate the mean
stkR <- function(sp, which.fun = "max") {
  fls <- list.files(paste0("./", sp, "/modPreds/"))
    these.files <- fls[grepl("_gbmcw_rast_trans.grd", fls)]
  stk <- stack()
  for (j in seq_along(these.files)) {
    r <- raster(paste0("./", sp, "/modPreds/", these.files[j]))
    stk <- stack(stk, r)
  }
  if (which.fun == "mean") {
    stk <- mean(stk)
  }
  if (which.fun == "max") {
    stk <- max(stk)
  }
  names(stk) <- sp
  
  # The means for MAPE and ROPE are combined
  if (sp == "MAPE") {
    fls2 <- list.files(paste0("./", "ROPE", "/modPreds/"))
      these.files2 <- fls2[grepl("_gbmcw_rast_trans.grd", fls2)]
    stk2 <- stack()
    for (k in seq_along(these.files2)) {
      r <- raster(paste0("./", "ROPE", "/modPreds/", these.files2[k]))
      stk2 <- stack(stk2, r)
    }
    if (which.fun == "mean") {
      stk2 <- mean(stk2)
    }
    if (which.fun == "max") {
      stk2 <- max(stk2)
    }
    
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
# 2. Run for all species

allSpMean <- do.call(stack, lapply(sp.list, stkR, which.fun = "mean"))
allSpMax<- do.call(stack, lapply(sp.list, stkR, which.fun = "max"))

# Look at the differences per species
plot(allSpMax - allSpMean, col = viridis(125))

# Mask
landmask <- raster("antarcticMask/mask.grd")
allSpMax <- mask(allSpMax, landmask)
allSpMean <- mask(allSpMean, landmask)

# Write output grid
spGridMean <- as.data.frame(rasterToPoints(allSpMean))
spGridMax <- as.data.frame(rasterToPoints(allSpMax))

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

# Projections
projW <- "+proj=longlat +datum=WGS84 +no_defs"
projP <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"

# ---------------------------
## Generate maps with AES

datMean <- meanR(spGridMean, species.groups)
datMax <- meanR(spGridMax, species.groups)

datMeanRaster <- rasterFromXYZ(datMean[ , c("x", "y", "MEAN")])
datMaxRaster <- rasterFromXYZ(datMax[ , c("x", "y", "MEAN")])

polyW90Mean <- aesPoly(datMeanRaster, 0.90)
polyW90Max <- aesPoly(datMaxRaster, 0.90)

trans.none <- rgb(0, 0, 0, alpha = 255/255, maxColorValue = 255)


# Project
crs(datMeanRaster) <- projW
crs(datMaxRaster) <- projW

datMeanRasterP <- projectRaster(datMeanRaster, res = 2500, crs = projP, over = F)
datMaxRasterP <- projectRaster(datMaxRaster, res = 2500, crs = projP, over = F)

polyW90MeanP <- st_transform(polyW90Mean, projP)
polyW90MaxP <- st_transform(polyW90Max, projP)

png("MapWeighted.png", width = 1000, height = 1000, units = "px", res = 150)
plot(datMaxRasterP - datMeanRasterP, col = viridis(125), main = "Habitat Importance",
     axes = FALSE,
     bty = "n",
     box = FALSE)
plot(polyW90MaxP, add = T, col = trans.none, border = "grey75")
plot(polyW90MeanP, add = T, col = trans.none, border = "grey100")
dev.off()


# Try normalizing each raster for plotting
rescale <- function(x) {
  (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))
}

datMeanRasterPnorm <- datMeanRasterP
values(datMeanRasterPnorm) <- rescale(values(datMeanRasterPnorm))

datMaxRasterPnorm <- datMaxRasterP
values(datMaxRasterPnorm) <- rescale(values(datMaxRasterPnorm))

png("MapDifference.png", width = 1500, height = 1500, units = "px", res = 150)
plot(datMaxRasterPnorm - datMeanRasterPnorm, col = rev(brewer.blues(125)), main = "Difference of normalized habitat importance\nmax(stages) - mean(stages)",
     axes = FALSE,
     bty = "n",
     box = FALSE)
plot(polyW90MaxP, add = T, col = trans.none, border = "black", lwd = 3)
plot(polyW90MaxP, add = T, col = trans.none, border = "#cc3311", lwd = 1.5)

plot(polyW90MeanP, add = T, col = trans.none, border = "black", lwd = 3)
plot(polyW90MeanP, add = T, col = trans.none, border = "white", lwd = 1.5)
dev.off()
