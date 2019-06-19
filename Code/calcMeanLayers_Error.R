## Mean layers

## Set working directory
setwd("~/RAATD_01/RAATD/")

library(raster)
library(viridis)
library(orsifronts)

# --------------------------------------------
# Function to load rasters for a given species and add error in quadrature
stkR <- function(sp) {
  
  # Get a file to crop the output # Not sure why the error rasters are bigger
  r.fls <- list.files(paste0("./", sp, "/modPreds/"))
  r.fls <- r.fls[grepl("_gbmcw_rast_trans.grd", r.fls)][1]
  tmplt <- raster(paste0("./", sp, "/modPreds/", r.fls))
  
  if (sp != "MAPE") {
    
  fls <- list.files("./predictionsError/")
    these.files <- fls[grepl("_gbm_error_rast.grd", fls)]
    these.files <- these.files[grepl(sp, these.files)]
  stk <- stack()
  for (j in seq_along(these.files)) {
    r <- raster(paste0("./predictionsError/", these.files[j]))
    stk <- stack(stk, r)
  }
  
  if(nlayers(stk) == 1) {
    stk <- stk
  } else {
    tmp <- calc(stk, fun=function(x){x * x})
    tmp <- calc(tmp, fun=function(x){sqrt(sum(x))})
    tmp <- tmp/nlayers(stk)
    stk <- tmp
  }
  
  names(stk) <- sp
  }
  
  #------------------------------
  # For MAPE/ROPE
  
  # The means for MAPE and ROPE are combined
  if (sp == "MAPE") {
    
    #------------------------
    ## First, calculate which layer is used where
    # ROPE
    fls2 <- list.files(paste0("./", "ROPE", "/modPreds/"))
      these.files2 <- fls2[grepl("_gbmcw_rast_trans.grd", fls2)]
    stk2 <- stack()
    for (k in seq_along(these.files2)) {
      r <- raster(paste0("./", "ROPE", "/modPreds/", these.files2[k]))
      stk2 <- stack(stk2, r)
    }
    stk2 <- mean(stk2)
    
    # MAPE
    fls3 <- list.files(paste0("./", "MAPE", "/modPreds/"))
      these.files3 <- fls3[grepl("_gbmcw_rast_trans.grd", fls3)]
    stk3 <- stack()
    for (k in seq_along(these.files3)) {
      r <- raster(paste0("./", "MAPE", "/modPreds/", these.files3[k]))
      stk3 <- stack(stk3, r)
    }
    stk3 <- mean(stk3)
    
    which.rope <- stk2 - stk3
    which.rope[which.rope < 0] <- NA
    
    which.mape <- stk3 - stk2
    which.mape[which.mape < 0] <- NA
    
    #------------------------
    ## Next, get the uncertainties as above
    
    # MAPE
    fls <- list.files("./predictionsError/")
    these.files <- fls[grepl("_gbm_error_rast.grd", fls)]
    these.files <- these.files[grepl("MAPE", these.files)]
    stk <- stack()
    for (j in seq_along(these.files)) {
      r <- raster(paste0("./predictionsError/", these.files[j]))
      stk <- stack(stk, r)
    }
    
    if(nlayers(stk) == 1) {
      stk <- stk
    } else {
      tmp <- calc(stk, fun=function(x){x * x})
      tmp <- calc(tmp, fun=function(x){sqrt(sum(x))})
      tmp <- tmp/nlayers(stk)
      stk <- tmp
    }
    
    mape.error <- stk
    
    # ROPE
    fls <- list.files("./predictionsError/")
    these.files <- fls[grepl("_gbm_error_rast.grd", fls)]
    these.files <- these.files[grepl("ROPE", these.files)]
    stk <- stack()
    for (j in seq_along(these.files)) {
      r <- raster(paste0("./predictionsError/", these.files[j]))
      stk <- stack(stk, r)
    }
    
    if(nlayers(stk) == 1) {
      stk <- stk
    } else {
      tmp <- calc(stk, fun=function(x){x * x})
      tmp <- calc(tmp, fun=function(x){sqrt(sum(x))})
      tmp <- tmp/nlayers(stk)
      stk <- tmp
    }
    
    rope.error <- stk
    
    #------------------------
    ## Then mask each and add them together
    
    which.mape <- crop(which.mape, mape.error)
    mape.error <- mask(mape.error, which.mape)
    
    which.rope <- crop(which.rope, rope.error)
    rope.error <- mask(rope.error, which.rope)
    
    stk <- merge(mape.error, rope.error)
    
    names(stk) <- "MAPE&ROPE"
  }
  
  #------------------------------
  
  stk <- crop(stk, tmplt)
  
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
# -----------------
# Run colony-weighted

allSp <- do.call(stack, lapply(sp.list, stkR))

plot(allSp, col = viridis(125))


# Mask
landmask <- raster("antarcticMask/mask.grd")

# Plot in polar stereographic
projection(allSp) <- "+proj=longlat +datum=WGS84"
prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"
allSpProj <- projectRaster(allSp, crs = prj)

pdf("./meanPredictions/allSpeciesColonyWeighted_Errors.pdf", paper = "a4r")
#par(mfrow=c(3,6), mar=c(0.5,0.5,0.5,3))
for(i in 1:dim(allSpProj)[3]){
  plot(allSpProj[[i]], col=viridis(125), axes = FALSE, main=names(allSpProj)[i], bty = "n", box = FALSE)
  #box(col = "white", lwd = 5)
  #plot(mp, add=T, col="grey")
  #plot(ofp, add=T, col="white")
}
dev.off()

# Write output grid
outSp <- as.data.frame(rasterToPoints(allSp))

# Load environmental data
load("./envarsAll/envdataAll.Rdata")

# Extract mean environmental conditions
env <- raster::extract(envdata, outSp[ , c("x", "y")])
outSp <- cbind(outSp, env)
saveRDS(outSp, "./meanPredictions/griddedScoresEnvarsTransformedColonyWeighted_Errors.RDS")