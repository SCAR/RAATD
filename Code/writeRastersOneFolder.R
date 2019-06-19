## Mean layers

## Set working directory
setwd("~/RAATD_01/RAATD/allSpRasters")

library(raster)

# --------------------------------------------
# Function to load raster layers for all species and write them to a single folder
stkR <- function(sp, colony.weighted = FALSE) {
  fls <- list.files(paste0("../", sp, "/modPreds/"))
  if (colony.weighted) {
    these.files <- fls[grepl("_gbmcw_rast_trans.grd", fls)]
  } else {
    these.files <- fls[grepl("_gbm_rast_trans.grd", fls)]
  }
  for (j in seq_along(these.files)) {
    r <- raster(paste0("../", sp, "/modPreds/", these.files[j]))
    writeRaster(r, these.files[j], format = "raster", overwrite = T)
  }
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
                "ROPE",
                "SOES",
                "WAAL",
                "WESE",
                "WHCP")

lapply(sp.list, stkR, colony.weighted = FALSE)
lapply(sp.list, stkR, colony.weighted = TRUE)
