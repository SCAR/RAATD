# Calculating uncertainty of predictions using bootstrap

# Ryan Reisinger & Ben Raymond

# April 2019

library(caret)
library(gbm)

library(raster)

# Libraries for parallel processing
library(iterators)
library(foreach)
library(parallel)
library(doParallel)

#--------------------------------------------------------
# Define species and stage

this.species <- "SOES"
this.stage <- "post-moult"

runInPar <- FALSE # Run in parallel?
ncores <- 5 # Number of cores, will use this number only if it is less than available cores

runInBatch <- TRUE # Run in batches to avoid out-of-memory problem? # Note parallel is not possible if this is true

this.species.low <- tolower(this.species)

#--------------------------------------------------------

# Number of bootstraps

K <- 50

#-------------------
# Species and stage specific stuff:

setwd(paste0("~/RAATD_01/RAATD/", this.species))
sp <- this.species

#-------------------
# Source functions
# source("/perm_storage/home/shared/github/raatd_modelling/caret/prepCaretBoot.R")
source("~/RAATD_01/RAATD/Code/prepCaretBoot.R")

# Data for modelling
# Load usage
load(paste0(this.species, "-usage-", this.stage, ".Rdata")) # Object should be called "usage", re-name if not (stage-specific example)

# Coerce the tbl back to data.frame
usage <- as.data.frame(usage)

#-------------------
# Data for prediction
# Load the raster values for prediction
load(paste0("./envdata/", this.species.low, "_envdata_", this.stage, ".Rdata"))
stk <- stack(envdata)
rm(envdata)
pred.dat <- rasterToPoints(stk)
rm(stk)

#-------------------
# Get the hyperparameters from the previously fitted GBM/BRT model
if(this.species == "ANFS" & this.stage == "post-moult") {
  gbmFit <- readRDS(paste0("fittedMods/", this.species, "_", this.stage, "_gbm_1.RDS"))
} else if (this.species == "SOES" & this.stage == "post-moult") {
  gbmFit <- readRDS(paste0("fittedMods/", this.species, "_", this.stage, "_gbm_1.RDS"))
} else {
  gbmFit <- readRDS(paste0("fittedMods/", this.species, "_", this.stage, "_gbm.RDS"))
}

gbmFitTune <- gbmFit$bestTune
rm(gbmFit)
gc()

# Set up train control
tc <- trainControl(method = "none",
                   number = 1,
                   search = "grid",
                   classProbs = TRUE,
                   allowParallel = TRUE,
                   summaryFunction = twoClassSummary,
                   verboseIter = TRUE,
                   sampling = "down",
                   indexOut = NULL)


#--------------------------------------------------------
## BRT Model
#--------------------------------------------------------

## Fit and predict K times

# Set up parallel processing if specificied above,
# otherwise run sequentially

if (!runInBatch) {
  
  if (runInPar) {
    ncr <- detectCores() - 1
    ncr <- min(ncores, ncr)
    cluster <- makeCluster(ncr, type = "FORK", outfile = "")
    registerDoParallel(cluster)
  } else {
    registerDoSEQ()
  }
  
  # Fit and predict
  out <- foreach (i=1:K,
                  .combine='cbind',
                  .inorder=FALSE,
                  .packages = c("caret", "gbm")) %dopar% {
                    
                    # Prep the data and sample with replacement for bootstrap
                    dat <- prepCaretBoot(usage = usage, sze = 0.5)
                    
                    # Fit the model
                    resp.col <- ncol(dat)
                    gbmMod <- caret::train(x = as.data.frame(dat[ ,c(8:(resp.col-2))]),
                                           y = dat[ , resp.col],
                                           method = "gbm",
                                           metric = "ROC",
                                           trControl = tc,
                                           tuneGrid = gbmFitTune)
                    
                    rm(dat)
                    gc()
                    
                    # Predict the model
                    pred.gbm <- caret::predict.train(gbmMod, newdata = pred.dat, type = "prob")
                    pred.gbm <- pred.gbm$Y
                    
                    pred.gbm
                  }
  
  # Save
  saveRDS(out, file = paste0("../predictionsError/", this.species, "_", this.stage, "_gbm_out.RDS"))
  
  # Stop the cluster
  if (runInPar) {
    stopCluster(cluster)
  }
  
}


# Run sub-batches if specified
if (runInBatch) {
  
  for (j in 1:50) {
    
    print(j)
    
    # Fit and predict
    # Prep the data and sample with replacement for bootstrap
    dat <- prepCaretBoot(usage = usage, sze = 0.5)
    
    # Fit the model
    resp.col <- ncol(dat)
    gbmMod <- caret::train(x = as.data.frame(dat[ ,c(8:(resp.col-2))]),
                           y = dat[ , resp.col],
                           method = "gbm",
                           metric = "ROC",
                           trControl = tc,
                           tuneGrid = gbmFitTune)
    
    rm(dat)
    gc()
    
    # Predict the model
    pred.gbm <- caret::predict.train(gbmMod, newdata = pred.dat, type = "prob")
    pred.gbm <- pred.gbm$Y
    
    out <- pred.gbm
    
    # Save
    saveRDS(out, file = paste0("../predictionsError/", this.species, "_", this.stage, "_gbm_out_", j, ".RDS"))
    
    rm(gbmMod, pred.gbm)
    gc()
    
  }
  
  # Gather the output
  foosa <- matrix(NA, nrow = nrow(pred.dat), ncol = 50)
  
  for (m in 1:50) {
    print(m)
    gthr <- readRDS(paste0("../predictionsError/", this.species, "_", this.stage, "_gbm_out_", m, ".RDS"))
    foosa[ , m] <- gthr
    rm(gthr)
  }
  
  # Save as for the normal, parallel-run output
  saveRDS(foosa, file = paste0("../predictionsError/", this.species, "_", this.stage, "_gbm_out.RDS"))
  
}
