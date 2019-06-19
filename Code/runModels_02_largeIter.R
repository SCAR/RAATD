# Modelling of usage data

# This code runs various habitat models on the ouput of the Marseille package (S. Wotherspoon)
# The Marseille package output includes environmental data extracted in a previous step
# Script is species and stage-specific

# Ryan Reisinger, Mark Hindell and Ian Jonsen

# last updated 2018-02-13

library(caret)
library(pdp)
library(gbm)

# Libraries for parallel processing
library(iterators)
library(foreach)
library(parallel)
library(doParallel)

for (i in 1:20) {
  
  print(i)

#--------------------------------------------------------
# Prepare

this.species <- "ANFS"
this.stage <- "post-moult"

tune <- FALSE # Tune the model, or use fixed parameters?
tooBig <- TRUE # If the data are too large, sample only a proportion thereof

#-------------------
# Species and stage specific stuff:

setwd(paste0("~/RAATD_01/RAATD/", this.species))

sp <- this.species

# Load usage
load(paste0(this.species, "-usage-", this.stage, ".Rdata")) # Object should be called "usage", re-name if not (stage-specific example)

#-------------------
# Coerce the tbl back to data.frame
usage <- as.data.frame(usage)

# Prepare data.frame for caret using the 'prepCaret' function
source("/perm_storage/home/shared/github/raatd_modelling/caret/prepCaret3.R")
dat <- prepCaret(usage = usage, toobig = tooBig, prp = 0.3) # toobig allows sumbsampling of massive data frames

#-------------------
# Choose only one presence and absence per cell
#dat <- dat[!duplicated(dat), ]

#-------------------
# Create CV folds using the 'foldCaret' function
# source("/perm_storage/home/shared/github/raatd_modelling/caret/foldCaretSpaceBlock.R")
# folds <- foldCaret(dat = dat, nm = 5, spaceBlock = FALSE, blockRes = NULL)

#-------------------
# Clean up
rm(usage)

#--------------------------------------------------------
#Models
#--------------------------------------------------------

#--------------------------------
# Load global caret parameters
source("/perm_storage/home/shared/github/raatd_modelling/caret/caretParameters.R")

#--------------------------------
# If not parameter tuninG

if (tune == FALSE) {
  gbmGrid <-  expand.grid(interaction.depth = 3, 
                          n.trees = 5000, 
                          shrinkage = 0.01,
                          n.minobsinnode = 20)
}

#--------------------------------
# Set up parallel processing

cluster <- makeCluster(detectCores() - 1) # leave 1 core for OS
registerDoParallel(cluster)

#-------------------
#1. BRT/GBM

resp.col <- ncol(dat)

system.time(
  gbmMod <- train(x = as.data.frame(dat[ ,c(8:(resp.col-2))]),
                        y = dat[ , resp.col],
                        method = "gbm",
                        metric = "ROC",
                        trControl = tc,
                        tuneGrid = gbmGrid)
)

gbmMod
summary(gbmMod)
varImp(gbmMod)
saveRDS(gbmMod, paste0("fittedMods/", this.species, "_", this.stage, "_gbm_", i, ".RDS"))
#gbmMod <- readRDS(paste0("fittedMods/", this.species, "_", this.stage, "_gbm.RDS"))

#--------------------------------------------------------
# Plot summaries

# source("../Code/poorManPartial.R")
# 
# # Warning: choosing to plot the pdps with pdp package (plot_pdp = TRUE) is VERY slow!
# # A plotmo plot can be produced (plot_plotmo = TRUE), which is far less accurate, but way faster
# poorManPartial(mod = gbmMod,
#                dat = dat,
#                plot_var = TRUE,
#                varimpFilename = paste0("varImpPlot_", this.species, "_", this.stage, ".pdf"),
#                plot_plotmo = TRUE,
#                plotmoFilename = paste0("plotmoPlot_", this.species, "_", this.stage, ".pdf"),
#                plot_pdp = FALSE,
#                partialFilename = paste0("partialPlot_", this.species, "_", this.stage, ".pdf"))


#--------------------------------
# Stop parallel
stopCluster(cluster)
registerDoSEQ()

#--------------------------------
# End

rm(list=ls(all=TRUE))

}