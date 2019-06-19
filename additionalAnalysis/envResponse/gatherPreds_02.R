## Gather data to plot single response curves

setwd("~/RAATD_01/RAATD/additionalAnalysis/envResponse/")

library(raster)

# Get the mean environmental data
env <- readRDS("~/RAATD_01/RAATD/meanPredictions/griddedScoresEnvarsUntransformed.RDS")
env <- env[ , c(1:2, 19:37)]

# List the rasters with gbm predictions
fls <- list.files("../../", full.names = T, recursive = T)
fls <- fls[grep("_gbm_rast.gri", fls)]

# Get species names
dx <- data.frame("files" = fls,
                 "sp" = substr(fls, 8, 11))

# Extract values for each species * stage,
# then take the mean for each species
sps <-unique(dx$sp)

# Don't run if redoing plots:

#---------------------
# bar <- env[ , c("x", "y")]
# 
# for (i in 1:length(sps)) {
#   this.sp <- sps[i]
#   print(this.sp)
#   
#   these.files <- dx[dx$sp == this.sp, "files"]
#   
#   foo <- env[ , c("x", "y")]
#   
#   # Get the gbm prediction from each file for given species
#   for (j in 1:length(these.files)) {
#     this.file <- raster(paste(these.files[j]))
#     hld <- extract(this.file, foo[ , c("x", "y")])
#     foo <- cbind(foo, hld)
#     rm(hld)
#   }
#   
#   if (ncol(foo) > 3) {
#     foo <- rowMeans(foo[ , -(1:2)])
#   } else {
#     foo <- foo[ , 3]
#   }
#   
#   bar <- cbind(bar, foo)
#   rm(foo)
#   
# }
# 
# # Re-name
# names(bar) <- c("x", "y", paste(sps))
# 
# # And add the envars
# bar <- cbind(env, bar[ , -(1:2)])
# 
# # Save
# saveRDS(bar, "untransformedPredsForResponsePlot.RDS")
#---------------------
# end of don't run

bar <- readRDS("untransformedPredsForResponsePlot.RDS")
#-----------------------------------------------------
## Plotting

library(ggplot2)
library(tidyr)
library(mgcv)


# library(extrafont)
# font_import()

# Transform data to long format
newbar <- gather(bar, species, score, ADPE:WHCP)
newbar <- gather(newbar, covariate, value, CHLA:WIND)

# Loop through covariates, plottin each
covars <- unique(newbar$covariate)

for (i in covars) {
  
  # i <- "CHLA"
  
  print(i)
  
  hld <- newbar[newbar$covariate == i, ]
  
  png(paste0(i, ".png"), height = 9/5, width = 7.2/4, units = "in", res = 300, family = "Arial")
  
  p <- ggplot(hld, aes(x = value, y = score, group = species)) +
    geom_smooth(colour = "black",
                size = 0.3,
                method = "bam", formula = y ~ s(x, bs = "tp"), se = FALSE) +
    geom_rug(size = 0.1, alpha = 0.1, colour = "gray10", sides = "b") +
    facet_wrap(.~covariate, ncol = 1, nrow = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    expand_limits(y = c(0, 1)) +
    theme_linedraw() +
    theme(text = element_text(family = "Arial",
                              size = 7),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          axis.title = element_blank())
  
  print(p)
  
  dev.off()
  
}