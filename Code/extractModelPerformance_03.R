## Load fitted models, extract various performance parameters, and write to file

# Ryan Reisinger
# Last modified: 2019-04-11

setwd("~/RAATD_01/RAATD/")

library(caret)
library(gbm)

meta.csv <- "/perm_storage/home/shared/github/raatd_data/metadata/SCAR_Metadata_2017_forWEBDAV.csv"
meta <- read.csv(meta.csv, stringsAsFactors = F)

species <- unique(meta$abbreviated_name)

allPerformance <- data.frame()
allImportance <- data.frame()

for (i in 1:length(species)) {
  this.species <- species[i]
  print(this.species)
  source(paste0("/perm_storage/home/shared/github/raatd_data/data_split/config/", this.species, "config.r"))
  
  resultsPerf <- data.frame()
  resultsImp <- data.frame()
  
  for (j in 1:nrow(stages)) {
    this.stage <- stages$stage[j]
    pth <- paste0(this.species, "/fittedMods/", this.species, "_", this.stage, "_gbm.RDS")
    
    if (this.species == "SOES" & this.stage == "post-moult") {
      pth <- paste0(this.species, "/fittedMods/", this.species, "_", this.stage, "_gbm_1.RDS")
    }
    if (this.species == "ANFS" & this.stage == "post-moult") {
      pth <- paste0(this.species, "/fittedMods/", this.species, "_", this.stage, "_gbm_1.RDS")
    }
    
    if (file.exists(pth)) {
      this.mod <- readRDS(pth)
      
      # Performance
      tuning <- this.mod$bestTune
      tuning$species <- this.species
      tuning$stage <- this.stage
      tuning$meanROC <- mean(this.mod$resample$ROC)
      tuning$sdROC <- sd(this.mod$resample$ROC)
      tuning$meanSens <- mean(this.mod$resample$Sens)
      tuning$meanSpec <- mean(this.mod$resample$Spec)
      tuning$nData <- nrow(this.mod$trainingData)
      tuning$nDatapoints <- length(this.mod$finalModel$data$y)
      resultsPerf <- rbind(resultsPerf, tuning)
      rm(tuning)
     
       # Variable importance
      importance <- varImp(this.mod)$importance
      importance$Variable <- row.names(importance)
      importance <- importance[order(importance$Overall, decreasing = TRUE), ]
      importance$Rank <- 1:nrow(importance)
      importance$species <- this.species
      importance$stage <- this.stage
      resultsImp <- rbind(resultsImp, importance)
      
      rm(importance)
      rm(this.mod)
    }
    
  }
  
  if (nrow(resultsPerf) > 0) {
    allPerformance <- rbind(allPerformance, resultsPerf)
    rm(resultsPerf)
    allImportance <- rbind(allImportance, resultsImp)
    rm(resultsImp)
    row.names(allImportance) <- NULL
  }
  
}

write.csv(allPerformance, "./modPerformance/modPerformance.csv", row.names = F)
write.csv(allImportance, "./modPerformance/modImportance.csv", row.names = F)

saveRDS(allPerformance, "./modPerformance/modPerformance.RDS")
saveRDS(allImportance, "./modPerformance/modImportance.RDS")

#-------------------------------------------------------------------

# Add species groups
group.1 <- c("ADPE", "ANPE", "CRAS", "EMPE", "WESE")
group.2 <- c("ANFS", "BBAL", "DMSA", "GHAL", "KIPE", "LMSA", "ROPE", "WAAL", "WHCP", "MAPE")
group.3 <- c("HUWH", "SOES")

allImportance$Group <- NA

allImportance[allImportance$species %in% group.1, "Group"] <- "Antarctic"
allImportance[allImportance$species %in% group.2, "Group"] <- "Subantarctic"
allImportance[allImportance$species %in% group.3, "Group"] <- "Both"

# Plot

library(ggplot2)

## Custom theme
theme_rr <- function () { 
  theme_bw(base_size=7, base_family="") %+replace% 
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA)
    )
}

#-------------------------------------------------------------------
## Model performance
allPerformance$Model <- paste(allPerformance$species, "-", allPerformance$stage, sep = " ")

pdf("./modPerformance/modPerformance.pdf", width= 3.5, height = 5)
p <- ggplot(data = allPerformance, aes(y = reorder(Model, meanROC, FUN = mean), x = meanROC)) + geom_point() +
  geom_segment(aes(y = Model, yend = Model, x = meanROC - sdROC, xend = meanROC + sdROC)) +
  theme_rr() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey"))
p
dev.off()

#-------------------------------------------------------------------
## Model performance v. sample size

pdf("./modPerformance/modPerformanceVSize.pdf", width= 4, height = 3.5, useDingbats = FALSE)
p <- ggplot(data = allPerformance, aes(y = meanROC, x = nDatapoints)) + geom_point() +
  geom_point() +
  geom_smooth(colour = "grey40") +
  theme_rr() +
  labs(x = "Number of data", y = "Mean ROC") +
  theme(panel.grid.major.y = element_blank())
p
dev.off()

#-------------------------------------------------------------------
## Variable importance
library(ggbeeswarm)

# All
pdf("./modPerformance/variableImportance.pdf", width= 5, height = 5)
p <- ggplot(data = allImportance, aes(x = reorder(Variable, Overall, FUN = median), y = Overall, color = Group)) +
  #geom_quasirandom(groupOnX = F) +
  geom_boxplot(outlier.shape = NA, fill = "grey", colour = "grey") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Group, ncol = 3) +
  scale_color_manual(values = c("#0077BB", "#EE7733", "#EE3377")) + guides(color = FALSE) +
  theme_rr() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) +
  coord_flip()
p
dev.off()

# Antarctic
pdf("./modPerformance/variableImportanceAntarctic.pdf", width= 3.5, height = 5)
p <- ggplot(data = allImportance[allImportance$Group == "Antarctic", ], aes(x = reorder(Variable, Overall, FUN = median), y = Overall)) +
  #geom_quasirandom(groupOnX = F) +
  geom_boxplot(outlier.shape = NA, fill = "grey", colour = "grey") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  geom_point(alpha = 0.3, color = "#0077BB") +
  theme_rr() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) +
  coord_flip()
p
dev.off()

# Subantarctic
pdf("./modPerformance/variableImportanceSubantarctic.pdf", width= 3.5, height = 5)
p <- ggplot(data = allImportance[allImportance$Group == "Subantarctic", ], aes(x = reorder(Variable, Overall, FUN = median), y = Overall)) +
  #geom_quasirandom(groupOnX = F) +
  geom_boxplot(outlier.shape = NA, fill = "grey", colour = "grey") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  geom_point(alpha = 0.3, color = "#EE3377") +
  theme_rr() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) +
  coord_flip()
p
dev.off()

# Both
pdf("./modPerformance/variableImportanceBoth.pdf", width= 3.5, height = 5)
p <- ggplot(data = allImportance[allImportance$Group == "Both", ], aes(x = reorder(Variable, Overall, FUN = median), y = Overall)) +
  #geom_quasirandom(groupOnX = F) +
  geom_boxplot(outlier.shape = NA, fill = "grey", colour = "grey") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  geom_point(alpha = 0.3, color = "#EE7733") +
  theme_rr() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) +
  coord_flip()
p
dev.off()



