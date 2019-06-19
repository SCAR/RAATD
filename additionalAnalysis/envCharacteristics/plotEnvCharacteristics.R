# Create density plots of used v. available habitat

setwd("~/RAATD_01/RAATD/additionalAnalysis/envCharacteristics/")

library(ggplot2)
library(tidyr)

source("/perm_storage/home/shared/github/raatd_modelling/caret/prepCaret3.R")
source("~/RAATD_01/RAATD/Code/function_splitViolin.R")

# Plotting theme
theme_ryan <- function () { 
  theme_bw(base_size=9, base_family="") %+replace% 
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_blank(),
      axis.ticks = element_line(colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.line = element_line(colour = "black")
    )
}

#--------------------------------------------------------
# Loop through species

these.species <- c("ADPE",
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

for (i in 1:length(these.species)) {
  
  # this.species <- "ADPE"
  this.species <- these.species[i]
  print(this.species)
  
  #-------------------
  # Get filenames
  fls <- list.files(paste0("~/RAATD_01/RAATD/", this.species), full.names = T)
  fls <- fls[grepl("-usage-", fls)]
  
  # For each stage
  for (j in 1:length(fls)) {
    
    # Load usage
    load(fls[j])
    
    # Coerce the tbl back to data.frame
    usage <- as.data.frame(usage)
    
    # Expand the dataframe
    dat <- prepCaret(usage = usage, toobig = FALSE, check.na = FALSE)
    
    # Get stage name
    this.stage <- unique(dat$stage)
    
    # Ignore dropped models
    run <- TRUE
    
    if (this.species == "ROPE" & this.stage == "early chick-rearing") {
      run <- FALSE
    }
    
    if (this.species == "WHCP" & this.stage == "late chick-rearing") {
      run <- FALSE
    }
    
    if (this.species == "GHAL" & this.stage == "Pre-breeding") {
      run <- FALSE
    }
    
    # Run if it's not a dropped model
    if (run) {
      
      # Add a dummy variable for labels
      dat$lab <- "Available | Used"
      
      # Check if CHLA is in output
      # then go to long
      if (length(which(grepl("CHLA", colnames(dat)))) == 0) {
        
        tst <- gather(data = dat,
                      key = "Covariate",
                      value = "Value",
                      "CURR",
                      "DEPTH",
                      "DEPTHg",
                      "dSHELF",
                      "EKE",
                      "ICE",
                      "ICEA",
                      "ICEsd",
                      "SAL",
                      "SHFLUX",
                      "SHFLUXsd",
                      "SSHa",
                      "SSHsd",
                      "SST",
                      "SSTg",
                      "VMIX",
                      "VMIXsd",
                      "WIND")
      } else {
        tst <- gather(data = dat,
                      key = "Covariate",
                      value = "Value",
                      "CHLA",
                      "CURR",
                      "DEPTH",
                      "DEPTHg",
                      "dSHELF",
                      "EKE",
                      "ICE",
                      "ICEA",
                      "ICEsd",
                      "SAL",
                      "SHFLUX",
                      "SHFLUXsd",
                      "SSHa",
                      "SSHsd",
                      "SST",
                      "SSTg",
                      "VMIX",
                      "VMIXsd",
                      "WIND")
      }
      
      # Create a filename
      fn <- paste0("./Plots/",
                   this.species,
                   "_",
                   this.stage,
                   ".pdf")
      
      # Plot to file
      pdf(fn, paper = "a4")
      p <- ggplot(data = tst, aes(x = lab, y = Value, fill = used)) +
        geom_split_violin() +
        scale_fill_manual(values = c("grey", "#0077BB"), guide = FALSE) +
        facet_wrap(~Covariate, ncol = 3, scales = "free_y") +
        labs(title = this.species, subtitle = this.stage) +
        theme_ryan()
      print(p)
      dev.off()
      
    }
    
    rm(usage)
    
  }
  
}
