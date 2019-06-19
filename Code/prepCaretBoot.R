# 'prepCaretBoot' function
# Prepares usage data for bootstrap modelling in caret package

# input:
# 'usage' = usage dataframe produced by the track simulation script

# Ryan Reisinger

# Last modified: 2018-01-29

prepCaretBoot <- function(usage, sze = 0.5){
  
  require(splitstackshape)
  
  #--------------------------------------------------------
  # The usage data is in the form of one row per cell,
  # it needs to be expanded row-wise for modelling
  # The expansion will then allow us to model:
  # given a cell was available, was it used (1) or not (0)
  
  #-------------------------------------------------------
  
  # Add row ids
  usage$idx <- 1:nrow(usage)
  
  # Expand each component separately
  usage$howmany <- usage$avail - usage$usage
  
  datAvail <- expandRows(dataset = usage, count = "howmany", drop = TRUE)
  datAvail$used <- "N"
  
  usage$howmany <- NULL
  
  datUse <- expandRows(dataset = usage, count = "usage", drop = FALSE)
  datUse$used <- "Y"
  
  # Omit NA values
  datAvail <- datAvail[complete.cases(datAvail), ]
  datUse <- datUse[complete.cases(datUse), ]
  
  # Sample with replacement (separately from used and available)
  datAvail <- datAvail[sample(nrow(datAvail), nrow(datAvail)*sze, replace = TRUE), ]
  datUse <- datUse[sample(nrow(datUse), nrow(datUse)*sze, replace = TRUE), ]
  
  dat <- rbind(datAvail, datUse)
  
  dat$used <- as.factor(dat$used)
  
  rm(datAvail, datUse)
  
  return(dat)
}