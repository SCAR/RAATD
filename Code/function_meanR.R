# Function calculates the mean of species in predefined groups,
# and then calculates the maximum of these group means for each cell

# Ryan Reisinger
# September 2018
# BR, Nov 2018 replace `these.groups=="ALL"` with `identical(these.groups, "ALL")`

meanR <- function(this.data, these.groups) {
  
  # this.data = a data.frame containing species columns, function returns this data.frame
  # with an additional column named MEAN (for compatibility, it has this uninformative name)
  
  # these.groups = a list of vectors of the species codes in each group,
  # OR 'ALL', which calculates the mean for all species.

  all_groups <- identical(these.groups, "ALL")
  
  if (all_groups) {
    print("Using all species")
  } else if (!is.list(these.groups)) {
    stop("these.groups must be a list of vectors")
  } else {
    sprintf("Using %s group(s)", length(these.groups))
  }
  
  if (all_groups) {
    these.groups <- list(c("ADPE",
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
                           "MAPE.ROPE",
                           "SOES",
                           "WAAL",
                           "WESE",
                           "WHCP"))
  }
  
  foo <- list()

  for (i in 1:length(these.groups)) {
    these.species <- unlist(these.groups[i])
    bar <- this.data[ , c(these.species)]
    foo[i] <- list(rowMeans(bar))
  }
  
  foo <- as.data.frame(do.call(cbind, foo))
  
  this.data$MEAN <- do.call(pmax, foo)
  
  return(this.data)
  
}
