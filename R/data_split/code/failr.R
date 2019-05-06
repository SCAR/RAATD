# Function which removes individuals which failed quality control from SSM fits

# Ryan Reisinger
# Last modified 2017-08-22

# 'fits' is an SSM fit object
# returns an object with the same class and structure, but with failures removed

failr <- function(fits){

fits.new <- list(id = NULL, ssm = NULL)
fits.new$id <- fits$id[which(fits$ssm != "failed QC round 1" & fits$ssm != "failed to converge")]
fits.new$ssm <- fits$ssm[which(fits$ssm != "failed QC round 1" & fits$ssm != "failed to converge")]
  return(fits.new)
}