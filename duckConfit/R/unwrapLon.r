##' unwrapLon
##'
##' Internal function not normally called by user
##'
##' @title unwrapLon
##' @export

unwrapLon <- function(lon, lmin = -180)
  cumsum(c(wrapLon(lon[1], lmin), wrapLon(diff(lon))))
