##' \code{sfilter()} fits a simple state-space model to fit to pre-filtered RAATD data and returns output
##'  as a tibble grouped by individual id or individual id and breeding stage.
##'
##' Wrapper function that takes data prepared by \code{prefilter()}, assigns appropriate
##' time steps for GPS, GLS and/or PTT devices & calls \code{ssmTMB::fit_ssm()} to do the
##' filtering. \code{ssmTMB} can be installed via \code{devtools::install_github("ianjonsen/ssmTMB")}.
##' Note, as \code{sfilter()} imports functions from the \code{ssmTMB} package, you must ensure you
##' have installed the \code{TMB} package and it's dependencies
##'
##' @title State-space filter RAATD tracking data
##' @param d input data from \code{prefilter()} as a tibble of individual tracks grouped by id
##' @param ts specify a list of time steps (in h) for gps, gls, & ptt datasets
##' @param ... additional arguments passed to \code{ssmTMB::fit_ssm()}. Two common arguments are
##' \code{span}, the degree of loess smoothing used to estimate initial values for the location states, and
##' \code{nu}, the degree of freedom to be used in the t-distributed error models for lon & lat
##' @return a list with 10 elements (see \code{?ssmTMB::fit_ssm()})
##'
##' @examples
##' \dontrun{
##' ## run prefilter with example royal penguin data & metadata
##' data(rope)
##'
##' pfd <- prefilter(
##'   dat = rope,
##'   metadata = meta,
##'   sp = "ROPE",
##'   min_obs = 30,
##'   min_days = 5,
##'   vmax = 10
##'   )
##'
##' ssm_by_id <- pfd %>%
##'   do(ssm = sfilter(., span = 0.4, nu = 5))
##'
##' ## re-filter tracks that failed to converge
##' ssm_by_id <- redo_sfilter(ssm_by_id, pfd, tries = 10)
##'
##' ## generate QC plots
##' ssm_by_id %>% qc_plot(sp = "rope")
##' }
##'
##' @importFrom dplyr mutate
##' @importFrom ssmTMB fit_ssm
##' @export

sfilter <- function(d,
                   ts = list(gps = 1, gls = 12, ptt = 2),
                   ...) {

  switch(unique(d$device_type),
         GPS = {
             ts = ts$gps
         },
         GLS = {
             ts = ts$gls
         },
         PTT = {
             ts = ts$ptt
         })

  ## when track straddles -180,+180 shift to 0,360
  f <- subset(d, keep)

  if(diff(range(f$lon)) > 300) {
    tmp <- f %>% mutate(lon = wrapLon(lon, lmin = 0))
    if(min(diff(tmp$lon)) < -310 || max(diff(d$lon)) > 310) {
      d[d$keep, ] <- f %>% mutate(lon = unwrapLon(lon, lmin = 0))
    }
    else{
      d[d$keep, ] <- tmp
    }
  }
options(warn = -1) ## prevents warnings from TMB re: NaNs in sqrt(diag(cov))
  out <- try(fit_ssm(d, subset = d$keep, tstep = ts / 24, ...), silent = TRUE)
options(warn = 0)

  ## if track on 0,360+ then shift back to -180,+180
  if(length(out) > 1) {
    out$predicted <- out$predicted %>% mutate(lon = wrapLon(lon, lmin = -180))
    out$fitted <- out$fitted %>% mutate(lon = wrapLon(lon, lmin = -180))
    out$data <- out$data %>% mutate(lon = wrapLon(lon, lmin = -180))
  }

  out
}

