##' \code{redo_sfilter()} iteratively re-fits the state-space model to tracks
##' that have failed to converge. It iteratively adjusts the values for
##' \code{span} and \code{nu} and makes \code{tries = n} attempts before giving up.
##'
##' redo ssm filter for tracks that failed to converge,
##'  trying incrementally different \code{span} & \code{nu} values.
##'  Up to \code{tries} re-filter attemps are made. If \code{common.tau = TRUE},
##'  function then searches for cases where hat(gamma) > 0.9, deployment < 30 d,
##'  device_type == PTT & refits sfilter. This may reduce over-smoothing of these
##'  short & typically sparsely observed tracks
##'
##' @title redo_sfilter
##' @param ssm_obj output from \code{sfilter()} - a \code{tibble} grouped by id or is and stage
##' @param data data object provided by \code{prefilter()}
##' @param s.inc increment for \code{span} argument
##' @param n.inc increment for \code{nu} argument
##' @param tries number of attempts to successfully filter all convergence failures before giving up
##' @param common.tau (logical) fit SSM with a single SD for measurement error in both lon & lat
##' @importFrom dplyr slice right_join select mutate ungroup group_by group_vars do
##' @export

redo_sfilter <-
  function(ssm_obj,
           data,
           s.inc = 0.1,
           n.inc = 1,
           tries = 5,
           common.tau = FALSE) {

    cf <-
      which(sapply(ssm_obj$ssm, function(x)
        length(x) == 1 || x$opt$conv != 0))
    cat(paste("\n", length(cf), " individuals failed to converge"))
    i <- 0
    spn <- 0.5
    n <- 10
    while (length(cf) > 0) {
      i <- i + 1
      spn <- spn + s.inc
      if (spn > 0.7)
        spn <- 0.7
      n = n + n.inc
      if (i == tries + 1)
        break
      cat(paste("\nattempt", i, "\n"))
      redo <- slice(ssm_obj, cf)

      if(length(group_vars(data)) == 1) {
        redo_ssm <- right_join(data, redo, by = "id") %>%
          select(-ssm) %>%
          ungroup() %>%
          group_by(id) %>%
          do(ssm = sfilter(., span = spn, nu = n))
      }
      else if(length(group_vars(data)) == 2) {
        redo_ssm <- right_join(data, redo, by = c("id", "stage")) %>%
          select(-ssm) %>%
          ungroup() %>%
          group_by(id, stage) %>%
          do(ssm = sfilter(., span = spn, nu = n))
      }
      ssm_obj$ssm[cf] <- redo_ssm$ssm
      cf <-
        which(sapply(ssm_obj$ssm, function(x)
          length(x) == 1 || x$opt$conv != 0))
      cat(paste("\n", length(cf), " individuals still not converged"))
    }
    if (common.tau) {
      ## deal with cases where gamma > 0.9, deployment < 30 days, device_type = PTT
      ##  try fitting with common.tau = TRUE
      hg <-
        which(
          sapply(ssm_obj$ssm, function(x)
            x$par["gamma", 1] > 0.9 &&
              difftime(max(x$data$date), min(x$data$date), units = "days") < 30 &&
              x$data$device_type[1] == "PTT")
        )
      cat(paste(
        "\n try estimating a common tau for lon & lat - ",
        length(hg),
        " cases\n"
      ))
      redo <- slice(ssm_obj, hg)
      if(length(group_vars(data)) == 1) {
        redo_ssm <- right_join(data, redo, by = "id") %>%
          select(-ssm) %>%
          ungroup() %>%
          group_by(id) %>%
          do(ssm = sfilter(., span = spn, nu = n))
      }
      else if( length(group_vars(data)) == 2) {
        redo_ssm <- right_join(data, redo, by = c("id", "stage")) %>%
          select(-ssm) %>%
          ungroup() %>%
          group_by(id, stage) %>%
          do(ssm = sfilter(., span = spn, nu = n))
      }
      ssm_obj$ssm[hg] <- redo_ssm$ssm
    }

    ssm_obj
  }
