##' Prefilter
##'
##' Prepare standardised track data for state-space filtering by performing
##' further standardisation and culling unacceptably short deployments
##'
##' The track data are fetched automatically from the RAATD repository by specifying
##' a species' abbreviated name.
##'
##' @title Prefilter
##' @param dat provide example royal penguin tracking data
##' @param metadata provide example royal penguin metadata
##' @param min_obs minimum number of observation records a deployment must have
##' @param min_days minimum number days a deployment must last
##' @param vmax maximum travel rate (m/s) for speed filter
##' @param sp for future functionality, not currently used
##' @param path2repo for future functionality, not currently used
##' @return writes a tibble to a .csv file with the following variables:
##' \item{\code{id}}{individual animal id}
##' \item{\code{date}}{POSIX date-time}
##' \item{\code{lc}}{location quality class}
##' \item{\code{lon}}{longitude}
##' \item{\code{lat}}{latitude}
##' \item{\code{device_type}}{tag type: GPS, PTT, or GLS}
##' \item{\code{keep}}{should record be used of ignored by SSM filter (logical)}
##'
##' @examples
##' \dontrun{
##' ## prefilter RAATD royal penguin tracks
##' ## load example data
##' data(rope)
##'
##' pfd <- prefilter(
##'   dat = rope,
##'   metadata = meta,
##'   min_obs = 30,
##'   min_days = 5,
##'   vmax = 10
##'   )
##' }
##'
##' @importFrom diveMove grpSpeedFilter
##' @importFrom lubridate ymd_hms
##' @importFrom readr read_csv cols col_character col_integer col_time col_double
##' @importFrom dplyr filter mutate select rename left_join group_by distinct arrange do
##' @importFrom tibble tibble
##' @export

prefilter <-
  function(dat,
           metadata,
           min_obs,
           min_days,
           vmax,
           sp = NULL,
           path2repo = NULL
           ) {
    if(is.null(dat)) stop("An input data frame must be supplied")
    if(is.null(metadata)) stop("Input metadata must be supplied")
    if(is.null(sp)) stop("A RAATD species 4-letter abbreivated name must be supplied")

    ## Functionality removed until RAATD data are published
    # ## load canonical metadata & tracking data files
    # if(is.null(fullpath)) {
    #   fp_meta <-
    #     file.path(
    #       path2repo,
    #       "raatd_data",
    #       "metadata",
    #       "SCAR_Metadata_2017_forWEBDAV.csv")
    #
    #   fp_raatd <-
    #     file.path(
    #       path2repo,
    #       "raatd_data",
    #       "data_raw_trimmed",
    #       paste("RAATD2017_", sp, ".csv", sep = "")
    #     )
    # }
    # else {
    #   fp_meta <- fullpath[1]
    #   fp_raatd <- fullpath[2]
    # }
    # ## merge metadata & subset to species == sp
    #
    #
    # m_sub <- read_csv(fp_meta) %>%
    #   filter(., abbreviated_name == sp)

    m_sub <- metadata
    dev_dat <- m_sub %>%
      select(individual_id, device_type) %>%
      rename(id = individual_id)

    ## find individual tracks to discard
    discard_ids <- m_sub %>%
      filter(keepornot == "Discard") %>%
      select(individual_id)

    ## Functionality removed until RAATD data are published
    # ## read in tracking dataset & merge with dev.dat to add correct Unit (tag/data) type
    # dat <-
    #   read_csv(
    #     fp_raatd,
    #     col_types = cols(
    #       individual_id = col_character(),
    #       breeding_stage = col_character(),
    #       year = col_integer(),
    #       month = col_integer(),
    #       day = col_integer(),
    #       time = col_time(format = ""),
    #       decimal_latitude = col_double(),
    #       decimal_longitude = col_double(),
    #       location_quality = col_character(),
    #       location_to_keep = col_integer()
    #     ),
    #     na = c("", "NA", "NaN")
    #   )

    ## remove 'hanging' locations at track starts & ends using
    ##    'location_to_keep' flag in canonical data this preserves
    ##     the "trimmed data" and discards the "trimmings"
    d <- dat %>%
      filter(location_to_keep == 1 |
               is.na(location_to_keep)) %>%
      mutate(date = ymd_hms(paste(paste(
        year, month, day, sep = "-"
      ), time), tz = "UTC")) %>%
      rename(id = individual_id) %>%
      rename(lc = location_quality) %>%
      rename(lon = decimal_longitude) %>%
      rename(lat = decimal_latitude) %>%
      left_join(., dev_dat, by = "id") %>%
      select(id, date, lc, lon, lat, device_type)

    ## remove individual tracks labelled as "Discard"
    d <- d %>%
      filter(!id %in% discard_ids$individual_id)

    ## assign lc B to all GLS data and to all PTT data where lc = NA
    ## assign lc 3 to all GPS data
    ## deal with common alternate PTT lc designations
    ## coerce all lc Z to lc B
    d <- d %>%
      mutate(lc = ifelse(is.na(lc) &
                           device_type != "GPS", "B", lc)) %>%
      mutate(lc = ifelse(is.na(lc) &
                           device_type == "GPS", "3", lc)) %>%
      mutate(lc = ifelse(lc == "-9" &
                           device_type == "PTT", "Z", lc)) %>%
      mutate(lc = ifelse(lc == "-3" &
                           device_type == "PTT", "Z", lc)) %>%
      mutate(lc = ifelse(lc == "-2" &
                           device_type == "PTT", "B", lc)) %>%
      mutate(lc = ifelse(lc == "-1" &
                           device_type == "PTT", "A", lc)) %>%
      mutate(lc = ifelse(lc == "Z", "B", lc))

    ## mutate lc to ordered factor
    d <- d %>%
      mutate(lc = factor(lc,
                         levels = c(3, 2, 1, 0, "A", "B"),
                         ordered = TRUE))

    d <- d %>%
      group_by(id)

    ## remove any duplicated time records, any nearly duplicate time records occur within 120 s
    ##  & order records by time within each individual track
    ##  keep track of # records removed
    d1 <- d %>%
      do(distinct(., date, .keep_all = TRUE)) %>%
      do(mutate(
        .,
        dup = difftime(date, lag(date), units = "secs") < 120 &
          device_type == "PTT"
      )) %>%
      do(arrange(., order(date))) %>%
      filter(.,!dup) %>%
      select(-dup)

    d1.rep <- nrow(d) - nrow(d1)
    cat(sprintf("%d duplicate time &/or near duplicate location/time records removed\n", d1.rep))

    ## remove locations in N hemisphere, eg. when tags turned on by manufacturer
    d2 <- d1 %>%
      do(filter(., lat < 10))
    d2.rep <- nrow(d1) - nrow(d2)
    cat(sprintf("%d records at > 10 Latitude N removed\n", d2.rep))

    ## remove deployments with less than min_obs
    d3 <- d2 %>%
      do(filter(., n() >= min_obs))
    d3.rep <- n_groups(d2) - n_groups(d3)
    cat(sprintf(
      paste(
        "%d individuals with fewer than",
        min_obs,
        "observations removed\n"
      ),
      d3.rep
    ))

    ## remove GPS deployments that last less than 1 day, other deployments les than min_days days
    f <- d3 %>%
      do(filter(
        .,
        device_type != "GPS",
        difftime(max(date), min(date), units = "days") >= min_days
      ))
    g <- d3 %>%
      do(filter(
        .,
        device_type == "GPS",
        difftime(max(date), min(date), units = "days") >= 1
      ))
    d4 <- rbind(g, f)
    d4.rep <- n_groups(d3) - n_groups(d4)
    cat(sprintf(
      paste(
        "%d individuals with fewer than",
        min_days,
        "deployment days removed\n"
      ),
      d4.rep
    ))

    ## flag extreme travel rate locations to be ignored at ssm filter stage
    d5 <- d4 %>%
      do(mutate(., keep = grpSpeedFilter(cbind(date, lon, lat), speed.thr = vmax)))
    d5.rep <- nrow(d4) - sum(d5$keep)
    cat(sprintf(
      paste(
        "\n%d records with travel rates >",
        vmax,
        "m/s will be ignored by SSM filter\n"
      ),
      d5.rep
    ))

    ## return pre-filtered data
    ind <- c(NA, NA, d3.rep, d4.rep, NA)
    ign <- c(NA, NA, NA, NA, d5.rep)

    ## Functionality removed until RAATD data are published
    # rep <- tibble(
    #   test = c(
    #     "duplicated_dates",
    #     "lat_above_10N",
    #     paste("obs_lessthan_", min_obs, sep = ""),
    #     paste("deploy_lessthan_", min_days, "_days", sep =
    #             ""),
    #     paste("obs_greaterthan_", vmax, "_ms", sep =
    #             "")
    #   ),
    #   obs_removed = c(d1.rep, d2.rep, NA, NA, NA),
    #   individuals_removed = c(NA, NA, d3.rep, d4.rep, NA),
    #   obs_flagged_to_ignore = c(NA, NA, NA, NA, d5.rep)
    # )
    #
    #   fp_rep <- file.path(
    #     path2repo,
    #     "raatd_data",
    #     "data_filtered",
    #     "data_prefiltered",
    #     paste(sp, "_prefilter_report.csv", sep = "")
    #   )
    #
    #   fp_pfout <- file.path(
    #     path2repo,
    #     "raatd_data",
    #     "data_filtered",
    #     "data_prefiltered",
    #     paste(sp, "_prefiltered_data.csv", sep = "")
    #   )
    #
    #  write_csv(rep, fp_rep)
    #  write_csv(d5, fp_pfout)

    d5
  }

