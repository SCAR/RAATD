track.csv <- "data_raw_trimmed/RAATD2017_LMSA.csv"
filter.rdata <- "data_filtered/filtered_1/lmsa_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for wanderers
stages <- data.frame(stage = c("Pre-breeding", "incubation", "early chick-rearing","late chick-rearing","post-breeding"),
                     start_day = c(274, 305, 1, 39, 152),
                     start_date = c("1 October", "1 November", "1 January", "8 February", "1 June"),
                     end_day = c(304, 365, 38, 151, 273),
                     end_date= c("31 October", "31 December", "7 February", "31 May", "30 September"))

#Buffers by device type for wchp
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
