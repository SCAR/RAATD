track.csv <- "data_raw_trimmed/RAATD2017_WAAL.csv"
filter.rdata <- "data_filtered/filtered_1/waal_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for wanderers
stages <- data.frame(stage=c("incubation", "chick-rearing","post-breeding"),
                     start_day= c(344, 59, 324),
                     start_date= c("10 December", "28 February", "25 November"),
                     end_day = c(60, 323, 343),
                     end_date= c("27 February", "24 November", "9 December"))

#Buffers by device type for 
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
