track.csv <- "data_raw_trimmed/RAATD2017_WESE.csv"
filter.rdata <- "data_filtered/filtered_1/wese_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#All WESE tracks are post-moult
stages <- data.frame(stage = c("no-stage"),
                     start_day = c(1),
                     start_date = c("1 January"),
                     end_day = c(365),
                     end_date = c("31 December"))

#Buffers by device type
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
