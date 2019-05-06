track.csv <- "data_raw_trimmed/RAATD2017_MAPE.csv"
filter.rdata <- "data_filtered/filtered_1/mape_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for wanderers
stages <- data.frame(stage = c("Pre-breeding", "incubation", "early chick-rearing","late chick-rearing","post-breeding"),
                     start_day = c(288, 335, 1, 32, 60),
                     start_date = c("15 October", "1 December", "1 January", "1 February", "1 March"),
                     end_day = c(334, 365, 31, 59, 287),
                     end_date= c("6 December", "31 December", "31 January", "28 February", "14 October"))

#Buffers by device type for wchp
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
