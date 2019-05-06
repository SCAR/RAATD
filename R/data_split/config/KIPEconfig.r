track.csv <- "data_raw_trimmed/RAATD2017_KIPE.csv"
filter.rdata <- "data_filtered/filtered_1/kipe_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for wanderers
stages <- data.frame(stage = c("Pre-breeding", "incubation", "chick-rearing"),
                     start_day = c(305, 356, 53),
                     start_date = c("1 November", "22 December", "22 February"),
                     end_day = c(355, 52, 304),
                     end_date= c("21 December", "21 February", "31 October"))

#Buffers by device type for wchp
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
