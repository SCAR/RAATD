track.csv <- "data_raw_trimmed/RAATD2017_GHAL.csv"
filter.rdata <- "data_filtered/filtered_1/ghal_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for wanderers
stages <- data.frame(stage = c("Pre-breeding", "incubation", "early chick-rearing","late chick-rearing","post-breeding"),
                     start_day = c(258, 294, 352, 32, 135),
                     start_date = c("15 Sep", "21 Oct", "18 Dec", "1 Feb", "15 May"),
                     end_day = c(293, 351, 31, 134, 257),
                     end_date= c("20 Oct", "17 Dec", "31 Jan", "14 May", "15 September"))

#Buffers by device type for wchp
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
