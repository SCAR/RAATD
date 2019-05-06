track.csv <- "data_raw_trimmed/RAATD2017_BBAL.csv"
filter.rdata <- "data_filtered/filtered_1/bbal_ssm_by_id.RDS"

bound <- c(-180,180,-80,-40)

#Pre-defined stages
stages <- data.frame(stage=c("arrival", "incubation", "chick-rearing","post-breeding"),
                     start_day= c(245, 275, 346, 107),
                     start_date= c("1 Sept", "1 Oct", "11 Dec", "15 Apr"),
                     end_day = c(274, 345, 106, 244),
                     end_date = c("30 Sept", "10 Dec", "14 Apr", "31 Aug"))

#Buffers by device
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
