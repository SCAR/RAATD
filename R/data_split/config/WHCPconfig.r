track.csv <- "data_raw_trimmed/RAATD2017_WHCP.csv"
filter.rdata <- "data_filtered/filtered_1/whcp_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages
stages <- data.frame(stage=c("Pre-breeding", "incubation", "early chick-rearing","late chick-rearing","post-breeding"),
                     start_day= c(250, 341, 14, 66, 104),
                     start_date= c("7 September", "7 December", "14 January", "7 March","14 April" ),
                     end_day = c(340, 13, 65, 103, 249),
                     end_date= c("6 December", "13 January", "6 March", "13 April", "6 September"))

#Buffers by device type for wchp
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
