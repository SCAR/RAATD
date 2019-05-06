track.csv <- "data_raw_trimmed/RAATD2017_DMSA.csv"
filter.rdata <- "data_filtered/filtered_1/dmsa_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for wanderers
stages <- data.frame(stage = c("Pre-breeding", "incubation", "early chick-rearing","late chick-rearing","post-breeding"),
                     start_day = c(244, 274, 15, 53, 151),
                     start_date = c("1 Sep", "1 Oct", "15 Jan", "22 Feb", "1 Jun"),
                     end_day = c(273, 14, 52, 150, 243),
                     end_date= c("30 Sep", "14 Jan", "21 Feb", "31 May", "31 Aug"))

#Buffers by device type for wchp
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
