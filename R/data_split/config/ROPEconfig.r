track.csv <- "data_raw_trimmed/RAATD2017_ROPE.csv"
filter.rdata <- "data_filtered/filtered_1/rope_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for wanderers
stages <- data.frame(stage = c("Pre-breeding", "incubation", "early chick-rearing", "post-breeding"),
                     start_day = c(244, 288, 335, 54),
                     start_date = c("1 Sep", "15 Oct", "1 Dec", "23 Feb"),
                     end_day = c(287, 334, 53, 243),
                     end_date= c("14 Oct", "30 November", "22 Feb", "31 Aug"))

#Buffers by device type for wchp
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))
