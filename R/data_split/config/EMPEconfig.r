track.csv <- "data_raw_trimmed/RAATD2017_EMPE.csv"
filter.rdata <- "data_filtered/filtered_1/empe_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for Emperor Penguins
stages <- data.frame(stage=c("courtship", "fem_incubation", "chick-rearing", "post-breeding"),
                     start_day= c(106, 135,205,15),
                     start_date= c("16 April", "15 May", "24 July", "15 January"),
                     end_day = c(134, 204, 14, 105),
                     end_date = c("14 May", "23 July", "14 Jan", "15 April"))

#Buffers by device type 
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 200000))
