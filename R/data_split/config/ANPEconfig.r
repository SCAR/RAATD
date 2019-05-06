track.csv <- "data_raw_trimmed/RAATD2017_ANPE.csv"
filter.rdata <- "data_filtered/filtered_1/anpe_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for Antarctic Petrels
stages <- data.frame(stage = c("Pre-breeding", "incubation", "chick-rearing","post-breeding"),
                     start_day = c(289, 332, 13, 62),
                     start_date = c("16 Oct", "28 Nov", "13 Jan", "3 Mar"),
                     end_day = c(331, 12, 61, 288),
                     end_date= c("27 Nov", "12 Jan", "2 Mar", "15 Oct"))

#Buffers by device type
buff <- data.frame(type= c("GPS", "PTT", "GLS"), tbuff=c(1000, 4000, 300000))