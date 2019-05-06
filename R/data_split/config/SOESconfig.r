track.csv <- "data_raw_trimmed/RAATD2017_SOES.csv"
filter.rdata <- "data_filtered/filtered_1/soes_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for southern elephant seals
stages <- data.frame(stage = c("post-breeding", "post-moult"),
                     start_day = c(284, 31),
                     start_date = c("11 Oct", "1 Feb"),
                     end_day = c(30, 283),
                     end_date = c("31 Jan", "10 Oct"))

#Buffers by device type 
buff <- data.frame(type = c("GPS", "PTT", "GLS"), tbuff = c(1000, 5000, 200000))
