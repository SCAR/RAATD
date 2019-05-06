track.csv <- "data_raw_trimmed/RAATD2017_ANFS.csv"
filter.rdata <- "data_filtered/filtered_1/anfs_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for Antarctic fur seals
stages <- data.frame(stage = c("breeding", "post-moult"),
                     start_day = c(335, 91),
                     start_date = c("1 Dec", "31 Mar" ),
                     end_day = c(92, 334),
                     end_date = c("1 Apr", "30 Nov"))

#Buffers by device type 
buff <- data.frame(type = c("GPS", "PTT", "GLS"), tbuff = c(1000, 4000, 50000))
