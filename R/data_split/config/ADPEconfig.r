track.csv <- "data_raw_trimmed/RAATD2017_ADPE.csv"
filter.rdata <- "data_filtered/filtered_1/adpe_ssm_by_id.RDS"

bound <- c(0,360,-80,-40)

#Pre-defined stages for Adelies
stages <- data.frame(stage = c("arrival", "incubation", "chick-rearing",  "pre-moult", "post-moult"),
                     start_day = c(289, 320, 354, 46, 92),
                     start_date = c("15 Oct", "15 Nov", "19 Dec", "15 Feb", "1 Apr"),
                     end_day = c(319, 353, 45, 91, 288),
                     end_date = c("14 Nov", "18 Dec", "14 Feb", "31 Mar", "14 Oct"))

#Buffers by device type for Adelies
buff <- data.frame(type = c("GPS", "PTT", "GLS"), tbuff = c(1000, 4000, 200000))
