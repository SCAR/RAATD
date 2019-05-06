# Manual tweaks and corrections of automatic stage date calculations

setwd("D:/RAATD/Data/GitHubClones/raatd_data")

stageMerger <- function(dat, id, becomes) {
  d <- dat[dat$id == id, ]
  newstart <- min(d$start)
  newend <- max(d$end)
  newrow <- data.frame(id = id, stage = becomes, start = newstart, end = newend)
  hold <- dat[dat$id != id, ]
  hold <- rbind(hold, newrow)
  return(hold)
}

# Double check that there are no QC1 fails in the list

fails <- read.csv("./data_filtered/data_for_QC/qcNotes.csv", stringsAsFactors = F)
fails <- fails[fails$qc_round == "qc1" & fails$response_action == "discard", "individual_id"]

#--------------------
#ADPE

ADPE <- readRDS("./data_split/stage_dates/ADPE_stages.RDS")

#Fix stage assignments
ADPE[ADPE$id == "DDU IDB GPS107 31122011", "stage"] <- "incubation"
ADPE[ADPE$id == "DDU IDE GPS126 30122011", "stage"] <- "incubation"
ADPE[ADPE$id == "H11CROZ0809", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "H3CROZ0809", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "H4CROZ0809", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "H6CROZ0506", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "H6CROZ0607", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "H7CROZ0607", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "H9CROZ0607", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "S4BIRD0203", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "unknown208", "stage"] <- "incubation"
ADPE[ADPE$id == "unknown250", "stage"] <- "pre-moult"
ADPE[ADPE$id == "unknown255", "stage"] <- "pre-moult"
ADPE[ADPE$id == "unknown289", "stage"] <- "pre-moult"
ADPE[ADPE$id == "ADPE-C0-RAATD", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "ADPE-C7-RAATD", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE1-RAATD", "stage"] <- "incubation"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE11-RAATD", "stage"] <- "incubation"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE12-RAATD", "stage"] <- "chick-rearing"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE15-RAATD", "stage"] <- "incubation"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE17-RAATD", "stage"] <- "pre-moult"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE19-RAATD", "stage"] <- "incubation"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE3-RAATD", "stage"] <- "pre-moult"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE9-RAATD", "stage"] <- "incubation"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE73-RAATD", "stage"] <- "incubation"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE60-RAATD", "stage"] <- "incubation"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE50-RAATD", "stage"] <- "incubation"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE5-RAATD", "stage"] <- "incubation"
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE4-RAATD", "stage"] <- "incubation"

#Merge stages
birds <- list("H11CROZ0607",
              "H1BIRD0607",
              "H3BIRD0607",
              "H5CROZ0506",
              "S1CROZ0001",
              "S2BIRD0203",
              "S2BIRD0304",
              "S3CROZ0203",
              "S5CROZ0203",
              "unknown253",
              "ADPE-D1-RAATD",
              "ADPE-dtsetBirdLife773-4-RAATD")
              
for (i in 1:length(birds)) {
  whichbird = birds[[i]]
  ADPE <- stageMerger(dat = ADPE, id = whichbird, becomes = "chick-rearing")
}

#Split stages

# ADPE-dtsetBirdLife910-ADPE24-RAATD
# last trip should be pre-moult
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE24-RAATD", ]
newrow <- data.frame(id = "ADPE-dtsetBirdLife910-ADPE24-RAATD",
                     stage = "pre-moult",
                     start = "2003-01-24 18:02:00",
                     end = "2003-03-08 05:02:00")
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE24-RAATD" & ADPE$stage == "chick-rearing", "end"] <- "2003-01-24 18:02:00"
ADPE <- rbind(ADPE, newrow)

# ADPE-dtsetBirdLife910-ADPE29-RAATD
# last trip should be pre-moult
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE29-RAATD", ]
newrow <- data.frame(id = "ADPE-dtsetBirdLife910-ADPE29-RAATD",
                     stage = "pre-moult",
                     start = "2005-01-14 23:51:00",
                     end = "2005-02-01 18:51:00")
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE29-RAATD" & ADPE$stage == "chick-rearing", "end"] <- "2005-01-14 23:51:00"
ADPE <- rbind(ADPE, newrow)

# ADPE-dtsetBirdLife910-ADPE31-RAATD
# last trip should be pre-moult
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE31-RAATD", ]
newrow <- data.frame(id = "ADPE-dtsetBirdLife910-ADPE31-RAATD",
                     stage = "pre-moult",
                     start = "2004-01-17 08:49:00",
                     end = "2004-02-08 13:49:00")
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE31-RAATD" & ADPE$stage == "chick-rearing", "end"] <- "2004-01-17 08:49:00"
ADPE <- rbind(ADPE, newrow)

# Kane
# last trip should be pre-moult
ADPE[ADPE$id == "Kane", ]
newrow <- data.frame(id = "Kane",
                     stage = "pre-moult",
                     start = "2003-01-15 21:34:43",
                     end = "2003-02-16 14:34:43")
ADPE[ADPE$id == "ADPE-dtsetBirdLife910-ADPE31-RAATD" & ADPE$stage == "chick-rearing", "end"] <- "2003-01-15 21:34:43"
ADPE <- rbind(ADPE, newrow)

ADPE <- ADPE[!(ADPE$id %in% fails), ]

saveRDS(ADPE, "./data_split/stage_dates/ADPE_stages.RDS")


#--------------------
#ANFS

ANFSdates <- readRDS("./data_split/stage_dates/ANFS_stages.RDS")

#'FM5' should all be post-moult
ANFSdates[ANFSdates$id == "FM5", ]
ANFSdates[ANFSdates$id == "FM5" & ANFSdates$stage == "post-moult", "start"] <- "2003-12-27 23:45:06"
ANFSdates <- ANFSdates[!(ANFSdates$id == "FM5" & ANFSdates$stage == "breeding"), ]

ANFSdates <- ANFSdates[!(ANFSdates$id %in% fails), ]

saveRDS(ANFSdates, "./data_split/stage_dates/ANFS_stages.RDS")


#--------------------
#BBAL

BBAL <- readRDS("./data_split/stage_dates/BBAL_stages.RDS")

#Fix Kerguelen post-breeding
for (i in 1:nrow(BBAL)) {
  if (substr(BBAL$id[i], 1, 14) == "BBAL-Kerguelen" & BBAL$stage[i] == "post-breeding") {
    yr <- format(BBAL$end[i], "%Y")
    tm <- paste0(yr, "-09-02 00:00:00")
    BBAL[i, "end"] <- as.character(tm)
  }
}

#Fix hanging stages
for (i in 1:nrow(BBAL)) {
  if (substr(BBAL$id[i], 1, 14) == "BBAL-Kerguelen" & BBAL$stage[i] == "incubation") {
    yr <- format(BBAL$start[i], "%Y")
    tm <- paste0(yr, "-10-12 00:00:00")
    strt <- BBAL[i, "start"]
    if (strt < tm) {
      strt <- tm
      BBAL[i, "start"] <- as.character(strt)
    }
  }
}

#Fixes for non-Kerguelen animals
BBAL[BBAL$id == "12119005_14257_2001", "stage"] <- "incubation"

#"12139117_20876_1999"
BBAL[BBAL$id == "12139117_20876_1999", ]
newrow <- data.frame(id = "12139117_20876_1999",
                     stage = "incubation",
                     start = "1999-12-20 05:54:25",
                     end = "2000-01-03 06:54:25")
BBAL[BBAL$id == "12139117_20876_1999" & BBAL$stage == "chick-rearing", "start"] <- "2000-01-03 06:54:25"
BBAL <- rbind(BBAL, newrow)
rm(newrow)

BBAL[BBAL$id == "12143525_1203_2005" & BBAL$stage == "post-breeding", "end"] <- "2006-09-14 22:12:39"
BBAL[BBAL$id == "12143525_1203_2005" & BBAL$stage == "arrival", "start"] <- "2006-09-14 22:12:39"

#unknown6
BBAL[BBAL$id == "unknown6" & BBAL$stage == "chick-rearing", "end"] <- "2004-02-16 10:23:44"
newrow <- data.frame(id = "unknown6",
                     stage = "post-breeding",
                     start = "2004-02-16 10:23:44",
                     end = "2004-02-23 16:23:44")
BBAL <- rbind(BBAL, newrow)
rm(newrow)

BBAL <- BBAL[!(BBAL$id %in% fails), ]

saveRDS(BBAL, "./data_split/stage_dates/BBAL_stages.RDS")


#--------------------
#DMSA

DMSA <- readRDS("./data_split/stage_dates/DMSA_stages.RDS")

#'105299'
DMSA[DMSA$id == "105299" & DMSA$stage == "late chick-rearing", "end"] <- "2011-04-10 15:36:47"
newrow <- data.frame(id = "105299",
                     stage = "post-breeding",
                     start = "2011-04-10 15:36:47",
                     end = "2011-05-03 11:36:47")
DMSA <- rbind(DMSA, newrow)
rm(newrow)

#'651_7340'
DMSA[DMSA$id == "651_7340" & DMSA$stage == "late chick-rearing", "end"] <- "2008-05-21 11:01:41"
newrow <- data.frame(id = "651_7340",
                     stage = "post-breeding",
                     start = "2008-05-21 11:01:41",
                     end = "2008-06-18 21:01:41")
DMSA <- rbind(DMSA, newrow)
rm(newrow)

#'651_7341'
DMSA[DMSA$id == "651_7341" & DMSA$stage == "late chick-rearing", "end"] <- "2008-04-16 13:32:55"
newrow <- data.frame(id = "651_7341",
                     stage = "post-breeding",
                     start = "2008-04-16 13:32:55",
                     end = "2008-06-04 23:32:55")
DMSA <- rbind(DMSA, newrow)
rm(newrow)


#'651_7341'
DMSA[DMSA$id == "651_7343" & DMSA$stage == "late chick-rearing", "end"] <- "2009-05-25 08:06:36"
newrow <- data.frame(id = "651_7343",
                     stage = "post-breeding",
                     start = "2009-05-25 08:06:36",
                     end = "2009-07-11 12:06:36")
DMSA <- rbind(DMSA, newrow)
rm(newrow)

#651_7344
DMSA[DMSA$id == "651_7344" & DMSA$stage == "late chick-rearing", "end"] <- "2009-05-15 18:15:36"
DMSA[DMSA$id == "651_7344" & DMSA$stage == "post-breeding", "start"] <- "2009-05-15 18:15:36"


#'651_7341'
DMSA[DMSA$id == "651_7346" & DMSA$stage == "post-breeding", "start"] <- "2009-06-08 07:22:13"
newrow <- data.frame(id = "651_7346",
                     stage = "late chick-rearing",
                     start = "2009-04-24 07:22:13",
                     end = "2009-06-08 07:22:13")
DMSA <- rbind(DMSA, newrow)
rm(newrow)

#'651_7348'
DMSA[DMSA$id == "651_7348" & DMSA$stage == "post-breeding", "start"] <- "2009-05-10 02:16:56"
newrow <- data.frame(id = "651_7348",
                     stage = "late chick-rearing",
                     start = "2009-04-16 14:16:56",
                     end = "2009-05-10 02:16:56")
DMSA <- rbind(DMSA, newrow)
rm(newrow)

#651_7349
DMSA[DMSA$id == "651_7349" & DMSA$stage == "late chick-rearing", "end"] <- "2009-05-22 08:14:58"
DMSA[DMSA$id == "651_7349" & DMSA$stage == "post-breeding", "start"] <- "2009-05-22 08:14:58"

#'98023'
DMSA[DMSA$id == "98023" & DMSA$stage == "late chick-rearing", "end"] <- "2011-04-28 10:48:46"
newrow <- data.frame(id = "98023",
                     stage = "post-breeding",
                     start = "2011-04-28 10:48:46",
                     end = "2011-06-30 08:48:46")
DMSA <- rbind(DMSA, newrow)
rm(newrow)

DMSA <- DMSA[!(DMSA$id %in% fails), ]

saveRDS(DMSA, "./data_split/stage_dates/DMSA_stages.RDS")


#--------------------
#EMPE
EMPE <- readRDS("./data_split/stage_dates/EMPE_stages.RDS")

EMPE[EMPE$id == "CI07321_1993", "stage"] <- "post-breeding"
EMPE[EMPE$id == "CR06886_1993", "stage"] <- "post-breeding"
EMPE[EMPE$id == "CR06889_1993", "stage"] <- "post-breeding"
EMPE[EMPE$id == "CR06890_1993", "stage"] <- "post-breeding"
EMPE[EMPE$id == "DAV2011_emp_f_x_3175", "stage"] <- "post-breeding"
EMPE[EMPE$id == "DDU2005_emp_a_f_19", "stage"] <- "post-breeding"
EMPE[EMPE$id == "DDU2005_emp_a_m_14", "stage"] <- "post-breeding"
EMPE[EMPE$id == "DDU2005_emp_a_m_15", "stage"] <- "post-breeding"
EMPE[EMPE$id == "DDU2005_emp_a_x_06", "stage"] <- "post-breeding"
EMPE[EMPE$id == "DDU2005_emp_a_x_16", "stage"] <- "post-breeding"
EMPE[EMPE$id == "DAV2010_emp_f_x_3010", "stage"] <- "post-breeding"
EMPE[EMPE$id == "DAV2011_emp_f_x_3178", "stage"] <- "post-breeding"
EMPE[EMPE$id == "MAW2000_emp_a_x_762", "stage"] <- "post-breeding"
EMPE[EMPE$id == "MAW2000_emp_a_x_764", "stage"] <- "post-breeding"
EMPE[EMPE$id == "MAW2000_emp_a_x_765", "stage"] <- "post-breeding"
EMPE[EMPE$id == "MAW2000_emp_a_x_766", "stage"] <- "post-breeding"
EMPE[EMPE$id == "MAW2000_emp_a_x_767", "stage"] <- "post-breeding"
EMPE[EMPE$id == "MAW2011_emp_f_x_3166", "stage"] <- "post-breeding"
EMPE[EMPE$id == "MAW2011_emp_f_x_3166", "stage"] <- "post-breeding"

# CW06885_1993
# last trip should post-breeding
EMPE[EMPE$id == "CW06885_1993", ]
newrow <- data.frame(id = "CW06885_1993",
                     stage = "post-breeding",
                     start = "1993-12-19 11:23:00",
                     end = "1994-01-04 09:23:00")
EMPE[EMPE$id == "CW06885_1993" & EMPE$stage == "chick-rearing", "end"] <- "1993-12-19 11:23:00"
EMPE <- rbind(EMPE, newrow)

# DDU2005_emp_a_f_07
# last trip should post-breeding
EMPE[EMPE$id == "DDU2005_emp_a_f_07", ]
newrow <- data.frame(id = "DDU2005_emp_a_f_07",
                     stage = "post-breeding",
                     start = "2005-11-24 12:24:06",
                     end = "2006-01-15 00:24:06")
EMPE[EMPE$id == "DDU2005_emp_a_f_07" & EMPE$stage == "chick-rearing", "end"] <- "2005-11-24 12:24:06"
EMPE <- rbind(EMPE, newrow)

# DDU2005_emp_a_f_18b
# last trip should post-breeding
EMPE[EMPE$id == "DDU2005_emp_a_f_18b", ]
newrow <- data.frame(id = "DDU2005_emp_a_f_18b",
                     stage = "post-breeding",
                     start = "2005-12-07 05:16:47",
                     end = "2006-01-15 09:16:47")
EMPE[EMPE$id == "DDU2005_emp_a_f_18b" & EMPE$stage == "chick-rearing", "end"] <- "2005-12-07 05:16:47"
EMPE <- rbind(EMPE, newrow)

# DDU2005_emp_a_x_03
# last trip should post-breeding
EMPE[EMPE$id == "DDU2005_emp_a_x_03", ]
newrow <- data.frame(id = "DDU2005_emp_a_x_03",
                     stage = "post-breeding",
                     start = "2005-12-01 20:25:17",
                     end = "2005-12-30 10:25:17")
EMPE[EMPE$id == "DDU2005_emp_a_x_03" & EMPE$stage == "chick-rearing", "end"] <- "2005-12-01 20:25:17"
EMPE <- rbind(EMPE, newrow)

# CW06886_1992
# last trip should post-breeding
EMPE[EMPE$id == "CW06886_1992", ]
newrow <- data.frame(id = "CW06886_1992",
                     stage = "post-breeding",
                     start = "1992-11-04 10:55:00",
                     end = "1992-11-16 16:55:00")
EMPE[EMPE$id == "CW06886_1992" & EMPE$stage == "chick-rearing", "end"] <- "1992-11-04 10:55:00"
EMPE <- rbind(EMPE, newrow)

# CW06886_1992
# last trip should post-breeding
EMPE[EMPE$id == "MAW2011_emp_f_x_3171", ]
newrow <- data.frame(id = "MAW2011_emp_f_x_3171",
                     stage = "post-breeding",
                     start = "2011-12-24 09:09:47",
                     end = "2012-01-12 05:09:47")
EMPE[EMPE$id == "MAW2011_emp_f_x_3171" & EMPE$stage == "chick-rearing", "end"] <- "2011-12-24 09:09:47"
EMPE <- rbind(EMPE, newrow)

EMPE <- EMPE[!(EMPE$id %in% fails), ]

saveRDS(EMPE, "./data_split/stage_dates/EMPE_stages.RDS")


#--------------------
#GHAL
GHAL <- readRDS("./data_split/stage_dates/GHAL_stages.RDS")

GHAL[GHAL$id == "12139225_1191_2005" & GHAL$stage == "late chick-rearing", "stage"] <- "post-breeding"

GHAL[GHAL$id == "A0404_20874_2000" & GHAL$stage == "incubation", "end"] <- "2000-12-23 10:02:34"
GHAL[GHAL$id == "A0404_20874_2000" & GHAL$stage == "early chick-rearing", "start"] <- "2000-12-23 10:02:34"

GHAL <- GHAL[!(GHAL$id %in% fails), ]

saveRDS(GHAL, "./data_split/stage_dates/GHAL_stages.RDS")


#--------------------
#KIPE

KIPE <- readRDS("./data_split/stage_dates/KIPE_stages.RDS")

KIPE[KIPE$id == "57332" & KIPE$stage == "Pre-breeding", "stage"] <- "incubation"
KIPE[KIPE$id == "Saanenland" & KIPE$stage == "incubation", "stage"] <- "chick-rearing"
KIPE[KIPE$id == "Tankini" & KIPE$stage == "incubation", "stage"] <- "chick-rearing"
KIPE[KIPE$id == "Ueli" & KIPE$stage == "incubation", "stage"] <- "chick-rearing"

KIPE <- KIPE[!(KIPE$id %in% fails), ]

saveRDS(KIPE, "./data_split/stage_dates/KIPE_stages.RDS")


#--------------------
#LMSA

LMSA <- readRDS("./data_split/stage_dates/LMSA_stages.RDS")

LMSA[LMSA$id == "650_7336", "stage"] <- "post-breeding"
LMSA[LMSA$id == "98026", "stage"] <- "post-breeding"
LMSA[LMSA$id == "105298", "stage"] <- "post-breeding"

# 650_7335
# last trip should be post-breeding
LMSA[LMSA$id == "650_7335", ]
newrow <- data.frame(id = "650_7335",
                     stage = "post-breeding",
                     start = "2009-01-09 07:39:24",
                     end = "2009-04-08 13:39:24")
LMSA[LMSA$id == "650_7335" & LMSA$stage == "late chick-rearing", "end"] <- "2009-01-09 07:39:24"
LMSA <- rbind(LMSA, newrow)

#Karine eds
LMSA[LMSA$id == "12143650_14257_2002", "stage"] <- "late chick-rearing"
LMSA[LMSA$id == "12145258_37467_2002", "stage"] <- "late chick-rearing"
LMSA[LMSA$id == "650_7333", "stage"] <- "post-breeding"

# 650_7335
# last trip should be post-breeding
newrow <- data.frame(id = "650_7335",
                     stage = "early chick-rearing",
                     start = "2008-12-24 10:11:29",
                     end = "2009-01-02 01:39:24")
LMSA[LMSA$id == "650_7335" & LMSA$stage == "late chick-rearing", "start"] <- "2009-01-02 01:39:24"
LMSA <- rbind(LMSA, newrow)

#'650_7337'
LMSA[LMSA$id == "650_7337" & LMSA$stage == "late chick-rearing", "end"] <- "2009-05-25 06:39:48"
LMSA[LMSA$id == "650_7337" & LMSA$stage == "post-breeding", "start"] <- "2009-05-25 06:39:48"

#'650_7338'
newrow <- data.frame(id = "650_7338",
                     stage = "early chick-rearing",
                     start = "2009-04-24 07:22:13",
                     end = "2009-05-25 09:22:13")
LMSA[LMSA$id == "650_7338" & LMSA$stage == "post-breeding", "start"] <- "2009-05-25 09:22:13"
LMSA <- rbind(LMSA, newrow)

#'A1084_20877_2002'
newrow <- data.frame(id = "A1084_20877_2002",
                     stage = "early chick-rearing",
                     start = "2003-01-05 05:31:02",
                     end = "2003-01-15 15:31:02")
LMSA[LMSA$id == "A1084_20877_2002" & LMSA$stage == "incubation", "end"] <- "2003-01-05 05:31:02"
LMSA <- rbind(LMSA, newrow)

#'LMSA-unknown15-44576-RAATD'
newrow <- data.frame(id = "LMSA-unknown15-44576-RAATD",
                     stage = "post-breeding",
                     start = "2004-01-07 03:28:04",
                     end = "2004-02-15 03:28:04")
LMSA[LMSA$id == "LMSA-unknown15-44576-RAATD" & LMSA$stage == "early chick-rearing", "end"] <- "2004-01-07 03:28:04"
LMSA <- rbind(LMSA, newrow)

LMSA <- LMSA[!(LMSA$id %in% fails), ]

saveRDS(LMSA, "./data_split/stage_dates/LMSA_stages.RDS")


#--------------------
#MAPE

MAPE <- readRDS("./data_split/stage_dates/MAPE_stages.RDS")

# Assign stages for Marion data based on names

for (i in 1:nrow(MAPE)) {
  if (substr(x = MAPE[i, "id"], start = 1, stop = 2) == "CR") {
    MAPE[i, "stage"] <- "early chick-rearing"
  }
  if (substr(x = MAPE[i, "id"], start = 1, stop = 3) == "Cre") {
    MAPE[i, "stage"] <- "late chick-rearing"
  }
}

# Stage merges
MAPE <- stageMerger(dat = MAPE, id = "MAPE-dtsetBirdLife751-22-RAATD", becomes = "early chick-rearing")
MAPE <- stageMerger(dat = MAPE, id = "MAPE-dtsetBirdLife751-A46-RAATD", becomes = "early chick-rearing")
MAPE <- stageMerger(dat = MAPE, id = "MAPE-dtsetBirdLife751-A5-RAATD", becomes = "early chick-rearing")
MAPE <- stageMerger(dat = MAPE, id = "M39_28947", becomes = "early chick-rearing")
MAPE <- stageMerger(dat = MAPE, id = "MAPE-dtsetBirdLife751-A8-RAATD", becomes = "early chick-rearing")
MAPE <- stageMerger(dat = MAPE, id = "MAPE-dtsetBirdLife751-H49-RAATD", becomes = "early chick-rearing")
MAPE <- stageMerger(dat = MAPE, id = "MAPE-dtsetBirdLife751-v4-RAATD", becomes = "early chick-rearing")
MAPE <- stageMerger(dat = MAPE, id = "MAPE-dtsetBirdLife751-X15-RAATD", becomes = "early chick-rearing")
MAPE <- stageMerger(dat = MAPE, id = "MAPE-dtsetBirdLife751-X17-RAATD", becomes = "early chick-rearing")

# Specific stage changes
MAPE[MAPE$id == "M30_44554" & MAPE$stage == "incubation", "stage"] <- "early chick-rearing"
MAPE[MAPE$id == "M30_44554" & MAPE$stage == "early chick-rearing", "stage"] <- "late chick-rearing"

# Stage reassignment
MAPE[MAPE$id == "M36_44551", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "M100_44551", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "M11_44557", "stage"] <- "early chick-rearing"
MAPE[MAPE$id == "M49_44557", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "M58_29817", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "M7_44556", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "M70_26139", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "M87_29987", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "M88_44552", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "M96_44548", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-14-RAATD", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-3.2-RAATD", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-4.4-RAATD", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-A35-RAATD", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-H36-RAATD", "stage"] <- "late chick-rearing"
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-H63-RAATD", "stage"] <- "early chick-rearing"
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-X16-RAATD_1999_2", "stage"] <- "early chick-rearing"
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-X32-RAATD", "stage"] <- "post-breeding"

# Stage splits

# MAPE-dtsetBirdLife743-Y33-RAATD
# last trip should be late chick-rearing
MAPE[MAPE$id == "MAPE-dtsetBirdLife743-Y33-RAATD", ]
newrow <- data.frame(id = "MAPE-dtsetBirdLife743-Y33-RAATD",
                     stage = "late chick-rearing",
                     start = "2003-01-25 00:11:00",
                     end = "2003-02-02 08:57:37")
MAPE[MAPE$id == "MAPE-dtsetBirdLife743-Y33-RAATD" & MAPE$stage == "early chick-rearing", "end"] <- "2003-01-25 00:11:00"
MAPE <- rbind(MAPE, newrow)


# MAPE-dtsetBirdLife751-A43-RAATD
# last trip should be late chick-rearing
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-A43-RAATD", ]
newrow <- data.frame(id = "MAPE-dtsetBirdLife751-A43-RAATD",
                     stage = "late chick-rearing",
                     start = "2005-01-29 23:03:00",
                     end = "2005-02-03 15:03:00")
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-A43-RAATD" & MAPE$stage == "early chick-rearing", "end"] <- "2005-01-29 23:03:00"
MAPE <- rbind(MAPE, newrow)


# MAPE-dtsetBirdLife751-13-RAATD
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-13-RAATD", ]
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-13-RAATD" & MAPE$stage == "early chick-rearing", "end"] <- "2003-01-27 00:32:00"
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-13-RAATD" & MAPE$stage == "late chick-rearing", "start"] <- "2003-01-27 00:32:00"


# MAPE-dtsetBirdLife751-A55-RAATD
# First trips should be early chick-rearing
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-A55-RAATD" & MAPE$stage == "late chick-rearing", "start"] <- "2005-02-06 07:21:00"
newrow <- data.frame(id = "MAPE-dtsetBirdLife751-A55-RAATD",
                     stage = "early chick-rearing",
                     start = "2005-02-03 07:21:00",
                     end = "2005-02-06 07:21:00")
MAPE <- rbind(MAPE, newrow)


# MAPE-dtsetBirdLife751-A57-RAATD
# First trips should be early chick-rearing
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-A57-RAATD" & MAPE$stage == "late chick-rearing", "start"] <- "2005-02-05 07:36:00"
newrow <- data.frame(id = "MAPE-dtsetBirdLife751-A57-RAATD",
                     stage = "early chick-rearing",
                     start = "2005-02-02 21:36:00",
                     end = "2005-02-05 07:36:00")
MAPE <- rbind(MAPE, newrow)


# MAPE-dtsetBirdLife751-A61-RAATD
# First trips should be early chick-rearing
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-A61-RAATD" & MAPE$stage == "late chick-rearing", "start"] <- "2005-02-06 04:39:00"
newrow <- data.frame(id = "MAPE-dtsetBirdLife751-A61-RAATD",
                     stage = "early chick-rearing",
                     start = "2005-02-04 08:39:00",
                     end = "2005-02-06 04:39:00")
MAPE <- rbind(MAPE, newrow)


# MAPE-dtsetBirdLife751-A62-RAATD
# First trips should be early chick-rearing
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-A62-RAATD" & MAPE$stage == "late chick-rearing", "start"] <- "2005-02-05 09:00:00"
newrow <- data.frame(id = "MAPE-dtsetBirdLife751-A62-RAATD",
                     stage = "early chick-rearing",
                     start = "2005-02-04 07:00:00",
                     end = "2005-02-05 09:00:00")
MAPE <- rbind(MAPE, newrow)


# MAPE-dtsetBirdLife751-H35-RAATD
# First trips should be early chick-rearing
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-H35-RAATD" & MAPE$stage == "late chick-rearing", "start"] <- "2004-01-24 06:30:00"
newrow <- data.frame(id = "MAPE-dtsetBirdLife751-H35-RAATD",
                     stage = "early chick-rearing",
                     start = "2004-01-22 06:30:00",
                     end = "2004-01-24 06:30:00")
MAPE <- rbind(MAPE, newrow)


# MAPE-dtsetBirdLife751-H38-RAATD
# First trips should be early chick-rearing
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-H38-RAATD" & MAPE$stage == "late chick-rearing", "start"] <- "2004-01-25 18:32:00"
newrow <- data.frame(id = "MAPE-dtsetBirdLife751-H38-RAATD",
                     stage = "early chick-rearing",
                     start = "2004-01-23 06:32:00",
                     end = "2004-01-25 18:32:00")
MAPE <- rbind(MAPE, newrow)


# MAPE-dtsetBirdLife751-H39-RAATD
# First trips should be early chick-rearing
MAPE[MAPE$id == "MAPE-dtsetBirdLife751-H39-RAATD" & MAPE$stage == "late chick-rearing", "start"] <- "2004-01-26 01:58:00"
newrow <- data.frame(id = "MAPE-dtsetBirdLife751-H39-RAATD",
                     stage = "early chick-rearing",
                     start = "2004-01-23 07:58:00",
                     end = "2004-01-26 01:58:00")
MAPE <- rbind(MAPE, newrow)

MAPE <- MAPE[!(MAPE$id %in% fails), ]

saveRDS(MAPE, "./data_split/stage_dates/MAPE_stages.RDS")


#--------------------
#ROPE

ROPEdates <- readRDS("./data_split/stage_dates/ROPE_stages.RDS")

#'r5' should be incubation
ROPEdates[ROPEdates$id == "r5", "stage"] <- "incubation"
ROPEdates[ROPEdates$id == "r6", "stage"] <- "incubation"
ROPEdates[ROPEdates$id == "r8", "stage"] <- "incubation"

ROPEdates <- ROPEdates[!(ROPEdates$id %in% fails), ]

saveRDS(ROPEdates, "./data_split/stage_dates/ROPE_stages.RDS")


#--------------------
#SOES

SOESdates <- readRDS("./data_split/stage_dates/SOES_stages.RDS")

#'CAR2013_sel_a_m_05' should all be post-breeding
SOESdates[SOESdates$id == "CAR2013_sel_a_m_05", ]
SOESdates[SOESdates$id == "CAR2013_sel_a_m_05" & SOESdates$stage == "post-breeding", "end"] <- "2014-03-16 00:25:00"
SOESdates <- SOESdates[!(SOESdates$id == "CAR2013_sel_a_m_05" & SOESdates$stage == "post-moult"), ]

#'CAR2013_sel_a_m_14' should all be post-breeding
SOESdates[SOESdates$id == "CAR2013_sel_a_m_14", ]
SOESdates[SOESdates$id == "CAR2013_sel_a_m_14" & SOESdates$stage == "post-breeding", "end"] <- "2014-04-12 03:43:00"
SOESdates <- SOESdates[!(SOESdates$id == "CAR2013_sel_a_m_14" & SOESdates$stage == "post-moult"), ]

#'CAR2013_sel_a_m_17' should all be post-breeding
SOESdates[SOESdates$id == "CAR2013_sel_a_m_17", ]
SOESdates[SOESdates$id == "CAR2013_sel_a_m_17" & SOESdates$stage == "post-breeding", "end"] <- "2014-04-02 03:39:00"
SOESdates <- SOESdates[!(SOESdates$id == "CAR2013_sel_a_m_17" & SOESdates$stage == "post-moult"), ]

#'CAR2013_sel_a_m_18' should all be post-breeding
SOESdates[SOESdates$id == "CAR2013_sel_a_m_18", ]
SOESdates[SOESdates$id == "CAR2013_sel_a_m_18", "stage"] <- "post-breeding"

#'CAR2013_sel_a_m_24' should all be post-breeding
SOESdates[SOESdates$id == "CAR2013_sel_a_m_24", ]
SOESdates[SOESdates$id == "CAR2013_sel_a_m_24", "stage"] <- "post-breeding"

#'ct46-69047-09' should all be post-moult
SOESdates[SOESdates$id == "ct46-69047-09", ]
SOESdates[SOESdates$id == "ct46-69047-09" & SOESdates$stage == "post-moult", "start"] <- "2009-01-25 19:38:17"
SOESdates <- SOESdates[!(SOESdates$id == "ct46-69047-09" & SOESdates$stage == "post-breeding"), ]

#'ft02-G-09' should all be post-moult
SOESdates[SOESdates$id == "ft02-G-09", ]
SOESdates[SOESdates$id == "ft02-G-09" & SOESdates$stage == "post-moult", "start"] <- "2008-12-24 10:54:55"
SOESdates <- SOESdates[!(SOESdates$id == "ft02-G-09" & SOESdates$stage == "post-breeding"), ]

SOESdates <- SOESdates[!(SOESdates$id %in% fails), ]

saveRDS(SOESdates, "./data_split/stage_dates/SOES_stages.RDS")


#--------------------
#WAAL
WAAL <- readRDS("./data_split/stage_dates/WAAL_stages.RDS")

#'14049642_37631_2006'
WAAL[WAAL$id == "14049642_37631_2006" & WAAL$stage == "chick-rearing", "end"] <- "2006-03-27 04:00:42"
newrow <- data.frame(id = "14049642_37631_2006",
                     stage = "post-breeding",
                     start = "2006-03-27 04:00:42",
                     end = "2006-07-07 00:00:42")
WAAL <- rbind(WAAL, newrow)

#'14052375_37629_2006'
WAAL[WAAL$id == "14052375_37629_2006" & WAAL$stage == "chick-rearing", "stage"] <- "post-breeding"

#'14052638 (154)_55167_2004'
WAAL[WAAL$id == "14052638 (154)_55167_2004" & WAAL$stage == "incubation", "stage"] <- "post-breeding"

#'unknown'
WAAL[WAAL$id == "unknown" & WAAL$stage == "chick-rearing", "stage"] <- "post-breeding"

#'waal_135F_5202390'
WAAL[WAAL$id == "waal_135F_5202390" & WAAL$stage == "incubation", "stage"] <- "post-breeding"

#'waal_501F_5098530'
WAAL[WAAL$id == "waal_501F_5098530" & WAAL$stage == "incubation", "stage"] <- "post-breeding"

#'WAAL-Crozet-BS11182 -9016-2008_WINTER1_analyse.txt'
WAAL[WAAL$id == "WAAL-Crozet-BS11182 -9016-2008_WINTER1_analyse.txt", ]
WAAL[WAAL$id == "WAAL-Crozet-BS11182 -9016-2008_WINTER1_analyse.txt" & WAAL$stage == "chick-rearing", "end"] <- "2010-04-13 22:13:06" 
WAAL <- WAAL[!(WAAL$id == "WAAL-Crozet-BS11182 -9016-2008_WINTER1_analyse.txt" & WAAL$stage == "post-breeding"), ]
newrow <- data.frame(id = "WAAL-Crozet-BS11182 -9016-2008_WINTER1_analyse.txt",
                     stage = "post-breeding",
                     start = "2010-04-13 22:13:06",
                     end = "2010-12-14 21:11:01")
WAAL <- rbind(WAAL, newrow)
rm(newrow)

#'WAAL-Crozet-BS11283-9070-2008_WINTER1_analyse.txt'
WAAL <- WAAL[!(WAAL$id == "WAAL-Crozet-BS11283-9070-2008_WINTER1_analyse.txt"), ]
start <- "2010-12-24 11:19:00"
end <- "2011-11-29 23:19:00"
newrow <- data.frame(id = "WAAL-Crozet-BS11283-9070-2008_WINTER1_analyse.txt",
                     stage = "post-breeding",
                     start = start,
                     end = end)
WAAL <- rbind(WAAL, newrow)

#Set to post-breeding
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS11979-9012-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS12318-5566-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS12505-5319-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS12636-9005-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS18847-5317-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS19239-1087-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS19280-9008-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS19469-5316-2007_WINTER2_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS19719-1116-2005_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS19748-9014-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS19913-9002-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20075-9071-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20090-10302-2009_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20090-5555-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20207-1118-2005_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20227-5557-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20233-5322-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20479-22663-2010_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20520-10367-2009_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20561-1120-2005_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20575-5565-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20693-9003-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20752-22670-2010_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS20776-10309-2009_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS21135-5320-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS21258-22664-2010_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS21302-9004-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS21365-10304-2009_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS21387-22669-2010_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS21458-24913-2011_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS21794-5561-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS22900-5549-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS22903-22681-2010_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS23203-5547-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS23207-10388-2009_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS23605-5559-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS25167-22677-2010_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS25167-22677-2010_WINTER2_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS25209-5562-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS25236-9011-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS25238-10308-2009_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS25844-5321-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS27614-22658-2010_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS27709-532-2006_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS27805 -9007-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS27805 -9007-2008_WINTER3_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS27813-10307-2009_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS5115-9072-2008_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS5728-5324-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS6229-5311-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS6975-5323-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS8581-10363-2009_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS8831-5554-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS9067-5569-2007_WINTER2_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS9466-22660-2010_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS9466-22660-2010_WINTER2_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Crozet-BS9508-5312-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS17369-24939-2011_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS17371-5334-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS20856-1099-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS20941-24936-2011_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS23408-5332-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS24207-24940-2011_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS24288-5346-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS26334-24930-2011_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS5526-5335-2007_WINTER1_analyse.txt", "post-breeding")
WAAL <- stageMerger(WAAL, "WAAL-Kerguelen-BS5538-5328-2007_WINTER1_analyse.txt", "post-breeding")

WAAL[WAAL$id == "WAAL-Crozet-BS12502-541-2006_WINTER1_analyse.txt" & WAAL$stage == "chick-rearing", "stage"] <- "post-breeding"
WAAL[WAAL$id == "WAAL-Crozet-BS19818-24912-2011_WINTER1_analyse.txt" & WAAL$stage == "chick-rearing", "stage"] <- "post-breeding"

#other fixes

#'WAAL-Crozet-BS19469-5316-2007_WINTER1_analyse.txt'
WAAL <- WAAL[!(WAAL$id == "WAAL-Crozet-BS19469-5316-2007_WINTER1_analyse.txt" & WAAL$stage == "chick-rearing"), ]
WAAL[WAAL$id == "WAAL-Crozet-BS19469-5316-2007_WINTER1_analyse.txt" & WAAL$stage == "post-breeding", "start"] <- "2009-03-06 14:28:46"

#'WAAL-Crozet-BS20054-542-2006_WINTER1_analyse.txt'
WAAL[WAAL$id == "WAAL-Crozet-BS20054-542-2006_WINTER1_analyse.txt" & WAAL$stage == "post-breeding", "end"] <- "2007-12-01 15:46:00"
WAAL[WAAL$id == "WAAL-Crozet-BS20054-542-2006_WINTER1_analyse.txt" & WAAL$stage == "incubation", "start"] <- "2007-12-01 15:46:00"

WAAL[WAAL$id == "WAAL-Crozet-BS20054-542-2006_WINTER2_analyse.txt" & WAAL$stage == "chick-rearing", "stage"] <- "incubation"
WAAL[WAAL$id == "WAAL-Crozet-BS20259-22676-2010_WINTER1_analyse.txt" & WAAL$stage == "chick-rearing", "stage"] <- "post-breeding"

#'WAAL-Crozet-BS20589-22667-2010_WINTER1_analyse.txt'
WAAL <- WAAL[!(WAAL$id == "WAAL-Crozet-BS20589-22667-2010_WINTER1_analyse.txt" & WAAL$stage == "chick-rearing"), ]
WAAL[WAAL$id == "WAAL-Crozet-BS20589-22667-2010_WINTER1_analyse.txt" & WAAL$stage == "post-breeding", "end"] <- "2011-12-08 13:17:59"

#"WAAL-Crozet-BS21107-5325-2007_WINTER1_analyse.txt"
WAAL <- WAAL[!(WAAL$id == "WAAL-Crozet-BS21107-5325-2007_WINTER1_analyse.txt" & WAAL$stage == "chick-rearing"), ]
WAAL <- WAAL[!(WAAL$id == "WAAL-Crozet-BS21107-5325-2007_WINTER1_analyse.txt" & WAAL$stage == "post-breeding"), ]
newrow <- data.frame(id = "WAAL-Crozet-BS21107-5325-2007_WINTER1_analyse.txt",
                      stage = "post-breeding",
                      start = "2007-12-01 03:29:01",
                      end = "2008-12-06 18:20:09")
WAAL <- rbind(WAAL, newrow)

#"WAAL-Crozet-BS21228-5574-2007_WINTER1_analyse.txt"
WAAL[WAAL$id == "WAAL-Crozet-BS21228-5574-2007_WINTER1_analyse.txt", ]

WAAL <- WAAL[WAAL$id != "WAAL-Crozet-BS21228-5574-2007_WINTER1_analyse.txt", ]
newrow1 <- data.frame(id = "WAAL-Crozet-BS21228-5574-2007_WINTER1_analyse.txt",
                                stage = "post-breeding",
                                start = "2008-11-30 09:50:43",
                                end = "2009-12-14 23:25:20")
newrow2 <- data.frame(id = "WAAL-Crozet-BS21228-5574-2007_WINTER1_analyse.txt",
                      stage = "incubation",
                      start = "2009-12-14 23:25:20",
                      end = "2009-12-26 21:50:43")
WAAL <- rbind(WAAL, newrow1, newrow2)
rm(newrow1, newrow2)

#"WAAL-Crozet-BS22630-10368-2009_WINTER1_analyse.txt"
WAAL <- WAAL[WAAL$id != "WAAL-Crozet-BS22630-10368-2009_WINTER1_analyse.txt", ]
newrow1 <- data.frame(id = "WAAL-Crozet-BS22630-10368-2009_WINTER1_analyse.txt",
                      stage = "post-breeding",
                      start = "2010-04-10 21:18:00",
                      end = "2010-11-13 08:02:53")
newrow2 <- data.frame(id = "WAAL-Crozet-BS22630-10368-2009_WINTER1_analyse.txt",
                      stage = "incubation",
                      start = "2010-11-13 08:02:53",
                      end = "2010-12-07 08:18:00")
WAAL <- rbind(WAAL, newrow1, newrow2)
rm(newrow1, newrow2)

#"WAAL-Crozet-BS22864-540-2006_WINTER1_analyse.txt"
WAAL[WAAL$id == "WAAL-Crozet-BS22864-540-2006_WINTER1_analyse.txt" & WAAL$stage == "incubation", "end"] <- "2008-01-02 20:28:00"

#"WAAL-Crozet-BS23091-10370-2009_WINTER1_analyse.txt"
WAAL <- WAAL[WAAL$id != "WAAL-Crozet-BS23091-10370-2009_WINTER1_analyse.txt", ]
newrow1 <- data.frame(id = "WAAL-Crozet-BS23091-10370-2009_WINTER1_analyse.txt",
                      stage = "post-breeding",
                      start = "2010-11-30 06:52:42",
                      end = "2011-12-13 08:25:35")
newrow2 <- data.frame(id = "WAAL-Crozet-BS23091-10370-2009_WINTER1_analyse.txt",
                      stage = "incubation",
                      start = "2011-12-13 08:25:35",
                      end = "2011-12-30 18:52:42")
WAAL <- rbind(WAAL, newrow1, newrow2)
rm(newrow1, newrow2)

#"WAAL-Crozet-BS25422-22665-2010_WINTER1_analyse.txt"
WAAL <- WAAL[WAAL$id != "WAAL-Crozet-BS25422-22665-2010_WINTER1_analyse.txt", ]
newrow1 <- data.frame(id = "WAAL-Crozet-BS25422-22665-2010_WINTER1_analyse.txt",
                      stage = "post-breeding",
                      start = "2011-12-14 08:18:52",
                      end = "2012-12-23 13:18:36")
newrow2 <- data.frame(id = "WAAL-Crozet-BS25422-22665-2010_WINTER1_analyse.txt",
                      stage = "incubation",
                      start = "2012-12-23 13:18:36",
                      end = "2012-12-31 20:18:52")
WAAL <- rbind(WAAL, newrow1, newrow2)
rm(newrow1, newrow2)

#"WAAL-Crozet-BS27642-22675-2010_WINTER1_analyse.txt"
WAAL <- WAAL[WAAL$id != "WAAL-Crozet-BS27642-22675-2010_WINTER1_analyse.txt", ]
newrow1 <- data.frame(id = "WAAL-Crozet-BS27642-22675-2010_WINTER1_analyse.txt",
                      stage = "post-breeding",
                      start = "2011-11-30 08:04:37",
                      end = "2012-12-12 17:13:03")
newrow2 <- data.frame(id = "WAAL-Crozet-BS27642-22675-2010_WINTER1_analyse.txt",
                      stage = "incubation",
                      start = "2012-12-12 17:13:03",
                      end = "2012-12-20 08:04:37")
WAAL <- rbind(WAAL, newrow1, newrow2)
rm(newrow1, newrow2)

#"WAAL-Kerguelen-BS17556-5342-2007_WINTER1_analyse.txt"
WAAL <- WAAL[!(WAAL$id == "WAAL-Kerguelen-BS17556-5342-2007_WINTER1_analyse.txt" & WAAL$stage == "chick-rearing"), ]
WAAL <- WAAL[!(WAAL$id == "WAAL-Kerguelen-BS17556-5342-2007_WINTER1_analyse.txt" & WAAL$stage == "incubation"), ]

WAAL[WAAL$id == "14052370 (152)_55166_2004" & WAAL$stage == "incubation", "stage"] <- "post-breeding"

WAAL <- WAAL[!(WAAL$id %in% fails), ]

saveRDS(WAAL, "./data_split/stage_dates/WAAL_stages.RDS")


#--------------------
#WHCP
WHCP <- readRDS("./data_split/stage_dates/WHCP_stages.RDS")

#'wcp_CRs_nest_45_gps240_female.csv'
WHCP[WHCP$id == "wcp_CRs_nest_45_gps240_female.csv" & WHCP$stage == "early chick-rearing", "end"] <- "2013-01-23 04:02:52"
newrow <- data.frame(id = "wcp_CRs_nest_45_gps240_female.csv",
                     stage = "post-breeding",
                     start = "2013-01-23 04:02:52",
                     end = "2013-02-01 12:02:52")
WHCP <- rbind(WHCP, newrow)

#'wcp_CRs_nest_46_gps_37_male.csv'
WHCP[WHCP$id == "wcp_CRs_nest_46_gps_37_male.csv" & WHCP$stage == "incubation", "end"] <- "2013-01-11 04:19:34"
newrow <- data.frame(id = "wcp_CRs_nest_46_gps_37_male.csv",
                     stage = "early chick-rearing",
                     start = "2013-01-11 04:19:34",
                     end = "2013-01-25 20:19:34")
WHCP <- rbind(WHCP, newrow)

#'wcp_CRs_nest_47_gps_243_male.csv'
WHCP[WHCP$id == "wcp_CRs_nest_47_gps_243_male.csv" & WHCP$stage == "early chick-rearing", "start"] <- "2013-01-18 03:10:02"
newrow <- data.frame(id = "wcp_CRs_nest_47_gps_243_male.csv",
                     stage = "incubation",
                     start = "2013-01-05 22:10:02",
                     end = "2013-01-18 03:10:02")
WHCP <- rbind(WHCP, newrow)

WHCP[WHCP$id == "WHCP_BI-HW04146-9295", "stage"] <- "late chick-rearing"
WHCP[WHCP$id == "WHCP_BI-HW04165-9298", "stage"] <- "late chick-rearing"

#Dump S Georgia post-breeding sections at spring equinox
for (i in 1:nrow(WHCP)) {
  if (substr(WHCP$id[i], 1, 21) == "WHCP-dtsetBirdLife439" & WHCP$stage[i] == "post-breeding") {
    yr <- format(WHCP$end[i], "%Y")
    tm <- paste0(yr, "-09-02 00:00:00")
    WHCP[i, "end"] <- as.character(tm)
  }
}

WHCP <- WHCP[!(WHCP$id %in% fails), ]

saveRDS(WHCP, "./data_split/stage_dates/WHCP_stages.RDS")