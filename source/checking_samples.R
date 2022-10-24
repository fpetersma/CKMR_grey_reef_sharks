## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Load packages
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(tidyverse)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Load in the data
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
csv_folder <- "../../Close-kin mark-recapture/Palmyra/csv/"
samples_issue <- read_csv(paste0(csv_folder, 
                                 "palmyra_tissue_samples_with_issues.csv"),
                          trim_ws = TRUE)
samples_good <- read_csv(paste0(csv_folder, 
                                "palmyra_tissue_samples_usable.csv"),
                         trim_ws = TRUE)
info_blacktip <- read_csv(paste0(csv_folder, 
                                 "palmyra_blacktip_reef_tissue_samples.csv"),
                          trim_ws = TRUE)
info_grey <- read_csv(paste0(csv_folder, 
                             "palmyra_greyreef_tissue_samples.csv"),
                      trim_ws = TRUE)


## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Combine things and create summary statistics
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Add issue column to samples_good to match other columns
samples_good$issue <- "nothing"

## Remove issue samples that are also included in good samples
samples_issue <- samples_issue[!(samples_issue$`fin exposed?` == "no" &
                                   samples_issue$issue == "lost_ethanol"), ]
## Combine the samples into one data frame
samples_all <- rbind(samples_good, samples_issue[, c(1,3,2)])
colnames(samples_all) <- c("fin_clip_id", "box", "issue")

## There are 1287 vials in total, instead of the 1288 on the form. Poteatoe-potaatoe.
length(unique(samples_all$fin_clip_id))
samples_all$fin_clip_id[duplicated(samples_all$fin_clip_id)] # these are the duplicates

## Add additional info
samples_all$info_sheet <- 0
samples_all$info_sheet[samples_all$fin_clip_id %in% info_grey$fin_clip] <- 1
samples_all$info_sheet[samples_all$fin_clip_id %in% info_blacktip$fin_clip] <- 2
samples_all$info_sheet[samples_all$fin_clip_id %in% info_blacktip$fin_clip &
                         samples_all$fin_clip_id %in% info_grey$fin_clip] <- 3

## There are 15 fin_clips that occur in both sheets. Why is that?
## These are the 4 duplicates (= 8 samples). The 7 that are left are:
samples_all[samples_all$info_sheet == 3 & samples_all$issue == "duplicate_id", ]

## Add additional info to good samples
samples_good$info_sheet <- 0
samples_good$info_sheet[samples_good$sample_id %in% info_grey$fin_clip] <- 1
samples_good$info_sheet[samples_good$sample_id %in% info_blacktip$fin_clip] <- 2
samples_good$info_sheet[samples_good$sample_id %in% info_blacktip$fin_clip &
                         samples_good$sample_id %in% info_grey$fin_clip] <- 3

sum(samples_good$info_sheet == 0) # 45 good samples have no info
sum(samples_good$info_sheet == 1) # 836 good samples have grey reef info 
sum(samples_good$info_sheet == 2) # 335 good samples have blacktip reef info
sum(samples_good$info_sheet == 3) # 6 good samples occur in grey and blacktip 

## Add additional info to issue samples
samples_issue$info_sheet <- 0
samples_issue$info_sheet[samples_issue$sample_id %in% info_grey$fin_clip] <- 1
samples_issue$info_sheet[samples_issue$sample_id %in% info_blacktip$fin_clip] <- 2
samples_issue$info_sheet[samples_issue$sample_id %in% info_blacktip$fin_clip &
                          samples_issue$sample_id %in% info_grey$fin_clip] <- 3

sum(samples_issue$info_sheet == 0) # 3 issue samples have no info
sum(samples_issue$info_sheet == 1) # 48 issue samples have grey reef info 
sum(samples_issue$info_sheet == 2) # 6 issue samples have blacktip reef info
sum(samples_issue$info_sheet == 3) # 8 issue samples occur in grey and blacktip 

## Evaluate info sheets
info_grey$fin_clip[info_grey$fin_clip == "?"] <- NA
sum(info_blacktip$fin_clip %in% info_grey$fin_clip)
sum(info_grey$fin_clip %in% info_blacktip$fin_clip)

info_grey$fin_clip[duplicated(info_grey$fin_clip)]
sum(info_grey$fin_clip != "", na.rm = T)
sum(info_grey$fin_clip %in% samples_all$fin_clip_id, na.rm = T)
sum(unique(info_grey$fin_clip) %in% samples_all$fin_clip_id, na.rm = T)


length(unique(info_grey$fin_clip)) # 911 - 1(NA) = 910 entries
length(unique(info_blacktip$fin_clip)) # 641  entries (no NAs or missing entries)

info_blacktip$fin_clip[duplicated(info_blacktip$fin_clip)]
sum(info_blacktip$fin_clip != "", na.rm = T)
sum(info_blacktip$fin_clip %in% samples_all$fin_clip_id, na.rm = T)
sum(unique(info_blacktip$fin_clip) %in% samples_all$fin_clip_id, na.rm = T)
