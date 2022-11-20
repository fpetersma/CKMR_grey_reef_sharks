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

info_grey <- subset(info_grey, !is.na(fin_clip))


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

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Add info to usable samples
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Keep the 836 good grey samples
grey_samples <- subset(samples_good, info_sheet == 1) 

## There are 847 entries for 836 fin clips
sum(info_grey$fin_clip %in% grey_samples$sample_id) 

## Extract duplicate ids. These could indicate recaptures, and should be checked manually
usable_grey <- subset(info_grey, 
                      subset = fin_clip %in% grey_samples$sample_id,
                      select = c(date, location, new_recapture, time, sex, 
                                 PCL, FL, TL, fin_clip))
dup_ids <- subset(info_grey, subset = duplicated(fin_clip), select = fin_clip)
# View(subset(usable_grey, fin_clip %in% dup_ids$fin_clip))

## Only keep the ones with at most 1 new_capture, and where the recaptures make 
## chronological sense.
bad_duplicates <- c("1243", "2385", "F235", "F429", "F78")
usable_grey <- subset(usable_grey, !(fin_clip %in% bad_duplicates))
usable_grey$recapture_within_sample <- "N"
usable_grey[usable_grey$fin_clip %in% dup_ids$fin_clip, "recapture_within_sample"] <- 
  usable_grey[usable_grey$fin_clip %in% dup_ids$fin_clip, "new_recapture"]
usable_grey$recapture_within_sample <- 
  as.integer(gsub("N", 0, 
                  gsub("R", 1, usable_grey$recapture_within_sample)))
## Manually correct for fin clip 2339
usable_grey[usable_grey$fin_clip == "2339", "recapture_within_sample"] <- c(0, 1)

usable_grey$POSIXct <- as.POSIXct(usable_grey$date, format = "%m/%d/%y", tz = "HST")
usable_grey$year <- format(usable_grey$POSIXct, "%Y")
usable_grey$month <- format(usable_grey$POSIXct, "%m")
usable_grey$day <- format(usable_grey$POSIXct, "%d")

usable_grey$id <- 1:nrow(usable_grey)

usable_grey <- dplyr::select(usable_grey, c(id, fin_clip, location, POSIXct, 
                                            time, year, month, day, 
                                            recapture_within_sample, sex, 
                                            TL, PCL, FL))
news <- subset(usable_grey, recapture_within_sample == 0)

## Sample according to scheme discussed with Len and Danielle
set.seed(27102022)
## Sample up to 19 from every location
first_batch <- usable_grey %>% 
  group_by(location) %>% 
  mutate(loc_freq = n()) %>% # add a frequency column
  sample_n(min(loc_freq[1], 19)) # sample all or up to 19
## Now there are 191 samples, which is 3 too many. Take these 3 from locations
## North_Forereef, South_Forereef, an Southside.
reduced <- first_batch %>% 
  filter(location %in% c("North_Forereef", 
                         "South_Forereef", 
                         "Southside")) %>% 
  group_by(location) %>% 
  sample_n(18)

## Add the three reduced sampling locations back in first_batch
first_batch <- rbind(filter(first_batch, !(location %in% c("North_Forereef", 
                                                         "South_Forereef", 
                                                         "Southside"))),
                     reduced)

## How is the sex/year representation
with(first_batch, table(sex, year)) # Looks great

## Map info on the locations:
## https://products.coastalscience.noaa.gov/collections/benthic/e58palmyra/#horizontalTab3

## :::::::::::::::::::::::
## Replacing 'bad' samples
## :::::::::::::::::::::::

## 2232 was too small (only 5mg)
alternative <- usable_grey %>% 
  filter(location == "North_Forereef") %>%          # only keep North_Forereef
  filter(!(fin_clip %in% first_batch$fin_clip)) %>% # remove the ones that were initially selected
  filter(!(fin_clip %in% c(2000, 2033))) %>%        # sample 2000 is also to small (7mg)
  sample_n(1)                                       # sample the replacement 2041

## Replace 2232 by 2041
first_batch[first_batch$fin_clip == 2232, 1:13] <- 
  usable_grey[usable_grey$fin_clip == 2041, 1:13]
## Replace 2301 by 
first_batch[first_batch$fin_clip == 2301, 1:13] <- 
  usable_grey[usable_grey$fin_clip == XXX, 1:13]



## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Look up the black tips for Darcy
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Keep the 836 good blacktip samples
blacktips <- subset(samples_good, info_sheet == 2) 
View(blacktips[blacktips$box == "23XX", ])

