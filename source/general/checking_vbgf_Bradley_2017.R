## ========================================================================== ##
## In this script, I redo the caluclations from Bradley 2017 in PLoS One      ##
## since it seems like the VBGF growth curve that they estimated makes no     ##
## sense. Sharks grow way too slow and mature after they ahve passed they     ##
## maximum age.                                                               ##
## ========================================================================== ##

## Load the data ============================================================ ##
df <- read.csv("C:/Users/felix/Downloads/data_bradley_2017_PLoS1.csv")

## Calculate the differences
df$diff_cm <- df$TL.2_.cm. - df$TL.1_.cm. # derive growth
df$Tag_Date <- as.POSIXct(df$Tag_Date, format = "%m/%d/%y") # dates as POSIXct
df$Recapture_Date <- as.POSIXct(df$Recapture_Date, format = "%m/%d/%y")
df$diff_date <- as.numeric(df$Recapture_Date - df$Tag_Date) # calculate elapsed time in days

## Plot the date to understand the problem (lots of negative growth!)
plot(as.numeric(df$diff_date), as.numeric(df$diff_cm))
abline(lm(df$diff_cm ~ df$diff_date), col = "red")

## Fit a VBGF model to the data ============================================= ##
## Use the fishmethods package
library(fishmethods)

## Francis *(1988) specification of VBGF via function grotag()
?grotag
## grotag() requires julian days, so convert dates to years since 1970-01-01 00:00
df$Tag_Julian <- as.numeric(df$Tag_Date) / 86400 / 365
df$Recapture_Julian <- as.numeric(df$Recapture_Date) / 86400 / 365

# Fit
vbgf <- grotag(
  L1 = df$TL.1_.cm.,
  L2 = df$TL.2_.cm.,
  T1 = df$Tag_Julian,
  T2 = df$Recapture_Julian,
  alpha = 100,
  beta = 130,
  design = list(nu = 1,    # v in paper, growth variability
                m = 0,     # m in paper, bias/mean of measurement error
                p = 0,     # p in paper, outlier probability/contamination
                sea = 0)   # u and w in paper, seasonal variation
)                     

vbgf$table # results match the paper
vbgf$VBparms

l_inf <- vbgf$VBparms[2,1]
r <- vbgf$VBparms[3,1]

## How to derive a0/t0? Assume l0 = 60cm (length at age 0), then we can derive 
l_0 = 60
a_0 = log(1 - l_0 / l_inf) / r # -8.27

## What is the max age?
# Maximum lifespan Tmax for sexes combined was then estimated as the time 
# required to attain >99% of TLinf as Tmax = 5*ln(2)*k-1"
## Comes from Fabens 1965 and bradley 2017
## In my notation it is 
a_max = 5 * log(2) / r # 62.5 years

## Plot this vbgf
ages <- 0:65
lengths <- CKMRcpp::vbgf(ages, a_0, r, l_inf)
plot(ages, lengths, type = "l")

## What if we remove the negative measurements? I know the measurement error is 
## supposed to scoop up the negative growth, but i hasn't really since the 
## predicted values still contain a lot of negative growth estimates... 

df <- df[df$diff_cm >= 0, ]

## Doesn't seem to change much. Overall, it seems that they did the fitting
## of the model to the data correctly, but that the data is just a weird sample.
## Honestly, considering how many sharks "shrunk", I am unsure whether to trust
## the data, which is worrying: what if the measurement error is non-random?

l_mature_male <- (boot::logit(0.5) + 26.02) / 0.21 # from clasper glm in bradley2017
l_mature_female <- 0.7 * 175.5 + 3.29 # Frisk, Miller, and Fogarty (2001) with TL_female from Bradley2017
alpha_m <- CKMRcpp::invvbgf(l_mature_male, a_0, r, l_inf) # 17.44 = 17
alpha_f <- CKMRcpp::invvbgf(l_mature_female, a_0, r, l_inf) # 18.5 = 19
## refitting the model gives L_inf = 175.26 and k = 0.045

## Use FSA package from fishR
library(FSA)

Francis <- vbFuns("Francis")
