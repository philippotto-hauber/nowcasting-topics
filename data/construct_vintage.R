rm(list = ls())
setwd("C:/Users/Philipp/Documents/GitHub/nowcasting-topics/data")
#_____________________________________________________#
#_This script ...
#_____________________________________________________#

# PACKAGES ----
library(lubridate)
library(dplyr)
library(tidyr)
library(bundesbank)

# FUNCTIONS ----
#_____________________________________________________#
#_f_interpolate (currently not required),
#_ roll_mean 
#_ bw_filter
#_ f_outl
#_____________________________________________________#

f_interpolate <- function(x)
# replace NA in x with linear interpolation
{
  for (t in seq(1, length(x)))
  {
    if (is.na(x[t]))
    {
      ind_nextobs <- 0
      t_NA = T
      while (t_NA == T)
      {
        
        if (is.na(x[t+ind_nextobs]))
          ind_nextobs <- ind_nextobs + 1
        else
          t_NA = F
        
      }
      x[t] <- 0.5 * (x[t-1] + x[t+ind_nextobs])
    }
  }
  return(x)
}

rollmean <- function(x, k){
  xroll <- array(NA, c(length(x)))
  for (t in seq(k, length(x)))
    xroll[t] <- mean(x[(t-k+1):t], na.rm = TRUE)
  
  return(xroll)
}

bw_filter <- function(y, bw)
{
  # compute un-normalized weights
  j <- seq(-bw, bw, 1) 
  omeg = (1 - (j/bw) ^ 2) ^ 2  
  
  # check for NA's in y
  n_NA <- sum(is.na(y)) 
  y_noNA <- y[(n_NA + 1) : length(y)]
  Nt <- length(y_noNA)
  
  # loop over t
  tau <- mat.or.vec(length(y_noNA), 1)
  for (t in 1 : length(y_noNA)) {
    # case distinction
    if (t <= bw) {
      
      indY <- c(1 : (2 * bw - (bw - t)))
      indOmeg <- c((bw - t + 1):(2 * bw)) 
      kap <- 1 / ( sum(omeg[indOmeg]))
      tau[ t ] <- kap * omeg[indOmeg] %*% y_noNA[indY] 
      
    } else if (t > (Nt - bw)) {
      
      indY <- c((t - bw) : Nt)
      indOmeg <- c(1 : (bw + 1 + (Nt - t)))
      kap <- 1 / ( sum( omeg[indOmeg] ) )
      tau[t] <- kap * omeg[indOmeg] %*% y_noNA[indY] 
      
    } else {
      
      indY <- seq(t - bw, t + bw, 1)
      indOmeg <- c( 1 : (2 * bw + 1))
      kap <- 1 / (sum(omeg[indOmeg]))
      tau[t] <- kap * omeg[indOmeg] %*% y_noNA[indY]  
    }
  }
  return(c(rep(NA, times = n_NA), tau))
}

f_outl <- function(y, aalpha)
{
  return(abs((y - median(y, na.rm = T))) > aalpha * IQR(y, na.rm = T))
}

# SET-UP ----
#_____________________________________________________#
#_specify vintage, 
#_forecast horizon,
#_number of topics to select
#_____________________________________________________#

vintage <- c("2010-01-30") # requires backcast!
date_h <- c("2010-06-30") # 2010Q2
Ntopics <- 10 # select Ntopics topics that have the highest correlation with quarterly GDP growth

# LOAD TOPICS ----
#_____________________________________________________#
#_select sample,
#_extend to 7-day week,
#_linearly interpolate to fill-in NA (commented out because not necessary)
#_____________________________________________________#

df_raw <- read.csv("./topics/daily_topics.csv")

# add date and quarter variable
df_raw %>%
  mutate(date = make_date(year = df_raw$year, 
                          month = df_raw$month, 
                          day = df_raw$day),
         quarter = ceiling(month / 3)) %>% 
  filter(date <= vintage) %>%
  select(date, year, quarter, month, day, everything()) -> df_topics 

# get rid of raw df
rm(df_raw)

# extend series with NA to 7-day week 
dates_tmp <- data.frame(date = seq(min(df_topics$date), 
                           max(df_topics$date), 
                           by = "days")
                )

dates_tmp %>% 
  mutate(year = year(date),
         month = month(date),
         quarter = ceiling(month / 3),
         day = day(date)) %>%
  merge(df_topics, by = c("date", "year", "quarter", "month", "day"), all.x = T) -> df_topics

# get rid of dates_tmp
rm(dates_tmp)

# col indices corresponding to topics
ind_topics <- which(grepl("T", names(df_topics)))

# linear interpolate topics to fill-in missings, storing pattern of NA => commented out because moving average should take care of missings
# ind_NA <- c()
# for (n in ind_topics)
# {
#   x <- df_topics[, n]
#   ind_NA <- cbind(ind_NA, is.na(x))
#   x <- f_interpolate(x)
#   df_topics[, n] <- x
# }
# 
# colnames(ind_NA) <- names(df_topics)[grepl("T", names(df_topics))]

# TRANSFORM TOPCIS ---- 
#_____________________________________________________#
#_rm outlier,
#_moving average,
#_detrend using biweight filter
#_reimpose NA pattern
#_____________________________________________________#

# remove outlier
ind_outl <- apply(df_topics[, ind_topics], c(2), f_outl, aalpha = 10)
dat <- df_topics[, ind_topics]
dat[ind_outl] <- NA

# store pattern of missings
ind_NA <- is.na(dat)
colnames(ind_NA) <- names(df_topics)[grepl("T", names(df_topics))]

# moving average
dat_ma <- apply(dat, c(2), rollmean, k = 60)

# detrend with biweight filter
dat_bw <- apply(dat_ma, c(2), bw_filter, bw = 1200)
dat_detrend <- dat_ma - dat_bw # de-trended topics

# reimpose NA pattern
dat_detrend_NA <- dat_detrend
dat_detrend_NA[ind_NA] <- NA

# store in df_topics_trafo
df_topics_trafo <- df_topics
df_topics_trafo[ind_topics] <- dat_detrend_NA


plot(df_topics$date, df_topics[, 7], 
     type = "l", col = "blue", 
     ylim = c(-0.02, 0.12), 
     main = "Topic T1", sub = "raw and detrended series",
     ylab = "",
     xlab = "")
lines(df_topics$date, df_topics_trafo[, 7], type = "l", col = "red")



# GDP DATA ----

