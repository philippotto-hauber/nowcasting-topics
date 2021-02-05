rm(list = ls())
setwd("C:/Users/Philipp/Documents/GitHub/nowcasting-topics/model")
#_____________________________________________________#
#_This script creates a dataframe containing the dates
#_from January 1st 1991 to December 31st 2018 along 
#_with variables indicating the year, quarter, month,
#_week, day and weekday of said dates. Furthermore,
#_it also calculates the indicator variables Xi_q, 
#_Xi_m, Xi_w that indicate the start of a new quarter,
#_month or week. These variables are used to generate
#_artifical data from a daily factor model mixing all
#_of the above frequencies. 
#_____________________________________________________#

# load packages
library(lubridate)
library(dplyr)
library(tidyr)

# date variable from 1.1.1991 to 31.12.2018
df <- data.frame(date = seq(
                            as.Date("1991-01-01"), 
                            as.Date("2018-12-31"), 
                            by = "day"
                            )
                  )

# add year, quarter, month, week, day, weekday vars  
df %>% 
  mutate(year = year(date),
         month = month(date),
         quarter = ceiling(month/3),
         week = isoweek(date),
         day = day(date),
         weekday = wday(date, label = T)) %>%
  select(date, weekday, year, quarter, month, week, day) -> df



# create Xi indicators that equal 0 at start of period and 1 elsewhere
df %>% 
  mutate(Xi_q = ifelse(month == 1 & day == 1, 0, 1),
         Xi_m = ifelse(day == 1, 0, 1),
         Xi_w = ifelse(weekday == "Mo", 0, 1)) -> df

# create index and dates vars containing the year on the first day, else NA
df$ind_plot <- NA
df$dates_plot <- NA
for (t in seq(1, nrow(df)))
{
  if (df$quarter[t] == 1 && df$month[t] == 1 && df$day[t] == 1)
  {
    df$dates_plot[t] = df$year[t]
    df$ind_plot[t] = t
  }
}

# export df to csv
write.csv(df, "dates_Xi_19912018.csv", row.names = F, na = "")
