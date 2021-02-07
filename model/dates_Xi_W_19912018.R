rm(list = ls())
setwd("C:/Users/Philipp/Documents/GitHub/nowcasting-topics/model")
#_____________________________________________________#
#_This script creates a dataframe containing the dates
#_from January 1st 1991 to December 31st 2018 along 
#_with variables indicating the year, quarter, month,
#_week, day and weekday of said dates. Furthermore,
#_it also calculates the indicator variables Xi_q, 
#_Xi_m, Xi_w that indicate the start of a new quarter,
#_month or week. In addition, the weights for quarterly
#_flow variables W_q_c and W_q_p are constructed. 
#_These variables are used to generate
#_artificial data from a daily factor model mixing all
#_of the above frequencies. 
#_____________________________________________________#

# load packages
library(lubridate)
library(dplyr)
library(tidyr)

#_____________________________________________________#
#_set up dataframe with year, quarter, ... vars
#_____________________________________________________#

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


#_____________________________________________________#
#_construct vars that indicate start of period
#_____________________________________________________#

# create Xi indicators that equal 0 at start of period and 1 elsewhere
df %>% 
  mutate(Xi_q = ifelse(month == 1 & day == 1, 0, 1),
         Xi_m = ifelse(day == 1, 0, 1),
         Xi_w = ifelse(weekday == "Mo", 0, 1)) -> df


#_____________________________________________________#
#_calculate weights for quarterly flow variables
#_____________________________________________________#

# days per quarter and average number of days per quarter over the entire sample
df %>% 
  select(year, quarter, day) %>% 
  group_by(year, quarter) %>% 
  summarise(n_days_q = n()) %>% 
  ungroup() -> df_n_days_q

df_n_days_q$n_days_avg_q = floor(mean(df_n_days_q$n_days_q))

# loop over quarters to construct W_q_c and W_q_p
years <- seq(1991, 1991)
quarters <- seq(1, 4)
t_prev <- 0
t <- 0
k <- df_n_days_q$n_days_avg_q[1]
df_W_dq <- data.frame()
for (y in years)
{
  print(y)
  for (q in quarters)
  {
    print(paste("q =", q))
    t_prev <- t
    print(paste("tprev =", t_prev))
    k_t <- df_n_days_q[df_n_days_q$year == y & df_n_days_q$quarter == q, "n_days_q", drop = T]
    print(paste("k_t =",k_t))
    t <- t + k_t
    print(paste("t =", t))
    s <- seq(t_prev + 1, t)
    #W_qd_c = k * (t + 1 - s) / k_t # see Banbura et al. (2011, p. 30, eqn 10)
    W_qd_c = (t + 1 - s) / k_t # quarterly GDP level  is the average of daily GDP in the respective quarter
    if (y == years[length(years)] & q == quarters[length(quarters)]) # W_q_p = 0 => Not necessary but conceptually clearer!
    {

      df_W_q <- rbind(df_W_qd, data.frame(year = y,
                                         quarter = q,
                                         day = seq(1, k_t),
                                         W_qd_c = W_qd_c,
                                         W_qd_p = 0
                                        )
                    )

    } else
    {
      df_W_qd <- rbind(df_W_q, data.frame(year = y,
                                         quarter = q,
                                         day = seq(1, k_t),
                                         W_dq_c = W_qd_c,
                                         W_dq_p = c(0, rev(W_qd_c[2:length(W_q_c)]))
                                         )
                      )
    }
  }
}

# merge with df
df <- left_join(df, df_W_q, by = c("year", "quarter", "day"))

#_____________________________________________________#
#_dates for MATLAB plots
#_____________________________________________________#

# index and dates vars containing the year on the first day, else NA 
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

#_____________________________________________________#
#_export df to csv
#_____________________________________________________#

write.csv(df, "dates_Xi_W_19912018.csv", row.names = F, na = "")
