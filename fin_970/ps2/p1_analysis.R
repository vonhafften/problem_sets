library(tidyverse)
library(readxl)
library(lubridate)
library(stargazer)

setwd("Documents/UW Madison/problem_sets/fin_970/ps2/")

data <- read_excel("p1_data.xlsx", skip = 1) %>%
  as_tibble() %>%
  filter(month(Date) == 12) %>% # reduce to annual data
  transmute(date = year(Date), pd, r, delta_d) %>% # pull out variables to use below
  mutate(delta_d_1 = lead(delta_d),
         delta_d_2 = (lead(delta_d) + lead(delta_d, 2L))/2,
         delta_d_3 = (lead(delta_d) + lead(delta_d, 2L) + lead(delta_d, 3L))/3,
         delta_d_4 = (lead(delta_d) + lead(delta_d, 2L) + lead(delta_d, 3L) + lead(delta_d, 4L))/4,
         delta_d_5 = (lead(delta_d) + lead(delta_d, 2L) + lead(delta_d, 3L) + lead(delta_d, 4L) + lead(delta_d, 5L))/5,
         r_1 = lead(r),
         r_2 = (lead(r) + lead(r, 2L))/2,
         r_3 = (lead(r) + lead(r, 2L) + lead(r, 3L))/3,
         r_4 = (lead(r) + lead(r, 2L) + lead(r, 3L) + lead(r, 4L))/4,
         r_5 = (lead(r) + lead(r, 2L) + lead(r, 3L) + lead(r, 4L) + lead(r, 5L))/5)

# dividend regressions
delta_d_1_lm <- lm(delta_d_1 ~ pd, data=data)
delta_d_2_lm <- lm(delta_d_2 ~ pd, data=data)
delta_d_3_lm <- lm(delta_d_3 ~ pd, data=data)
delta_d_4_lm <- lm(delta_d_4 ~ pd, data=data)
delta_d_5_lm <- lm(delta_d_5 ~ pd, data=data)

# return regressions
r_1_lm <- lm(r_1 ~ pd, data=data)
r_2_lm <- lm(r_2 ~ pd, data=data)
r_3_lm <- lm(r_3 ~ pd, data=data)
r_4_lm <- lm(r_4 ~ pd, data=data)
r_5_lm <- lm(r_5 ~ pd, data=data)

sink("p1_delta_d_table.tex")
stargazer(delta_d_1_lm, delta_d_2_lm, delta_d_3_lm, delta_d_4_lm, delta_d_5_lm, 
          omit.stat = c("f", "ser", "adj.rsq"),
          title = "Dividends")
sink()

sink("p1_r_table.tex")
stargazer(r_1_lm, r_2_lm, r_3_lm, r_4_lm, r_5_lm, 
          omit.stat = c("f", "ser", "adj.rsq"),
          title = "Market Returns")
sink()

# for comparison in part 4
mean(data$pd)


