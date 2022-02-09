# Alex von Hafften
# ECON 810: Advanced Macro
# Carter Braxton

# Problem set 2

library(haven)
library(tidyverse)
library(plm)
library(sscore)

setwd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps2/")

psid <- read_dta("../ps1/pequiv_long.dta")

sample <- psid %>% 
  transmute(id = x11101LL,
            year,
            age = d11101,
            cohort = year - age,
            income = i11103,
            is_seo = x11104LL == 12,
            hours = e11101)%>%
  filter(year >= 1978,
         year <= 1997,
         !is_seo, 
         age >= 25,
         age <= 25 + 34,
         income > 0,
         income < exp(12),
         hours > 50*36)

#hist(sample$income)

#hist(sample$hours, breaks = 60)
#abline(v=50*36, col ="red")

sample_panel <- pdata.frame(sample, index = c("id", "year")) 

sample_panel$income_lag <- lag(sample_panel$income)
sample_panel$income_change <- (sample_panel$income - sample_panel$income_lag)/sample_panel$income_lag

# limits to income changes that are less than doubling.
sample_panel <- sample_panel %>%
  filter(income_change < 1)

hist(sample_panel$income_change)

# equals 0.0735215
mean(sample_panel$income_change, na.rm =TRUE, trim = 0.01)
