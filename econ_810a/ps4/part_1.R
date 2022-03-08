library(haven)
library(tidyverse)
library(moments)
library(plm)

# model periods
setwd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/ps4/")

psid <- read_dta("../ps1/pequiv_long.dta")

# sample
# between 1978 and 1997
# between 25 and 55 (=25+30)
# income between 1000 and 150,000
sample <- psid %>% 
  transmute(id = x11101LL,
            year,
            age = d11101,
            cohort = year - age,
            income = i11103,
            is_seo = x11104LL == 12) %>%
  filter(year >= 1978,
         year <= 1997,
         !is_seo, 
         age >= 25,
         age <= 25+35,
         income > 1000) %>%
  mutate(log_income = log(income))

# Estimate fixed effects for age and cohort
ols <- sample %>%
  lm(log_income ~ factor(age) + factor(year) - 1, data = .)

summary(ols)

# set all years to base year to remove year effects
sample$year <- 1978

# use ols to predict log_income
sample$log_income_hat <- predict(ols, sample) + ols$residuals

annual_summary <- sample %>%
  group_by(age) %>%
  summarize(mean = mean(log_income_hat),
            sd = sd(log_income_hat),
            skewness = skewness(log_income_hat),
            kurtosis = kurtosis(log_income_hat)) %>%
  pivot_longer(-age) %>%
  mutate(name = factor(name, levels = c("mean", "sd", "skewness", "kurtosis")))

png("part_1.png", height = 360, width = 540)

annual_summary %>%
  ggplot(aes(x=age, y=value)) + 
  geom_line() +
  facet_wrap(~name, scales = "free_y")

dev.off() 

