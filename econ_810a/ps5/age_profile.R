library(haven)
library(tidyverse)
library(plm)

# model periods
paste0((0:12)*6, " - ", 5+(0:12)*6)

setwd("/Users/vonhafften/Documents/uw_madison/problem_sets/econ_810a/ps5/")

psid <- read_dta("../ps1/pequiv_long.dta")

# sample
# between 1978 and 1997
# between 25 and 60 (=25+35)
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
         age >= 24,
         age <= 77,
         income > 0)

# Estimate fixed effects for age and cohort
ols <- sample %>%
  lm(log(income) ~ factor(age) + factor(cohort), data = .)

# pull out age and cohort fixed effects
age_profile <- tibble(age = 24:77,
                      age_fe = ols$coefficients[1] + c(0, unname( ols$coefficients[2:54]))) %>%
  mutate(model_age = floor(age / 6)) %>%
  group_by(model_age) %>%
  summarize(kappa = mean(age_fe))

write_csv(age_profile, "age_profile.csv")

# plot age and cohort fixed effects
png("age.png")
plot(y = age_profile$kappa, x = 4:12, xlab = "Model Age", ylab = "log(income)", main = "Life-Cycle Components")
dev.off()
