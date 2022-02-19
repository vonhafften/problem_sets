# Alex von Hafften
# ECON 810: Advanced Macro
# Professor Carter Braxton

# Problem Set 2

library(haven)
library(tidyverse)
library(plm)
library(fixest)

setwd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps3/")

psid <- read_dta("../ps1/pequiv_long.dta")

# initial filtering
sample <- psid %>% 
  transmute(id = x11101LL,
            year,
            age = d11101,
            cohort = year - age,
            income = i11103,
            is_seo = x11104LL == 12,
            hours = e11101,
            size = d11106) %>%
  filter(year >= 1978,
         year <= 1997,
         !is_seo, 
         age >= 25,
         age <= 65)

# for nicely behaved lags
sample_panel <- pdata.frame(sample, index = c("id", "year")) 

# hours lags
sample_panel$hours_lag_1 <- plm::lag(sample_panel$hours, k=1L)
sample_panel$hours_lag_2 <- plm::lag(sample_panel$hours, k=2L)
sample_panel$hours_lag_3 <- plm::lag(sample_panel$hours, k=3L)

# HH size lags
sample_panel$size_lag <- plm::lag(sample_panel$size, k=1L)

# income lags/leads
sample_panel$income_lead_5 <- plm::lead(sample_panel$income, k=5L)
sample_panel$income_lead_4 <- plm::lead(sample_panel$income, k=4L)
sample_panel$income_lead_3 <- plm::lead(sample_panel$income, k=3L)
sample_panel$income_lead_2 <- plm::lead(sample_panel$income, k=2L)
sample_panel$income_lead_1 <- plm::lead(sample_panel$income, k=1L)
sample_panel$income_lag_1 <- plm::lag(sample_panel$income, k=1L)
sample_panel$income_lag_2 <- plm::lag(sample_panel$income, k=2L)
sample_panel$income_lag_3 <- plm::lag(sample_panel$income, k=3L)
sample_panel$income_lag_4 <- plm::lag(sample_panel$income, k=4L)

# Definition of treatment

sample_2 <- sample_panel %>%
  filter(hours_lag_1 >= 50*36,
         hours_lag_2 >= 50*36,
         hours_lag_3 >= 50*36,
         size == size_lag) %>%
  mutate(treatment = (hours < (hours_lag_1 + hours_lag_2 + hours_lag_3)/3 * .75))

summary(sample_2$treatment)

# Reorganize data by each lead/lag
sample_lag_4 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income_lag_4,
            treatment_lag_4 = as.integer(treatment),
            k=-4) %>%
  filter(!is.na(income))


sample_lag_3 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income_lag_3,
            treatment_lag_3 = as.integer(treatment),
            k=-3) %>%
  filter(!is.na(income))

sample_lag_2 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income_lag_2,
            treatment_lag_2 = as.integer(treatment),
            k=-2)%>%
  filter(!is.na(income))


sample_lag_1 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income_lag_1,
            treatment_lag_1 = as.integer(treatment),
            k=-1)%>%
  filter(!is.na(income))


sample_0 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income,
            treatment_0 = as.integer(treatment),
            k=0)%>%
  filter(!is.na(income))


sample_lead_1 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income_lead_1,
            treatment_lead_1 = as.integer(treatment),
            k=1)%>%
  filter(!is.na(income))


sample_lead_2 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income_lead_2,
            treatment_lead_2 = as.integer(treatment),
            k=2)%>%
  filter(!is.na(income))

sample_lead_3 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income_lead_3,
            treatment_lead_3 = as.integer(treatment),
            k=3)%>%
  filter(!is.na(income))

sample_lead_4 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income_lead_4,
            treatment_lead_4 = as.integer(treatment),
            k=4)%>%
  filter(!is.na(income))

sample_lead_5 <- sample_2 %>%
  transmute(id,
            year,
            age,
            cohort,
            income = income_lead_5,
            treatment_lead_5 = as.integer(treatment),
            k=5)%>%
  filter(!is.na(income))

# recombine data into single dataset

analysis_sample <- bind_rows(sample_lag_4, sample_lag_3,sample_lag_2,sample_lag_1, sample_lag_0, sample_lead_1, sample_lead_2, sample_lead_3, sample_lead_4, sample_lead_5)

analysis_sample[is.na(analysis_sample)] <- 0

distributed_lags <- feols(income ~ treatment_lag_4 + treatment_lag_3 + treatment_lag_2 + treatment_lag_1 + treatment_0 + treatment_lead_1 + treatment_lead_2 + treatment_lead_3 + treatment_lead_4 + treatment_lead_5 | c(k, id), data = analysis_sample)

summary(distributed_lags)

png("part_1.png",width = 480, height = 320)
plot(x = -4:5, y = distributed_lags$coefficients, ylim = c(-14000, 0), col = "blue", main = "Earnings Loss Around Job Loss", ylab = "$", xlab = "Period Relative to Job Loss")
lines(x = -4:5, y = distributed_lags$coefficients, col = "blue")
abline(h=0)
dev.off()
