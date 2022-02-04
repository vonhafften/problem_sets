library(haven)
library(tidyverse)
library(plm)

setwd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps1/")

psid <- read_dta("pequiv_long.dta")

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
         age >= 25,
         age <= 25 + 34,
         income > 0)

quantile(sample$income, probs = c(0.01, 0.99), na.rm=TRUE)
         
sample <- sample %>%
  filter(income >= 1200,
         income <= 170000) 

# Estimate fixed effects for age and cohort
ols <- sample %>%
  lm(log(income) ~ factor(age) + factor(cohort), data = .)

# pull out age and cohort fixed effects
age_profile <- c(0, unname( ols$coefficients[2:35]))
cohort_profile <- c(0, unname(ols$coefficients[36:89]))

# plot age and cohort fixed effects
png("age.png")
plot(y = age_profile, x = 25:(25+34), xlab = "Age", ylab = "log(income)", main = "Life-Cycle Components")
dev.off()

png("cohort.png")
plot(y = cohort_profile, x = 1918:1972, xlab = "Cohort", ylab = "log(income)", main = "Cohort Components")
dev.off()

# save age profile for use in model
write_csv(tibble(log_income = age_profile + ols$coefficients[1]), "age_profile.csv")

########################################################################################
# estimate variance of persistent and transitory shock
########################################################################################

sample_2 <- cbind(sample, y = ols$residuals) 
sample_3 <- pdata.frame(sample_2, index = c("id", "year")) # use panel data structure for lags and leads

rho <- 0.97 # AR parameter

delta_y_tilde <- sample_3$y - rho*lag(sample_3$y) # find pseudo difference

# Formulas slide 79 of carters presentation
var_epsilon <- -(1/rho)*cov(delta_y_tilde, lead(delta_y_tilde), use = "complete.obs")
var_zeta <- (1/rho)*cov(delta_y_tilde, rho^2 * lag(delta_y_tilde) + rho*delta_y_tilde + lead(delta_y_tilde), use = "complete.obs")

print(var_epsilon)
print(var_zeta)

