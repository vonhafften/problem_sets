log using analysis, replace

* Alex von Hafften
* Problem set 1
* ECON 717A: Applied Economics

* clear workspace
clear

* install user defined functions (if needed)
ssc install outreg2

* change working directory
cd "/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_717a/ps1/"

* open dataset
use "Field et al (2010) Analysis Sample"

********************************************************************************
* problem #1 - drop missing values
********************************************************************************

drop if missing(HH_Income)
drop if miss_Client_Age == 1
drop if miss_Client_Married == 1
drop if miss_Client_Education == 1

********************************************************************************
* problem #2 - estimate linear probability model with homoskedastic standard errors
********************************************************************************

* define list of variable in usual_covariate
local covariates " Client_Married Client_Education HH_Size HH_Income muslim Hindu_SC_Kat Treated"

* estimate LPM
regress taken_new Client_Age `covariates'
outreg2 using p2_table, tex(frag) replace

********************************************************************************
* problem #3 - estimate linear probability model with heteroskedastic standard errors
********************************************************************************

regress taken_new Client_Age `covariates'
outreg2 using p3_table, tex(frag) replace addtext(Comma Robust SEs, No)
regress taken_new Client_Age `covariates', robust
outreg2 using p3_table, tex(frag) append addtext(Comma Robust SEs, Yes)

********************************************************************************
* problem #4 - predicted probabilities
********************************************************************************

predict taken_new_hat_lpm
histogram taken_new_hat_lpm
graph export p4_figure.png, replace
count if taken_new_hat_lpm > 1
count if taken_new_hat_lpm < 0

********************************************************************************
* problem #5 - weighted least squares
********************************************************************************

vwls taken_new `covariates'

* Doesn't work due to small sample, but works for more parsimonious model
regress taken_new Client_Age Client_Education HH_Size  Treated, robust
outreg2 using p5_table, tex(frag) addtext(Weighted By, None) replace
vwls taken_new Client_Age Client_Education HH_Size  Treated 
outreg2 using p5_table, tex(frag) addtext(Weighted By, Variance) append

********************************************************************************
* problem #6 - probit and logit
********************************************************************************

regress taken_new Client_Age `covariates', robust
outreg2 using p6_table, tex(frag) replace addtext(Model, LPM)
probit taken_new Client_Age `covariates'
outreg2 using p6_table, tex(frag) append addtext(Model, Probit)
logit taken_new Client_Age `covariates'
outreg2 using p6_table, tex(frag) append addtext(Model, Logit)


********************************************************************************
* problem #7 - mean partial derivatives of Client_Age
********************************************************************************

* LPM 
regress taken_new Client_Age `covariates'
* mean partial derivative is just the LPM coefficient = -.0000283

* probit 
* part a - using dprobit
dprobit taken_new Client_Age `covariates'
* mean partial is the value that dprobit outputs = .0000382

* part b - analytical derivative
probit taken_new Client_Age `covariates'
* get the linear prediction based on probit
predict taken_new_hat_xb, xb
* using formula phi(xb)*b_j
gen Client_Age_Partial_a = normalden(taken_new_hat_xb) * e(b)[1,1]
summarize Client_Age_Partial_a
* mean partial derivative is the mean = .0000382

* part c - numerically calculating marginal effects
* predict probability based on probit
predict taken_new_hat_probit, pr
gen taken_new_hat_probit_epsilon = normal(taken_new_hat_xb + 0.001*e(b)[1,1])

* compute numerical derivative
gen Client_Age_Partial_n =  (taken_new_hat_probit_epsilon - taken_new_hat_probit) / 0.001
summarize Client_Age_Partial_n
* mean partial derivative is the mean = 0

* part d - using margins
probit taken_new Client_Age `covariates'
margins , dydx(Client_Age) atmeans
* mean partial derivative is the mean =  0.0000382

* logit
logit taken_new Client_Age `covariates'
margins , dydx(Client_Age) atmeans
* mean partial derivative is the mean =   -.0000525

********************************************************************************
* problem #8 - LPM with quadratic age
********************************************************************************

* baseline
regress taken_new `covariates', robust
outreg2 using p8_table, tex(frag) replace

* create quadratic transformations of age
gen Client_Age_2 = Client_Age^2 
gen Client_Age_3 = Client_Age^3 
gen Client_Age_4 = Client_Age^4 

* estimate lpm with quadratic transformations of age
regress taken_new Client_Age Client_Age_2 Client_Age_3 Client_Age_4 `covariates', robust
outreg2 using p8_table, tex(frag) append
predict taken_new_hat_lpm_q

* check for observations outside 0, 1
count if taken_new_hat_lpm_q > 1
count if taken_new_hat_lpm_q < 0
histogram taken_new_hat_lpm_q
graph export p8_figure.png, replace

* create quadratic transformations with epsilon
gen Client_Age_epsilon_2 = Client_Age_epsilon^2 
gen Client_Age_epsilon_3 = Client_Age_epsilon^3 
gen Client_Age_epsilon_4 = Client_Age_epsilon^4 
regress taken_new Client_Age_epsilon Client_Age_epsilon_2 Client_Age_epsilon_3 Client_Age_epsilon_4 `covariates'
predict taken_new_hat_lpm_q_epsilon

* numerically compute derivative
gen Client_Age_Partial_q_n =  (taken_new_hat_lpm_q - taken_new_hat_lpm_q_epsilon) / 0.001
summarize Client_Age_Partial_q_n
* mean is 2.32e-07

********************************************************************************
* problem #9 - LRI
********************************************************************************

* baseline probit
probit taken_new Client_Age `covariates'
* ll_hat is  -238.07469

* probit only with constant
probit taken_new
* ll_0 is -240.23429

* lri = 1 - ll_hat/ll_0 = 0.00898955765

********************************************************************************
* problem #10 - Prediction rate
********************************************************************************

* using cutoff of 50 percent
gen predicted_over_50 =  taken_new_hat_probit > .5
tab taken_new predicted_over_50, nofreq cell

* using unconditional probability as cutoff
tab taken_new
* unconditional probability = .1673
gen predicted_over_up =  taken_new_hat_probit > .1673
tab taken_new predicted_over_up, nofreq cell

********************************************************************************
* problem #11 - In sample vs. out-of-sample prediction
********************************************************************************

gen estimation_sample = imidlineid < 1400

* estimate probit on subsample
probit taken_new Client_Age `covariates' if estimation_sample
predict taken_new_hat_probit_11, pr

* using cutoff of 50 percent
gen predicted_over_50_11 =  taken_new_hat_probit_11 > .5
tab taken_new predicted_over_50_11 if estimation_sample, cell nofreq
tab taken_new predicted_over_50_11 if !estimation_sample, cell  nofreq

* using unconditional probability as cutoff
gen predicted_over_up_11 =  taken_new_hat_probit_11 > .1673
tab taken_new predicted_over_up_11 if estimation_sample, cell nofreq
tab taken_new predicted_over_up_11 if !estimation_sample, cell nofreq

********************************************************************************
* problem #12 - Interaction terms
********************************************************************************

probit taken_new Client_Age `covariates'
outreg2 using p12_table, tex(frag) replace

gen married_muslim = Client_Married * muslim
probit taken_new Client_Age `covariates' married_muslim
outreg2 using p12_table, tex(frag) append

********************************************************************************
* problem #13 - Interaction effects
********************************************************************************

* compute interaction effect without accounting for terms in Ai and Norton (2003)
margins , dydx(married_muslim)
* interaction effect estimate is -.0671875

* compute interaction effect by hand accounting for terms in Ai and Norton (2003)
* follows logic from lecture notes
predict index_hat, xb
* predicted index with both dummies zero
gen index_hat_0 = index_hat - Client_Married * e(b)[1,2] - muslim * e(b)[1,6] - married_muslim *e(b)[1, 9]
* predicted index with both married one and muslim zero
gen index_hat_01 = index_hat_0 + e(b)[1,2]
* predicted index with both muslim one and married zero
gen index_hat_02 = index_hat_0 + e(b)[1,6]
* predicted index with both dummies zero
gen index_hat_012 = index_hat_0 + e(b)[1,2] + e(b)[1,6] + e(b)[1,9]
gen finite_difference = (normal(index_hat_012) - normal(index_hat_02)) - (normal(index_hat_01) - normal(index_hat_0))
summarize finite_difference
* interaction effect estimate is -.0656881

********************************************************************************
* problem #14 - Interaction effects variance
********************************************************************************

* see summarize table from problem #13.

********************************************************************************
* problem #15 - Heteroskedasticity test
********************************************************************************

* compute residuals
regress taken_new Client_Age `covariates' 
predict residuals_p, residuals

* regress squared residuals on usual covariates.
gen residuals_p_2 = residuals_p^2
regress residuals_p_2 Client_Age `covariates' 
outreg2 using p15_table, tex(frag) replace


********************************************************************************
* problem #16 - hetprob
********************************************************************************

hetprob taken_new Client_Age `covariates', het(Client_Age Client_Education)
outreg2 using p16_table, tex(frag) replace

log close

translate analysis.smcl analysis.log, linesize(255) replace
