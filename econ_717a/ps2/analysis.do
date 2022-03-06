* ECON 717A: Applied Econometrics
* Problem set 2
* Professor: Jeff Smith
* Alex von Hafften
* Matching and weighting

* clear workspace
clear

* install user defined functions (if needed)
ssc install outreg2
ssc install psmatch2
ssc install texdoc
ssc install stddiff

* change working directory
cd "/Users/vonhafften/Documents/UW Madison/problem_sets/econ_717a/ps2/"

* open dataset
use "Economics 717 Spring 2022 NSW Data.dta"

********************************************************************************
* problem #0 - drop observation from the PSID
********************************************************************************

drop if sample == 3

********************************************************************************
* problem #1 - treatment effect from experimental data
********************************************************************************

gen age_2 = age^2

regress re78 treated, robust
outreg2 using table_1, tex(frag) replace

regress re78 treated age age_2 educ black hisp married nodegree re74 re75, robust
outreg2 using table_1, tex(frag) append

********************************************************************************
* problem #2 - drop the experimental treatment group
********************************************************************************

drop if treated == 1  & sample==1 

********************************************************************************
* problem #3 - estimate propensity scores
********************************************************************************

gen in_control = (sample ==1) 

probit in_control age age_2 educ black hisp married nodegree
outreg2 using table_3, tex(frag) replace addstat(Comparison group obs. completely determined, e(N_cdf), Control group obs. completely determined,  e(N_cds))
predict pscorea, pr

probit in_control age age_2 educ black hisp married nodegree re74 re75
outreg2 using table_3, tex(frag) append addstat(Comparison group obs. completely determined, e(N_cdf), Control group obs. completely determined,  e(N_cds))
predict pscoreb, pr

********************************************************************************
* problem #4 - compare pscorea and pscoreb descriptive statistics
********************************************************************************

est clear
estpost tabstat pscorea, by(in_control) c(stat) stat(mean, sd, min, median, max, count)
esttab, cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) p50(fmt(%13.2fc)) max(fmt(%13.2fc))  count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Median" "Max" "N")
esttab using "table_4a.tex", replace cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) p50(fmt(%13.2fc)) max(fmt(%13.2fc))  count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Median" "Max" "N")

est clear
estpost tabstat pscoreb, by(in_control) c(stat) stat(mean, sd, min, median, max, count)
esttab, cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) p50(fmt(%13.2fc)) max(fmt(%13.2fc))  count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Median" "Max" "N")
esttab using "table_4b.tex", replace cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) p50(fmt(%13.2fc)) max(fmt(%13.2fc))  count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Median" "Max" "N")

********************************************************************************
* problem #5 - histograms
********************************************************************************

* create bins based on pscorea and pscoreb
egen binsa=cut(pscorea), at(0(.05)1) icodes
egen binsb=cut(pscoreb), at(0(.05)1) icodes

* histogram of pscorea
graph bar (count) pscorea, over(in_control) over(binsa, label(nolab)) asyvars title("Histogram of pscorea")
graph export figure_5a.png, replace

* histogram of pscoreb
graph bar (count) pscoreb, over(in_control) over(binsb, label(nolab)) asyvars title("Histogram of pscoreb")
graph export figure_5b.png, replace

* histogram of pscorea wo lowest bin
graph bar (count) pscorea if binsa > 0, over(in_control) over(binsa, label(nolab)) asyvars title("Histogram of pscorea wo lowest bin")
graph export figure_5a_2.png, replace

* histogram of pscoreb wo lowest bin
graph bar (count) pscoreb if binsb > 0, over(in_control) over(binsb, label(nolab)) asyvars title("Histogram of pscoreb wo lowest bin")
graph export figure_5b_2.png, replace
 
********************************************************************************
* problem #6 - Nearest neighbor wo replacement
********************************************************************************

global table_number = 6

est clear 
eststo: psmatch2 in_control, noreplacement outcome(re78) pscore(pscorea) neighbor(1) common

global att_coarse = r(att)
global att_coarse_se = r(seatt)

esttab, se

global unm_coef = r(coefs)[1, 1]
global unm_se = r(coefs)[1, 2]

psmatch2 in_control, noreplacement outcome(re78) pscore(pscoreb) neighbor(1) common

global att_fine = r(att)
global att_fine_se = r(seatt)

texdoc do table_maker

********************************************************************************
* problem #7 - nearest neighbor w replacement
********************************************************************************

global table_number = 7

est clear 
eststo: psmatch2 in_control, outcome(re78) pscore(pscorea) neighbor(1) common

global att_coarse = r(att)
global att_coarse_se = r(seatt)

esttab, se

global unm_coef = r(coefs)[1, 1]
global unm_se = r(coefs)[1, 2]

psmatch2 in_control, outcome(re78) pscore(pscoreb) neighbor(1) common

global att_fine = r(att)
global att_fine_se = r(seatt)

texdoc do table_maker

********************************************************************************
* problem #8 - standardized differences
********************************************************************************

* raw data
stddiff re74 re75, by(in_control)

* rich pscore nearest neighbor - re74
psmatch2 in_control, outcome(re74) pscore(pscoreb) neighbor(1) common

* rich pscore nearest neighbor - re74
psmatch2 in_control, outcome(re75) pscore(pscoreb) neighbor(1) common

********************************************************************************
* problem #9 - gaussian kernel matching
********************************************************************************

global table_number = 9

psmatch2 in_control, kernel outcome(re78) kerneltype(normal) pscore(pscoreb) bwidth(0.02) common

global att_low = r(att)
global att_low_se = r(seatt)

psmatch2 in_control, kernel outcome(re78) kerneltype(normal) pscore(pscoreb) bwidth(0.2) common

global att_med = r(att)
global att_med_se = r(seatt)

psmatch2 in_control, kernel outcome(re78) kerneltype(normal) pscore(pscoreb) bwidth(2.0) common

global att_hi = r(att)
global att_hi_se = r(seatt)

texdoc do table_maker_2

********************************************************************************
* problem #10 - local linear matching
********************************************************************************

global table_number = 10

psmatch2 in_control, llr outcome(re78) pscore(pscoreb) bwidth(0.02) common

global att_low = r(att)
global att_low_se = r(seatt)

psmatch2 in_control, llr outcome(re78) pscore(pscoreb) bwidth(0.2) common

global att_med = r(att)
global att_med_se = r(seatt)

psmatch2 in_control, llr outcome(re78) pscore(pscoreb) bwidth(2.0) common

global att_hi = r(att)
global att_hi_se = r(seatt)

texdoc do table_maker_2

********************************************************************************
* problem #11 - linear regression of in_control
********************************************************************************

regress re78 in_control age age_2 educ black hisp married nodegree re74 re75, robust
outreg2 using table_11a, tex(frag) replace

predict y_hat_11

est clear
estpost tabstat y_hat_11, by(in_control) c(stat) stat(mean, sd, min, median, max, count)
esttab, cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) p50(fmt(%13.2fc)) max(fmt(%13.2fc))  count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Median" "Max" "N")
esttab using "table_11b.tex", replace cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) p50(fmt(%13.2fc)) max(fmt(%13.2fc))  count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Median" "Max" "N")

********************************************************************************
* problem #12 - linear regression of in_control out-of-sample
********************************************************************************

regress re78 age age_2 educ black hisp married nodegree re74 re75 if in_control == 0, robust 
outreg2 using table_12a, tex(frag) replace

predict y_hat_12

est clear
estpost tabstat y_hat_12, by(in_control) c(stat) stat(mean, sd, min, median, max, count)
esttab, cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) p50(fmt(%13.2fc)) max(fmt(%13.2fc))  count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Median" "Max" "N")
esttab using "table_12b.tex", replace cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) p50(fmt(%13.2fc)) max(fmt(%13.2fc))  count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Median" "Max" "N")

********************************************************************************
* problem #13 - inverse probability weighting
********************************************************************************

* get number of treated
count if in_control == 1
scalar n_1 = r(N)

* get number of untreated
count if in_control == 0
scalar n_0 = r(N)

* get unconditional probability of treatment
scalar p_hat = n_1/(n_0 + n_1)

* get y * d
gen y_d = re78 * in_control
summarize y_d
scalar y_d_sum = r(sum)

* get second term wo rescaling
gen wo_rescale = (1 - p_hat) / p_hat * pscoreb * re78 * (1 - in_control)/(1 - pscoreb)
summarize wo_rescale
scalar wo_rescale_sum = r(sum)

* get second term w rescaling
gen rescaling = pscoreb * (1-in_control)/(1-pscoreb) 
summarize rescaling
scalar rescaling = 1/n_0 * r(sum)
gen w_rescale = (1/rescaling) * pscoreb * re78 * (1 - in_control)/(1-pscoreb)
summarize w_rescale
scalar w_rescale_sum = r(sum)

* treatment effect on treated ipw w and wo rescaling
display 1/n_1 * y_d_sum - 1/n_0 * wo_rescale_sum
display 1/n_1 * y_d_sum - 1/n_0 * w_rescale_sum








