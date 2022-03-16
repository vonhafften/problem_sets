log using analysis, replace

* ECON 717A: Applied Econometrics
* Problem Set 3
* Professor: Jeff Smith
* Alex von Hafften
* Diff-in-diff

* clear workspace
clear

* install user defined functions (if needed)
ssc install outreg2
ssc install asreg


* change working directory
cd "/Users/vonhafften/Documents/UW Madison/problem_sets/econ_717a/ps3/"

* open dataset
use "Economics 717 Miron and Tetelbaum Data"

********************************************************************************
* problem #1 - create panel
********************************************************************************

xtset state year

********************************************************************************
* problem #2 - create treatment indicator
********************************************************************************

gen mlda21 = (mlda == 21)
tab mlda21

********************************************************************************
* problem #3 - naive treatment estimate
********************************************************************************

regress rate18_20ht mlda21, robust
outreg2 using table_3, addtext(Fixed Effects, None, Clusters, None) tex(frag) replace

* summary statistics by year
tab year, sum(rate18_20ht)
tab year, sum( mlda21)

* run ar regression by state
gen rate18_20ht_ar = .
gen mlda21_ar = .

sort state year 
by state: gen rate18_20ht_lag = rate18_20ht[_n-1]
by state: gen mlda21_lag = mlda21[_n-1]

egen group = group(state)
summ group

forvalues i = 1/`r(max)' {
	reg rate18_20ht rate18_20ht_lag if `i'==group
	replace rate18_20ht_ar = _b[rate18_20ht_lag] if `i'==group
	
	reg mlda21 mlda21_lag if `i'==group
	replace mlda21_ar = _b[mlda21_lag] if `i'==group
}

* if mlda never changes, then set ar coefficient to 1.
replace mlda21_ar = 1 if mldayr == 0

summ rate18_20ht_ar
summ mlda21_ar

********************************************************************************
* problem #4 - naive treatment estimate with state or year fes
********************************************************************************

regress rate18_20ht mlda21 i.state, robust
outreg2 using table_4_state, addtext(Fixed Effects, State, Clusters, None) tex(frag) replace

regress rate18_20ht mlda21 i.year, robust
outreg2 using table_4_year, addtext(Fixed Effects, Year, Clusters, None) tex(frag) replace

********************************************************************************
* problem #5 - naive treatment estimate 
* with state and year fes and clustered ses
********************************************************************************

regress rate18_20ht mlda21 i.year i.state, robust cluster(state)
outreg2 using table_5, keep(mlda21) addtext(Fixed Effects, State and Year, Clusters, States) nocons tex(frag) replace

********************************************************************************
* problem #6 - treatment estimate 
* with state and year fes and unclustered ses
********************************************************************************

regress rate18_20ht mlda21 i.year i.state, robust
outreg2 using table_6, keep(mlda21) addtext(Fixed Effects, State and Year, Clusters, None) nocons tex(frag) replace

********************************************************************************
* problem #7 - treatment estimate 
* with state and year fes and unclustered ses
* before and including 1990
********************************************************************************

regress rate18_20ht mlda21 i.year i.state if year <= 1990, robust cluster(state) 
outreg2 using table_7, keep(mlda21)  addtext(Fixed Effects, State and Year, Clusters, States, Period, Before 1990)  nocons tex(frag) replace

********************************************************************************
* problem #8 - placebo test
********************************************************************************

* placebo treatment indicator equals one in states treated in 1987 and in years 1982 or later
gen placebo82 = (mldayr == 1987) & (year >= 1982)

* in_placebo82_sample is 1 for states with 21 drinking ages always (mldayr == 0) and that switched in 1987 (mldayr == 1987)
gen in_placebo82_sample = ((mldayr == 0) | (mldayr == 1987)) & (year < 1987)

regress rate18_20ht placebo82 i.year i.state if in_placebo82_sample == 1, robust cluster(state)
outreg2 using table_8, keep(placebo82) addtext(Fixed Effects, State and Year, Clusters, States) tex(frag) nocons replace

********************************************************************************
* problem #9 - mi and md
********************************************************************************

gen in_mi_sample = (mldayr == 0) | (state == 23)
gen mlda21_mi = (mlda21 == 1) & (state == 23)
regress rate18_20ht mlda21_mi i.year i.state if in_mi_sample, robust cluster(state) 
outreg2 using table_9, keep(mlda21_mi mlda21_md)  ctitle(Michigan) addtext(Fixed Effects, State and Year)  tex(frag) nocons replace

gen in_md_sample = (mldayr == 0) | (state == 21)
gen mlda21_md = (mlda21 == 1) & (state == 21)
regress rate18_20ht mlda21_md i.year i.state if in_md_sample, robust cluster(state) 
outreg2 using table_9, keep(mlda21_mi mlda21_md)  ctitle(Maryland) addtext(Fixed Effects, State and Year)  tex(frag) nocons append

********************************************************************************
* problem #10 - early vs late treatment
********************************************************************************

gen mlda21_14 = mlda21 & (year < mldayr + 4)
gen mlda_later = mlda21 & (year >= mldayr + 4)

regress rate18_20ht mlda21_14 mlda_later i.year i.state, robust cluster(state) 
outreg2 using table_10, keep(mlda21_14 mlda_later)  addtext(Fixed Effects, State and Year)  tex(frag) nocons replace


********************************************************************************
* create and translate log file
********************************************************************************

log close

translate analysis.smcl analysis.log, linesize(255) replace





