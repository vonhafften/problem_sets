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
outreg2 using table_3, tex(frag) replace addstat(Failures completely determined, e(N_cdf), Successes completely determined,  e(N_cds))
predict pscorea, pr

probit in_control age age_2 educ black hisp married nodegree re74 re75
outreg2 using table_3, tex(frag) append addstat(Failures completely determined, e(N_cdf), Successes completely determined,  e(N_cds))
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

texdoc init "table_7.tex", replace

/*tex
\begin{table}[h!]
	\label{table_7}
\begin{center}
\begin{tabular}{rll}
\toprule
Test & Test & Test \\
\hline
tex*/

psmatch2 in_control, noreplacement outcome(re78) pscore(pscorea) neighbor(1) common

list if _support == 0

psmatch2 in_control, noreplacement outcome(re78) pscore(pscoreb) neighbor(1) common
list if _support == 0

* figure out outreg-ing

********************************************************************************
* problem #7 - Nearest neighbor w replacement
********************************************************************************

texdoc init "table_7.tex", replace


psmatch2 in_control, outcome(re78) pscore(pscorea) neighbor(1) common
list if _support == 0

psmatch2 in_control, outcome(re78) pscore(pscoreb) neighbor(1) common
list if _support == 0

