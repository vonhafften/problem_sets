* Corporate Finance FIN 971B
* Taught by Oliver Levine
* Problem set 2
* Q-Theory Regressions
* By Alex von Hafften

cd "/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_971b/ps2"

* Install user-defined functions
ssc install asrol
ssc install winsor2
ssc install reghdfe
ssc install ftools
ssc install extremes 
ssc install outreg2
ssc install estout
ssc install egenmore
ssc install xtewreg


*******************************************************************
* prepare PCE data
*******************************************************************

clear

import delimited PCEPILFE

* sets 2004 december the base year 
gen inflation_adjustment = pcepilfe / 87.311 * 100
gen monthly_date = mofd(date(date, "YMD"))

save inflation, replace

*******************************************************************
* Main compustat datawork starts here
*******************************************************************

* compustat.dta is from the "FHP (1988) Replication" data query on WRDS.

clear

use compustat
gen id = real(GVKEY)
gen date = datadate
gen sic_n = real(sic)

* filtering
keep if indfmt=="INDL" & datafmt=="STD" & popsrc=="D" & consol=="C"

* positive book assts
keep if at > 0 

* no financials or utilities
keep if (sic_n < 6000 | sic_n > 6999)
keep if (sic_n < 4900 | sic_n > 4999)

* drop if missing id or year
keep if GVKEY != "" & fyear != .

* drop if non us
keep if fic == "USA"

* drop if missing stock data
keep if prcc_f!=. &csho!=.

* This drops 3,857 observations that have duplicates firm id (GVKEY) and time id (fyear)
sort id fyear
quietly by id fyear:  gen dup = cond(_N==1,0,_n)
drop if dup

* tells stata it's a panel data set
xtset id fyear 
gen monthly_date = mofd(date)

* merge on inflation adjustment
merge m:1 monthly_date using inflation, keepusing(inflation_adjustment)

*******************************************************************
* Create variables 
*******************************************************************

* based on the EW 2012 appendix
gen investment = capx/at
gen cash_flow = (ib + dpc) / at
gen tobin_q = (at + prcc_f * csho - ceq - txdb)/ppegt


* trims variables between p1 and p99
winsor2 investment, trim
winsor2 cash_flow, trim
winsor2 tobin_q, trim

* create lags
gen investment_tr_l = investment_tr[_n-1]
gen cash_flow_tr_l = cash_flow_tr[_n-1]
gen tobin_q_tr_l = tobin_q_tr[_n-1]

* bin based on leverage
gen book_leverage = (dlc + dltt) / at
egen leverage_bin = xtile(book_leverage), by (fyear) nq(5)

* sa index from HP 2010
gen at_2004 = at * inflation_adjustment/100
gen size = min(log(at_2004), log(4500))
sort GVKEY fyear 
by GVKEY: gen age = min(_n, 37)
gen sa_index = (-0.737*size) + (0.043 * size^2) - (0.040*age)
egen sa_index_bin = xtile(sa_index), nq(5)

*******************************************************************
* Summary statistics
*******************************************************************

est clear
estpost tabstat investment_tr tobin_q_tr cash_flow_tr, c(stat) stat(mean, median, sd, count)
esttab, cells("mean(fmt(%13.2fc)) p50(fmt(%13.2fc)) sd(fmt(%13.2fc)) count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "Median" "SD" "N")
esttab using "table_0.tex", replace cells("mean(fmt(%13.2fc)) p50(fmt(%13.2fc)) sd(fmt(%13.2fc)) count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "Median" "SD" "N")

*******************************************************************
* Problem 1
*******************************************************************

* without cash flow or firm fes
reg investment_tr tobin_q_tr_l
outreg2 using table_1, tex(frag) replace nocons addtext(Firm FEs, No)


* without cash flow. with firm fes
areg investment_tr tobin_q_tr_l, absorb(id)
outreg2 using table_1, tex(frag) append nocons addtext(Firm FEs, Yes)


* with cash flow. without firm fes
reg investment_tr tobin_q_tr_l cash_flow_tr
outreg2 using table_1, tex(frag) append nocons addtext(Firm FEs, No)


* with cash flow. with firm fes
areg investment_tr tobin_q_tr_l cash_flow_tr, absorb(id)
outreg2 using table_1, tex(frag) append nocons addtext(Firm FEs, Yes)

*******************************************************************
* Problem 2
*******************************************************************

* with cash flow. without firm fes
reg investment_tr tobin_q_tr_l cash_flow_tr_l
outreg2 using table_2, tex(frag) replace nocons  addtext(Firm FEs, No)


* with cash flow. with firm fes
areg investment_tr tobin_q_tr_l cash_flow_tr_l, absorb(id) 
outreg2 using table_2, tex(frag) append nocons addtext(Firm FEs, Yes)

*******************************************************************
* Problem 3
*******************************************************************

gen tobin_q_tr_l_2 = tobin_q_tr_l^2

* without cash flow or firm fes
reg investment_tr tobin_q_tr_l tobin_q_tr_l_2
outreg2 using table_3, tex(frag) replace nocons addtext(Firm FEs, No)


* without cash flow. with firm fes
areg investment_tr tobin_q_tr_l tobin_q_tr_l_2, absorb(id)
outreg2 using table_3, tex(frag) append nocons addtext(Firm FEs, Yes)

*******************************************************************
* Problem 4
*******************************************************************

* without cash flow or firm fes
reg investment_tr tobin_q_tr_l, robust
outreg2 using table_4, tex(frag) replace nocons addtext(Firm FEs, No)


* without cash flow. with firm fes
areg investment_tr tobin_q_tr_l, absorb(id) robust 
outreg2 using table_4, tex(frag) append nocons addtext(Firm FEs, Yes)


* with cash flow. without firm fes
reg investment_tr tobin_q_tr_l cash_flow_tr, robust
outreg2 using table_4, tex(frag) append nocons addtext(Firm FEs, No)


* with cash flow. with firm fes
areg investment_tr tobin_q_tr_l cash_flow_tr, absorb(id) robust
outreg2 using table_4, tex(frag) append nocons addtext(Firm FEs, Yes)

*******************************************************************
* Problem 5
*******************************************************************

* see calculation in latex code

*******************************************************************
* Problem 6
*******************************************************************

areg investment_tr tobin_q_tr_l cash_flow_tr if leverage_bin == 1, absorb(id) robust 
outreg2 using table_6, tex(frag) replace nocons addtext(Firm FEs, Yes)


areg investment_tr tobin_q_tr_l cash_flow_tr if leverage_bin == 2, absorb(id) robust 
outreg2 using table_6, tex(frag) append nocons addtext(Firm FEs, Yes)


areg investment_tr tobin_q_tr_l cash_flow_tr if leverage_bin == 3, absorb(id) robust 
outreg2 using table_6, tex(frag) append nocons addtext(Firm FEs, Yes)

areg investment_tr tobin_q_tr_l cash_flow_tr if leverage_bin == 4, absorb(id) robust 
outreg2 using table_6, tex(frag) append nocons addtext(Firm FEs, Yes)

areg investment_tr tobin_q_tr_l cash_flow_tr if leverage_bin == 5, absorb(id) robust 
outreg2 using table_6, tex(frag) append nocons addtext(Firm FEs, Yes)

*******************************************************************
* Problem 7
*******************************************************************

areg investment_tr tobin_q_tr_l cash_flow_tr if sa_index_bin == 1, absorb(id) robust 
outreg2 using table_7, tex(frag) replace nocons addtext(Firm FEs, Yes)

areg investment_tr tobin_q_tr_l cash_flow_tr if sa_index_bin == 2, absorb(id) robust 
outreg2 using table_7, tex(frag) append nocons addtext(Firm FEs, Yes)

areg investment_tr tobin_q_tr_l cash_flow_tr if sa_index_bin == 3, absorb(id) robust 
outreg2 using table_7, tex(frag) append nocons addtext(Firm FEs, Yes)

areg investment_tr tobin_q_tr_l cash_flow_tr if sa_index_bin == 4, absorb(id) robust 
outreg2 using table_7, tex(frag) append nocons addtext(Firm FEs, Yes)

areg investment_tr tobin_q_tr_l cash_flow_tr if sa_index_bin == 5, absorb(id) robust 
outreg2 using table_7, tex(frag) append nocons addtext(Firm FEs, Yes)

*******************************************************************
* Problem 8
*******************************************************************

xtewreg investment_tr tobin_q_tr_l , maxdeg(5)
outreg2 using table_8, tex(frag) replace nocons addtext(Firm FEs, No)

xtewreg investment_tr tobin_q_tr_l cash_flow_tr , maxdeg(5) mis(2) cent(set)
outreg2 using table_8, tex(frag) append nocons addtext(Firm FEs, No)


