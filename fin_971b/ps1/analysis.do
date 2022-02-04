* Corporate Finance FIN 971B
* Taught by Oliver Levine
* Problem set 1
* Replication of Lemmon, Roberts, and Zender (2018)
* By Alex von Hafften

cd "/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_971b/ps1"

* Install user-defined functions
ssc install asrol
ssc install winsor2
ssc install reghdfe
ssc install ftools
ssc install extremes 
ssc install outreg2
ssc install estout

*******************************************************************
* computes equity volatility from CRPS data
*******************************************************************

clear

* crsp.dta is from the "LRZ (2018) Extension" data query on WRDS.
use crsp
gen monthly_date = mofd(date)

xtset PERMNO monthly_date, monthly

drop if PRC < 0
drop if missing(PRC)

* computes stock return
gen s_return = (PRC - L.PRC)/L.PRC

* computes rolling standard deviation of stock return
bys PERMNO: asrol s_return, stat(sd) window(date 36)

* saves dataset to be merged with compustat
save crsp_2, replace

*******************************************************************
* Main compustat datawork starts here
*******************************************************************

* compustat.dta is from the "Lemmon, Roberts, and Zender (2018) Replication" data query on WRDS.

clear

use compustat
gen PERMNO = LPERMNO
gen date = datadate

keep if indfmt=="INDL" & datafmt=="STD" & popsrc=="D" & consol=="C"

* This drops 14 observations that have duplicates firm id (LPERMNO fyear) and time id (fyear)
gen id = LPERMNO
sort id fyear
quietly by id fyear:  gen dup = cond(_N==1,0,_n)
drop if dup

* tells stata it's a panel data set
xtset id fyear 

* merge on the stock return volatility
merge 1:1 PERMNO date using crsp_2, keepusing(s_return_sd36)

*******************************************************************
* Create variables 
*******************************************************************

* based on the data appendix
gen debt = dlc + dltt
gen b_lev = debt/at
gen size = log(at)
gen profit = oibdp/at
gen m_equity = prcc_f * cshpri
gen m_lev = debt / (debt + m_equity)
gen m_b = (m_equity + debt + pstkl - txditc) / at
gen collat = (invt + ppent) / at
gen z_score = (3.3*pi + sale + 1.4 * re + 1.2 *(act - lct) )/at
gen tang = ppent/at
gen l_sales = log(sale)

* others
gen intangibles = intan/at
gen pays_dv = dv > 0

* creates survivor dummy
by id, sort: egen n_periods = count(_n)
gen survivor = n_periods >= 20

* creates initial leverage
bys id: gen i_b_lev = b_lev[1]
bys id: gen i_m_lev = m_lev[1]

* computes cash flow volatility - rolling sd of cashflow in past 3 years
bys id: asrol oibdp, stat(sd) window(fyear 3) min(3)
gen cf_vol = oibdp_sd3 / at

* create roa volatility from crsp data
gen roa_vol = m_equity / (m_equity + debt) * s_return_sd36

* trims variables between p1 and p99
winsor2 b_lev, trim
winsor2 m_lev, trim
winsor2 l_sales, trim
winsor2 m_b, trim
winsor2 profit, trim
winsor2 tang, trim
winsor2 cf_vol, trim
winsor2 intangibles, trim

* normalizes (i.e. substract mean and divide by std dev)  for each regression analysis
* coefficients can be interpretted as a one standard deviation change in the independent var

norm i_b_lev, method(zee)
norm i_m_lev, method(zee)
norm l_sales_tr, method(zee)
norm m_b_tr, method(zee)
norm profit_tr, method(zee)
norm tang_tr, method(zee)
norm cf_vol_tr, method(zee)
norm intangibles_tr, method(zee)
norm roa_vol, method(zee)

*******************************************************************
* filter based on the beginning of "Part 1. Data and Sample Selection"
*******************************************************************

drop if missing(at)
drop if m_lev > 1
drop if m_lev < 0
drop if b_lev > 1
drop if b_lev < 0

* drops the outliers
gen outlier = b_lev_tr+m_lev_tr+l_sales_tr+m_b_tr+profit_tr+tang_tr+cf_vol_tr+intangibles_tr
drop if missing(outlier)

* creates median book leverage by industry
bys gind fyear, sort: egen ind_b_lev = median(b_lev_tr)
norm ind_b_lev, method(zee)

*******************************************************************
* Replication of table 1
*******************************************************************

* Makes variable names pretty for the table

gen Book_Leverage = b_lev_tr
gen Market_Leverage = m_lev_tr
gen log_Sales = l_sales_tr
gen Market_to_Book = m_b_tr
gen Profitability = profit_tr
gen Tangibility = tang_tr
gen Cash_Flow_Volatility = cf_vol_tr
gen Median_Industry_Book_Leverage = ind_b_lev
gen Dividend_Payer = pays_dv
gen Intangible_Assets = intangibles_tr

est clear
estpost tabstat Book_Leverage Market_Leverage log_Sales Market_to_Book Profitability Tangibility Cash_Flow_Volatility Median_Industry_Book_Leverage Dividend_Payer Intangible_Assets, c(stat) stat(mean, median, sd, count)
esttab, cells("mean(fmt(%13.2fc)) p50(fmt(%13.2fc)) sd(fmt(%13.2fc)) count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "Median" "SD" "N")
esttab using "table_1_a.tex", replace cells("mean(fmt(%13.2fc)) p50(fmt(%13.2fc)) sd(fmt(%13.2fc)) count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "Median" "SD" "N")
  
est clear
estpost tabstat Book_Leverage Market_Leverage log_Sales Market_to_Book Profitability Tangibility Cash_Flow_Volatility Median_Industry_Book_Leverage Dividend_Payer Intangible_Assets, c(stat) stat(mean, median, sd, count), if survivor
esttab, cells("mean(fmt(%13.2fc)) p50(fmt(%13.2fc)) sd(fmt(%13.2fc)) count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "Median" "SD" "N")
esttab using "table_1_b.tex", replace cells("mean(fmt(%13.2fc)) p50(fmt(%13.2fc)) sd(fmt(%13.2fc)) count(fmt(%13.0fc))") nonumber nomtitle nonote noobs label collabels("Mean" "Median" "SD" "N")

*******************************************************************
* Replication of table 2
*******************************************************************

* Makes variable names pretty for the table

gen Initial_Book_Leverage = zee_i_b_lev
replace log_Sales = zee_l_sales_tr
replace Market_to_Book = zee_m_b_tr
replace Profitability = zee_profit_tr
replace Tangibility = zee_tang_tr
replace Cash_Flow_Volatility = zee_cf_vol_tr
replace Median_Industry_Book_Leverage = zee_ind_b_lev
gen ROA_Volatility = zee_roa_vol

* Panel A
* All firm. Dependent variable book leverage
reg  Book_Leverage Initial_Book_Leverage, robust cluster(id)
outreg2 using table_2_a, tex(frag) replace nocons

areg Book_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility, absorb(fyear) robust cluster(id)
outreg2 using table_2_a, tex(frag) append nocons

areg Book_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage Cash_Flow_Volatility Dividend_Payer, absorb(fyear) robust cluster(id)
outreg2 using table_2_a, tex(frag) append nocons

* All firm. Dependent variable market leverage
reg  Market_Leverage Initial_Book_Leverage, robust cluster(id)
outreg2 using table_2_a, tex(frag) append nocons

areg Market_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility, absorb(fyear) robust cluster(id)
outreg2 using table_2_a, tex(frag) append nocons

areg Market_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage Cash_Flow_Volatility Dividend_Payer, absorb(fyear) robust cluster(id)
outreg2 using table_2_a, tex(frag) append nocons

* Panel B
* Survivors. Dependent variable book leverage
reg  Book_Leverage Initial_Book_Leverage, robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) replace nocons

areg Book_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility, absorb(fyear) robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append nocons

areg Book_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage Cash_Flow_Volatility Dividend_Payer, absorb(fyear) robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append nocons


* Survivors. Dependent variable market leverage
reg  Market_Leverage Initial_Book_Leverage, robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append nocons

areg Market_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility, absorb(fyear) robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append nocons

areg Market_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage Cash_Flow_Volatility Dividend_Payer, absorb(fyear) robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append nocons

*******************************************************************
* Replication of table 2 Panel A but with firm FEs instead of initial leverage
*******************************************************************

* Panel A
* All firm. Dependent variable book leverage
areg  Book_Leverage, absorb(id) cluster(id)
outreg2 using table_2_a_fe, tex(frag) replace nocons

reghdfe Book_Leverage log_Sales Market_to_Book Profitability Tangibility, absorb(fyear id) vce(cluster id)
outreg2 using table_2_a_fe, tex(frag) append nocons

reghdfe Book_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage Cash_Flow_Volatility Dividend_Payer, absorb(fyear id) vce(cluster id)
outreg2 using table_2_a_fe, tex(frag) append nocons

* All firm. Dependent variable market leverage
areg  Market_Leverage, absorb(id) cluster(id)
outreg2 using table_2_a_fe, tex(frag) append nocons

reghdfe Market_Leverage log_Sales Market_to_Book Profitability Tangibility, absorb(fyear id) vce(cluster id)
outreg2 using table_2_a_fe, tex(frag) append nocons

reghdfe Market_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage Cash_Flow_Volatility Dividend_Payer, absorb(fyear id) vce(cluster id)
outreg2 using table_2_a_fe, tex(frag) append nocons

*******************************************************************
* Replication of table 2 Panel A with return on assets from crsp 
*******************************************************************

areg Book_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage Cash_Flow_Volatility Dividend_Payer, absorb(fyear) robust cluster(id)
outreg2 using table_2_a_crsp, tex(frag) replace nocons

areg Book_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage ROA_Volatility Dividend_Payer, absorb(fyear) robust cluster(id)
outreg2 using table_2_a_crsp, tex(frag) append nocons

areg Market_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage Cash_Flow_Volatility Dividend_Payer, absorb(fyear) robust cluster(id)
outreg2 using table_2_a_crsp, tex(frag) append nocons

areg Market_Leverage Initial_Book_Leverage log_Sales Market_to_Book Profitability Tangibility Median_Industry_Book_Leverage ROA_Volatility Dividend_Payer, absorb(fyear) robust cluster(id)
outreg2 using table_2_a_crsp, tex(frag) append nocons

*******************************************************************
* Comparison of market equity from compustat vs crsp
*******************************************************************

clear

use compustat

gen PERMNO = LPERMNO
gen date = datadate

* in compustat: LPERMNO and datadate
* in crsp: PERMNO and date

merge 1:1 PERMNO date using crsp, keepusing(PRC SHROUT)

drop if PRC < 0
drop if missing(PRC)
drop if LINKDT > datadate
drop if LINKENDDT < datadate
keep if LINKTYPE == "LU" | LINKTYPE == "LC"
keep if LINKPRIM == "C"

gen Market_Equity_Compustat = prcc_f * cshpri
gen Market_Equity_CRPS = PRC * SHROUT / 1000

drop if missing(Market_Equity_Compustat)

reg Market_Equity_Compustat Market_Equity_CRPS
outreg2 using m_equity_reg, tex(frag) replace

scatter Market_Equity_Compustat Market_Equity_CRPS
graph export m_equity_scatter.png, width(600) height(450) replace

predict residual_m_equity, residuals

extremes residual_m_equity conm fyear

