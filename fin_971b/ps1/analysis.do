* Corporate Finance FIN 971B
* Taught by Oliver Levine
* Problem set 1
* Replication of Lemmon, Roberts, and Zender (2018)
* By Alex von Hafften

cd "/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_971b/ps1"

ssc install asrol
ssc install winsor2
ssc install reghdfe
ssc install ftools
ssc install extremes 
ssc install outreg2

* crsp.dta is from the "LRZ (2018) Extension" data query on WRDS.
* compute equity volatility from CRPS

clear

use crsp
gen monthly_date = mofd(date)

xtset PERMNO monthly_date, monthly

drop if PRC < 0
drop if missing(PRC)

gen s_return = (PRC - L.PRC)/L.PRC
bys PERMNO: asrol s_return, stat(sd) window(date 36)

save crsp_2, replace

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

* creates cash flow volatility
bys id: asrol oibdp, stat(sd) window(fyear 3) min(3)
gen cf_vol = oibdp_sd3 / at

* create roa volatility from crsp data
gen roa_vol = m_equity / (m_equity + debt) * s_return_sd36

* trims variables
winsor2 b_lev, trim
winsor2 m_lev, trim
winsor2 l_sales, trim
winsor2 m_b, trim
winsor2 profit, trim
winsor2 tang, trim
winsor2 cf_vol, trim
winsor2 intangibles, trim

* normalize (i.e. substract mean and divide by std dev for each regression var)
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

tabstat b_lev_tr m_lev_tr l_sales_tr m_b_tr profit_tr tang_tr cf_vol_tr ind_b_lev pays_dv intangibles_tr, stat(mean, median, sd, count)

tabstat b_lev_tr m_lev_tr l_sales_tr m_b_tr profit_tr tang_tr cf_vol_tr ind_b_lev pays_dv intangibles_tr, stat(mean, median, sd, count), if survivor

* To do
* - Export summary statistics tables

*******************************************************************
* Replication of table 2
*******************************************************************

* Panel A
* All firm. Dependent variable book leverage
reg  b_lev_tr zee_i_b_lev, robust cluster(id)
outreg2 using table_2_a, tex(frag) replace

areg b_lev_tr zee_i_b_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr, absorb(fyear) robust cluster(id)
outreg2 using table_2_a, tex(frag) append

areg b_lev_tr zee_i_b_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr zee_ind_b_lev zee_cf_vol_tr pays_dv, absorb(fyear) robust cluster(id)
outreg2 using table_2_a, tex(frag) append

* All firm. Dependent variable market leverage
reg  m_lev_tr zee_i_m_lev, robust cluster(id)
outreg2 using table_2_a, tex(frag) append

areg m_lev_tr zee_i_m_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr, absorb(fyear) robust cluster(id)
outreg2 using table_2_a, tex(frag) append

areg m_lev_tr zee_i_m_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr zee_ind_b_lev zee_cf_vol_tr pays_dv, absorb(fyear) robust cluster(id)
outreg2 using table_2_a, tex(frag) append

* Panel B
* Survivors. Dependent variable book leverage
reg  b_lev_tr zee_i_b_lev, robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) replace

areg b_lev_tr zee_i_b_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr, absorb(fyear) robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append

areg b_lev_tr zee_i_b_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr zee_ind_b_lev zee_cf_vol_tr pays_dv, absorb(fyear) robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append


* Survivors. Dependent variable market leverage
reg  m_lev_tr zee_i_m_lev, robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append

areg m_lev_tr zee_i_m_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr, absorb(fyear) robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append

areg m_lev_tr zee_i_m_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr zee_ind_b_lev zee_cf_vol_tr pays_dv, absorb(fyear) robust cluster(id), if survivor
outreg2 using table_2_b, tex(frag) append

*******************************************************************
* Replication of table 2 Panel A but with firm FEs instead of initial leverage
*******************************************************************

* Panel A
* All firm. Dependent variable book leverage
areg  b_lev_tr, absorb(id) robust cluster(id)
outreg2 using table_2_a_fe, tex(frag) replace

reghdfe b_lev_tr zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr, absorb(fyear id) vce(cluster id)
outreg2 using table_2_a_fe, tex(frag) append

reghdfe b_lev_tr zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr zee_ind_b_lev zee_cf_vol_tr pays_dv, absorb(fyear id) vce(cluster id)
outreg2 using table_2_a_fe, tex(frag) append

* All firm. Dependent variable market leverage
areg  m_lev_tr, absorb(id) robust cluster(id)
outreg2 using table_2_a_fe, tex(frag) append

reghdfe m_lev_tr zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr, absorb(fyear id) vce(cluster id)
outreg2 using table_2_a_fe, tex(frag) append

reghdfe m_lev_tr zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr zee_ind_b_lev zee_cf_vol_tr pays_dv, absorb(fyear id) vce(cluster id)
outreg2 using table_2_a_fe, tex(frag) append

*******************************************************************
* Replication of table 2 Panel A with return on assets from crsp 
*******************************************************************

areg b_lev_tr zee_i_b_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr zee_ind_b_lev roa_vol pays_dv, absorb(fyear) robust cluster(id)
outreg2 using table_2_a_crsp, tex(frag) replace

areg m_lev_tr zee_i_m_lev zee_l_sales_tr zee_m_b_tr zee_profit_tr zee_tang_tr zee_ind_b_lev roa_vol pays_dv, absorb(fyear) robust cluster(id)
outreg2 using table_2_a_crsp, tex(frag) append

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

gen m_equity_compustat = prcc_f * cshpri
gen m_equity_crsp = PRC * SHROUT / 1000

drop if missing(m_equity_compustat)

reg m_equity_compustat m_equity_crsp
outreg2 using m_equity_reg, tex(frag) replace

scatter m_equity_compustat m_equity_crsp
graph export m_equity_scatter.png, width(600) height(450)

predict residual_m_equity, residuals

extremes residual_m_equity conm fyear

* To do 
* - outlier tables

