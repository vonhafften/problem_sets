* ECON 717A Final Project
* Alex von Hafften

* clear workspace
clear

* install user defined functions (if needed)
ssc install outreg2

cd "/Users/vonhafften/Documents/UW Madison/problem_sets/econ_717a/project/"

********************************************************************************
* commercial bank data
********************************************************************************

use "jlnb4ujafkjqwza4"

gen BHCK2170 = RCON2170
gen BHCK3210 = RCON3210
replace RSSD9001 = 1027004
drop if RSSD9999 < 20180930

save "commercial_bank.dta"

clear

********************************************************************************
* BHC data
********************************************************************************

use "ct2lkiyvusf2codh"

append using "commercial_bank.dta"

gen leverage = (BHCK2170 - BHCK3210)/BHCK2170

clear
