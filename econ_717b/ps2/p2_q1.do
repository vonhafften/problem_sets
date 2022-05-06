* Applied Econometrics
* ECON 717B
* Taught by Matt Wiswall
* Problem set 2
* By Alex von Hafften

cd "/Users/vonhafften/Documents/UW Madison/problem_sets/econ_717b/ps2"

* Install user-defined functions
ssc install outreg2

**********************************************
* part a - generate data
**********************************************

clear all

set seed 456

* unobservables
drawnorm epsilon, n(2000) sd(0.5)
drawnorm a,   n(2000) sd(4)
drawnorm eta, n(2000) sd(1)

* observables
drawnorm z_1, n(2000) sd(0.1)
drawnorm z_2, n(2000) sd(25)
generate z_3  = runiform()

* schooling and wages
generate s = 3 * a + z_1 + z_2 + eta
generate ln_w = 1 + 0.05 * s + 0.1 * a + epsilon

**********************************************
* part b - ols
**********************************************

regress ln_w s
outreg2 using table_1_p2_q2, tex(frag) replace ctitle(OLS)

**********************************************
* part b - 2sls
**********************************************

ivregress 2sls ln_w (s = z_1)
estat firststage
outreg2 using table_1_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_1)

ivregress 2sls ln_w (s = z_2)
estat firststage
outreg2 using table_1_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_2)

ivregress 2sls ln_w (s = z_3)
estat firststage
outreg2 using table_1_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_3)

ivregress 2sls ln_w (s = z_1 z_2)
estat firststage
outreg2 using table_1_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_1 z_2)

ivregress 2sls ln_w (s = z_2 z_3)
estat firststage
outreg2 using table_1_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig))addtext(Instruments, z_2 z_3)

ivregress 2sls ln_w (s = z_1 z_3)
estat firststage
outreg2 using table_1_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_1 z_3)

ivregress 2sls ln_w (s = z_1 z_2 z_3)
estat firststage
outreg2 using table_1_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_1 z_2 z_3)

**********************************************
* part e - see write-up
**********************************************

**********************************************
* part f - see write-up
**********************************************

clear all

set seed 456

* unobservables
drawnorm epsilon, n(500000) sd(0.5)
drawnorm a,   n(500000) sd(4)
drawnorm eta, n(500000) sd(1)

* observables
drawnorm z_1, n(500000) sd(0.1)
drawnorm z_2, n(500000) sd(25)
generate z_3  = runiform()

* schooling and wages
generate s = 3 * a + z_1 + z_2 + eta
generate ln_w = 1 + 0.05 * s + 0.1 * a + epsilon

* ols
regress ln_w s
outreg2 using table_2_p2_q2, tex(frag) replace ctitle(OLS)

* 2sls
ivregress 2sls ln_w (s = z_1)
estat firststage
outreg2 using table_2_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_1)

ivregress 2sls ln_w (s = z_2)
estat firststage
outreg2 using table_2_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_2)

ivregress 2sls ln_w (s = z_3)
estat firststage
outreg2 using table_2_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_3)

ivregress 2sls ln_w (s = z_1 z_2)
estat firststage
outreg2 using table_2_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_1 z_2)

ivregress 2sls ln_w (s = z_2 z_3)
estat firststage
outreg2 using table_2_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_2 z_3)

ivregress 2sls ln_w (s = z_1 z_3)
estat firststage
outreg2 using table_2_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_1 z_3)

ivregress 2sls ln_w (s = z_1 z_2 z_3)
estat firststage
outreg2 using table_2_p2_q1, tex(frag) append ctitle(2SLS) addstat(F-Statistic, r(mineig)) addtext(Instruments, z_1 z_2 z_3)
