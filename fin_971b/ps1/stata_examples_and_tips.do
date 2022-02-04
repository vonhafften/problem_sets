/* This file summarizes useful Stata coding tips */
* Author: Helene Hall
* Date: June 19, 2020

/*----- Inputting called variables -----*/

* you can call a stata .do file within another stata file. Suppose there exists another file called
*	subfile.do. You can input values to the stata file like so: 
		do subfile.do input1 input2 

		* then, in subfile.do, at the top of the file you would define these inputs as locals and use them in the file's evaluation with:
		local input1 `1'
		local input2 `2'



/*----- File Header Tips -----*/
clear all 
set more off 

* ssc install a program 
capture ssc install texsave
capture ssc install freduse
capture ssc install labmask 

* Define paths
global folder "path/to/folder/code"
cd $folder 

global input  	"../input"
global output	"../output"
global temp		"../temp"
global code		"../code"
* generate a subfolder in output corresponding to today's date
local today : display %tdCY-N-D date("`c(current_date)'", "DMY")
global output_today "$output/`today'"
capture mkdir $output_today

* Define other globals
global restriction "no"


/* create programs */

	* define program to export basic summary table
	capture program drop output_sum_stats
	program define output_sum_stats
		args var_for_summary output_path
		eststo: estpost tabstat `var_for_summary', statistics(mean sd min p1 p10 p25 p50 p75 p90 p99 max)
	    esttab using `output_path', replace nodepvars noobs nonumbers ///
	        cell("mean (fmt(%9.3f)) sd min p1 p10 p25 p50 p75 p90 p99 max") mlabels(none) 
	    eststo clear 
	end

* Define a global using a conditional statement
global conditional = cond("$type_variable" == "color", "red orange yellow green blue purple", "1 2 3 4 5 6")


* import stata data
use $input/stata_file.dta, clear 
* import csv data 
import delimited using $input/csv_file.csv, delimiter(",") varn(1) stringcols(2) clear
* import excel data
import excel using $input/excel_file.xlsx, clear


/*----- Calculation Tips -----*/
* generate a variable that is a substring of another variable (useful for string dates you need to convert to datetime)
	gen year = substr(date_string, 1, 4)
	gen month = substr(date_string, 5, 2)
	foreach var of varlist year month {
		destring `var', replace
	}

* replace characters in a string variable
	replace string_variable = subinstr(string_variable, "-", "", .)

* rename variables to the first row of a dataset 
	foreach var of varlist * {
		rename `var' `=`var'[1]'
	}

* generate aggregation by a grouping
	bys type date: egen sum_variable = total(variable)

* collapse data 
	collapse (sum) variable (firstnm) variable2 newname=variable3, by(date type)

* generate a numerical id variable from a string variable
	encode type, gen(type_code)

* specify the data as a panel
	xtset type_code date 

* if else statement example
	if "$restriction" == "no" {
		drop if keep == .
	}
	else if "$restriction" == "yes" {
		keep if keep == 1
	}	
	else {
		di "Incorrect Global RESTRICTION specification."
	}

* drop duplicates of all variables
	duplicates drop
* tag duplicates of all variables 
	duplicates tag type date, gen(duplicates)

* drop duplicates of a subset of variables 
	duplicates drop type date, force
* tag duplicates of a subset of variables 
	duplicates tag type date, gen(duplicates_type_date)

* merging data 
	merge 1:1 date id using $input/merging_data, keepusing(variable_merge1 variable_merge2) keep(1 2 3) nogen 

* use capture to perform a command when it applies (if it doesnt, you won't get an error, the line just won't be applied)
	cap drop if variable1 == "NULL" 	// so if variable1 is actually numeric, this will still continue onto the next portion of code without giving an error

* for loop examples:
	
	* loop over variable list
		foreach var of varlist * {
			replace `var' = `var' / 1e+3
		}
		foreach var of varlist *_pctile {
			egen mean_`var' = mean(`var')
		}

	* loop oveer strings/words
		foreach color in red orange yellow {
			gen `color'_amt = amount if color_var == "`color'"
		}

	* loop over numbers
		local counter = 0
		forval i=1/10 {
			local counter = `counter' + `i'
		}
		di `counter'

* select all variables that are not in a subset of the variables 
	ds name date, not 
	foreach var of varlist `r(varlist)' {
		display "`var'"
	}

* have a column for nicknames and one cell can have a list (amanda, mandy, etc.) - convert to one name per cell
		*https://sites.google.com/site/facebooknamelist/namelist
	import delimited using $input/data.csv, clear

	split variable1, parse(,) gen(other)
	gen id = _n 
	reshape long other, i(id) j(name)
	drop name variable1
	drop if mi(other)
	rename other variable_new
	order variable_new, before(id)

* have a variable (like a name formatted as last name, first name) and want to separate into different columns
gen seppos = strpos(name, ",")
gen last_name = substr(name, 1, seppos - 1)
gen first_name = substr(name, seppos+2, .)

/*----- Summarizing Data Tips -----*/

* general variable summary 
	sum variable1 
* detailed variable summary 
	sum variable1, de

* tabulate the values of a variable 
	tab variable1

* determine the unique values of a variable
	unique variable1

* generate a local variable with the values of a variable 
	levelsof variable1, local(levels_of_var1)






/*----- Plotting Figures -----*/

* label variables (these will be used in the legend if legend labels aren't specified)
label variable y_variable "Amount, Dollars"

* twoway line graph
	tw (line y_variable date if type == "Type1", color(navy)) ///
		(line y_variable date if type == "Type2", color(cranberry) axis(2)), xtitle("X Axis Title") ///
		ylabel(-2(1)2) xlabel(`=tm(2000m1)'(1)`=tm(2000m12)', format(%tmMonCCYY)) ///  can also format as %tmCY
		ytitle("Y Axis Title", margin(r+2)) ytitle("Y Axis Title", axis(2) margin(l-2)) ///
		title("Graph Title") graphregion(color(white)) ///
		legend(label(1 "Type1") label(2 "Type2") span symxsize(*0.5) rows(1) pos(1)) scale(1.2)
	graph export $output/graph_name.png, replace

* twoway Area Graph 
	tw (area y_variable date, lwidth(thin) color(navy) fint(60)) ///
		(area y_variable_next date, lwidth(thin) color(navy) fint(70)) if date<tm(2020m1), ///
		xtitle("X title") ylabel(0(2)10) ytitle("Y title", margin(r+2)) ///
		ylabel(0(20)100, axis(2)) ytitle("Y Title2", axis(2) margin(l-2)) ///
		tlabel(`=m(2000m1)'(2)`=m(2001m12)',format(%tmCY)) tmtick(`=m(2000m1)'(1)`=m(2001m12)') ///
		legend(symxsize(*.5)) scale(0.7) 
