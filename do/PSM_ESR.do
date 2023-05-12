
*** 
* Author: Wataru Kodama
* Purpose: GSR analysis
* Date: Nov 4, 2020
* Update: Dec 19, 2021
* Note: 
***

global root "stata"
global root2 ""
global xvars "rainfed upland ce_method gender education farm_exper hhsize organization ln_distance"


**************************************************************************
* Preparation 
**************************************************************************
use "https://raw.githubusercontent.com/Wataru-jp/GSR-Impact-assessment-in-Mozambique/main/data/GSR_QOpen.dta", clear
label var yield "Yield (kg/ha)"
label var tcost "Total cost (MZN/ha)"
label var area "Area (ha)"
label var seed_input "Seed input"
label var fert_input "Fertilizer input"
label var paid_labor_cost "Hired labor input"
label var irrigation "Irrigated lowland"
label var rainfed "Rainfed lowland"
label var upland "Upland"
label var ce_method "Transplanted rice"
label var soil1 "Clay soil"
label var soil2 "Loam soil"
label var soil3 "Other soil"
label var gender "Male household head (=1, if yes)"
label var education "Education (year)"
label var farm_exper "Farm experience (year)"
label var hhsize "Household size (\#)"
label var organization "Organization (=1, if yes)"
label var distance "Distance from seed source (minutes)"
label var distance2 "Distance from extension office (minutes)"

gen cost = (tcost/yield)*1000
label var cost "Cost efficiency (MZN/1000kg)"
gen ln_yield=log(yield)
label var ln_yield "Yield (Kg/ha), log"
gen ln_cost = log(tcost/yield+1)
label var ln_cost "Cost efficiency (MZN/kg), log"
gen ln_distance=log(distance+1)
label var ln_distance "Distance (seed source), log"
gen ln_distance2=log(distance2+1)
label var ln_distance2 "Distance (extension office), log"
gen ln_seed=log(seed_cost+1)
label var ln_seed "Seed, log"
gen ln_fert=log(fert_cost+1)
label var ln_fert "Fertilizer, log"
gen ln_pest=log(pest_cost+1)
label var ln_pest "Pesticide, log"
gen ln_labor=log(paid_labor_cost+1)
label var ln_labor "Hired labor, log"
gen ln_area=log(area)
label var ln_area "Area, log"
tempfile GSR
save `GSR', replace


**************************************************************************
*Propensity Score Matching
**************************************************************************
* Gaza province 
use `GSR', clear
keep if province=="Gaza"
psmatch2 GSR rainfed upland soil1 ce_method gender education farm_exper hhsize, n(1) logit noreplacement com
gen pair = _id if _treated==0
replace pair = _n1 if _treated==1
bysort pair: egen paircount = count(pair)
gen match=(paircount==2)
rename _pscore pscore1
rename match match1
keep hhid pscore match
tempfile temp1
save `temp1', replace

* Nampula province 
use `GSR', clear
keep if province=="Nampula"
psmatch2 GSR rainfed upland soil1 soil2 ce_method gender education farm_exper hhsize, n(1) logit noreplacement com
gen pair = _id if _treated==0
replace pair = _n1 if _treated==1
bysort pair: egen paircount = count(pair)
gen match=(paircount==2)
rename _pscore pscore2
rename match match2
keep hhid pscore match
tempfile temp2
save `temp2', replace

* Sofala province 
use `GSR', clear
keep if province=="Sofala"
psmatch2 GSR rainfed upland soil1 soil2 ce_method gender education farm_exper hhsize, n(1) logit noreplacement com
gen pair = _id if _treated==0
replace pair = _n1 if _treated==1
bysort pair: egen paircount = count(pair)
gen match=(paircount==2)
replace match=1 if GSR==0 & _pscore>0.5
rename _pscore pscore3
rename match match3
keep hhid pscore match
tempfile temp3
save `temp3', replace

* Combine samples
use `GSR', clear
merge 1:1 hhid using `temp1'
drop _merge
merge 1:1 hhid using `temp2'
drop _merge
merge 1:1 hhid using `temp3'
drop _merge
gen pscore = pscore1
replace pscore = pscore2 if pscore==.
replace pscore = pscore3 if pscore==.
gen match = match1
replace match = match2 if match==.
replace match = match3 if match==.
drop pscore1 pscore2 pscore3 match1 match2 match3
tempfile temp
save `temp', replace



**************************************************************************
*Endogenous switching regression: Yield
**************************************************************************
keep if match==1
*****************
* ALL
*****************
eststo ESR1: movestay ln_yield ln_seed ln_fert ln_labor ln_area rainfed upland ce_method soil1 soil2, ///
select(GSR = ln_distance ln_distance2) nolog iterate(100)
estadd local wald = round(e(chi2), .001)
estadd local indp = round(e(chi2_c), .001)
mspredict y1_T_TG, yc1_1
mspredict y1_UT_TG, yc2_1
mspredict y1_T_UTG, yc1_2
mspredict y1_UT_UTG, yc2_2
mspredict IMR1_1, mills1
mspredict IMR1_2, mills2

*****************
* Gaza
*****************
eststo ESR2: movestay ln_yield ln_seed ln_fert ln_labor ln_area rainfed ce_method soil1 ///
	if Gaza == 1, ///
select(GSR = ln_distance ln_distance2) nolog iterate(100)
estadd local wald = round(e(chi2), .001)
estadd local indp = round(e(chi2_c), .001)
mspredict y1_T_TG1, yc1_1
mspredict y1_UT_TG1, yc2_1
mspredict y1_T_UTG1, yc1_2
mspredict y1_UT_UTG1, yc2_2
mspredict IMR1_11, mills1
mspredict IMR1_21, mills2

*****************
* Nampula
*****************
eststo ESR3: movestay ln_yield ln_seed ln_labor ln_area rainfed upland ce_method soil1 soil2 ///
	if Nampula == 1, ///
select(GSR = ln_distance ln_distance2) nolog iterate(100)
estadd local wald = round(e(chi2), .001)
estadd local indp = round(e(chi2_c), .001)
mspredict y1_T_TG2, yc1_1
mspredict y1_UT_TG2, yc2_1
mspredict y1_T_UTG2, yc1_2
mspredict y1_UT_UTG2, yc2_2
mspredict IMR1_12, mills1
mspredict IMR1_22, mills2

*****************
* Sofala
*****************
eststo ESR4: movestay ln_yield ln_seed ln_labor ln_area rainfed upland ce_method soil1 soil2 ///
	if Sofala == 1, ///
select(GSR = ln_distance ln_distance2) nolog iterate(100)
estadd local wald = round(e(chi2), .001)
estadd local indp = round(e(chi2_c), .001)
mspredict y1_T_TG3, yc1_1
mspredict y1_UT_TG3, yc2_1
mspredict y1_T_UTG3, yc1_2
mspredict y1_UT_UTG3, yc2_2
mspredict IMR1_13, mills1
mspredict IMR1_23, mills2

**************************************************************************
*Endogenous switching regression: Cost efficiency
**************************************************************************
*****************
* ALL
*****************
eststo ESR5: movestay ln_cost ln_seed ln_fert ln_labor ln_area rainfed upland ce_method soil1 soil2, ///
select(GSR = ln_distance ln_distance2) nolog iterate(100)
estadd local wald = round(e(chi2), .001)
estadd local indp = round(e(chi2_c), .001)
mspredict y2_T_TG, yc1_1
mspredict y2_UT_TG, yc2_1
mspredict y2_T_UTG, yc1_2
mspredict y2_UT_UTG, yc2_2
mspredict IMR2_1, mills1
mspredict IMR2_2, mills2

*****************
* Gaza
*****************
eststo ESR6: movestay ln_cost ln_seed ln_fert ln_labor ln_area rainfed ce_method soil1 ///
	if Gaza == 1, ///
select(GSR = ln_distance ln_distance2) nolog iterate(100)
estadd local wald = round(e(chi2), .001)
estadd local indp = round(e(chi2_c), .001)
mspredict y2_T_TG1, yc1_1
mspredict y2_UT_TG1, yc2_1
mspredict y2_T_UTG1, yc1_2
mspredict y2_UT_UTG1, yc2_2
mspredict IMR2_11, mills1
mspredict IMR2_21, mills2

*****************
* Nampula
*****************
eststo ESR7: movestay ln_cost ln_seed ln_labor ln_area rainfed upland ce_method soil1 soil2 ///
	if Nampula == 1, ///
select(GSR = ln_distance ln_distance2) nolog iterate(100)
estadd local wald = round(e(chi2), .001)
estadd local indp = round(e(chi2_c), .001)
mspredict y2_T_TG2, yc1_1
mspredict y2_UT_TG2, yc2_1
mspredict y2_T_UTG2, yc1_2
mspredict y2_UT_UTG2, yc2_2
mspredict IMR2_12, mills1
mspredict IMR2_22, mills2

*****************
* Sofala
*****************
eststo ESR8: movestay ln_cost ln_seed ln_labor ln_area rainfed upland ce_method soil1 soil2 ///
	if Sofala == 1, ///
select(GSR = ln_distance ln_distance2) nolog iterate(100)
estadd local wald = round(e(chi2), .001)
estadd local indp = round(e(chi2_c), .001)
mspredict y2_T_TG3, yc1_1
mspredict y2_UT_TG3, yc2_1
mspredict y2_T_UTG3, yc1_2
mspredict y2_UT_UTG3, yc2_2
mspredict IMR2_13, mills1
mspredict IMR2_23, mills2

***************
* Results
***************
* combine provinces
forvalues n = 1/2 {
	foreach var in y`n'_T_TG y`n'_UT_TG y`n'_T_UTG y`n'_UT_UTG IMR`n'_1 IMR`n'_2 {
		replace `var'1 = . if Gaza == 0
		replace `var'2 = . if Nampula == 0
		replace `var'3 = . if Sofala == 0
	}
}

* ATT and ATU
forvalues n = 1/2 {
	gen ATT_y`n' = y`n'_T_TG - y`n'_UT_TG
	gen ATU_y`n' = y`n'_T_UTG - y`n'_UT_UTG
	gen pATT_y`n' = .
	gen pATU_y`n' = .
	qui su ATT_y`n'
	local att = r(mean) 
	qui su y`n'_UT_TG
	replace pATT_y`n' = `att'/r(mean)*100 if GSR == 1
	qui su ATU_y`n'
	local atu = r(mean) 
	qui su y`n'_UT_UTG
	replace pATU_y`n' = `atu'/r(mean)*100 if GSR == 0
	// province specific
	forvalues i = 1/3 {
		gen ATT_y`n'_`i' = y`n'_T_TG`i' - y`n'_UT_TG`i'
		gen ATU_y`n'_`i' = y`n'_T_UTG`i' - y`n'_UT_UTG`i'
		gen pATT_y`n'_`i' = .
		gen pATU_y`n'_`i' = .
		// % change
		qui su ATT_y`n'_`i'
		local att = r(mean) 
		qui su y`n'_UT_TG`i'
		replace pATT_y`n'_`i' = `att'/r(mean)*100 if GSR == 1 & province2 == `i'
		qui su ATU_y`n'_`i'
		local atu = r(mean) 
		qui su y`n'_UT_UTG`i'
		replace pATU_y`n'_`i' = `atu'/r(mean)*100 if GSR == 1 & province2 == `i'
	}	
}
 
su ATT_y1 ATU_y1 ATT_y1_1 ATU_y1_1 ATT_y1_2 ATU_y1_2 ATT_y1_3 ATU_y1_3
su pATT_y1 pATU_y1 pATT_y1_1 pATU_y1_1 pATT_y1_2 pATU_y1_2 pATT_y1_3 pATU_y1_3
su ATT_y2 ATU_y2 ATT_y2_1 ATU_y2_1 ATT_y2_2 ATU_y2_2 ATT_y2_3 ATU_y2_3
su pATT_y2 pATU_y2 pATT_y2_1 pATU_y2_1 pATT_y2_2 pATU_y2_2 pATT_y2_3 pATU_y2_3 

forvalues n = 1/2 {
	su y`n'_T_TG y`n'_UT_TG y`n'_T_UTG y`n'_UT_UTG 
	su y`n'_T_TG1 y`n'_UT_TG1 y`n'_T_UTG1 y`n'_UT_UTG1
	su y`n'_T_TG2 y`n'_UT_TG2 y`n'_T_UTG2 y`n'_UT_UTG2
	su y`n'_T_TG3 y`n'_UT_TG3 y`n'_T_UTG3 y`n'_UT_UTG3
}
forvalues n = 1/2 {
	su IMR`n'_11 IMR`n'_12 IMR`n'_13 if GSR ==1
	su IMR`n'_21 IMR`n'_22 IMR`n'_23 if GSR ==0
}
forvalues n = 1/2 {
	ttest IMR`n'_11 == 0
	ttest IMR`n'_12 == 0
	ttest IMR`n'_13 == 0
}
forvalues n = 1/2 {
	ttest IMR`n'_21 == 0
	ttest IMR`n'_22 == 0
	ttest IMR`n'_23 == 0
}
stop 

**************************************************************************
* List of tables
**************************************************************************
use "https://raw.githubusercontent.com/Wataru-jp/GSR-Impact-assessment-in-Mozambique/main/data/GSR_QOpen.dta", clear
// create variables 
gen cost = (tcost/yield)*1000
gen ln_yield=log(yield)
gen ln_cost = log(tcost/yield + 1)
gen ln_distance=log(distance + 1)
gen ln_distance2=log(distance2 + 1)
gen improved = (type==3)
gen traditional = (type==2)
gen GSR2 = 1 - GSR
// labels 
label var ln_yield "Yield, log"
label var ln_cost "Cost efficiency, log"
label var yield "Yield"
label var cost "Cost efficiency"
label var area "Area"
label var seed_input "Seed input"
label var fert_input "Fertilizer input"
label var paid_labor_cost "Hired labor input"
label var irrigation "Irrigated lowland"
label var rainfed "Rainfed lowland"
label var upland "Upland"
label var ce_method "Transplanted rice"
label var soil1 "Clay soil"
label var soil2 "Loam soil"
label var soil3 "Other soil"
label var gender "Male household head"
label var education "Education"
label var farm_exper "Farm experience"
label var hhsize "Household size"
label var gender "Male household head"
label var education "Education"
label var farm_exper "Farm experience"
label var hhsize "Household size"
label var distance "Distance from seed source"
label var distance2 "Distance from extension office"
label var ln_distance "Distance: seed source, log"
label var ln_distance2 "Distance: extension office, log"
label var GSR "GSR variety"
label var improved "Improved non-GSR variety"
label var traditional "Traditional variety"
tempfile GSR
save `GSR', replace

// PSM 
* Gaza province 
use `GSR', clear
keep if province=="Gaza"
psmatch2 GSR rainfed upland soil1 ce_method gender education farm_exper hhsize, n(1) noreplacement com
gen pair = _id if _treated==0
replace pair = _n1 if _treated==1
bysort pair: egen paircount = count(pair)
gen match=(paircount==2)
rename _pscore pscore1
rename match match1
keep hhid pscore match
tempfile temp1
save `temp1', replace
* Nampula province 
use `GSR', clear
keep if province=="Nampula"
psmatch2 GSR rainfed upland soil1 soil2 ce_method gender education farm_exper hhsize, n(1) logit noreplacement com
gen pair = _id if _treated==0
replace pair = _n1 if _treated==1
bysort pair: egen paircount = count(pair)
gen match=(paircount==2)
rename _pscore pscore2
rename match match2
keep hhid pscore match
tempfile temp2
save `temp2', replace
* Sofala province 
use `GSR', clear
keep if province=="Sofala"
psmatch2 GSR rainfed upland soil1 soil2 ce_method gender education farm_exper hhsize, n(1) logit noreplacement com
gen pair = _id if _treated==0
replace pair = _n1 if _treated==1
bysort pair: egen paircount = count(pair)
gen match=(paircount==2)
replace match=1 if GSR==0 & _pscore>0.5
rename _pscore pscore3
rename match match3
keep hhid pscore match
tempfile temp3
save `temp3', replace
* Combine samples
use `GSR', clear
merge 1:1 hhid using `temp1'
drop _merge
merge 1:1 hhid using `temp2'
drop _merge
merge 1:1 hhid using `temp3'
drop _merge
gen pscore = pscore1
replace pscore = pscore2 if pscore==.
replace pscore = pscore3 if pscore==.
gen match = match1
replace match = match2 if match==.
replace match = match3 if match==.
drop pscore1 pscore2 pscore3 match1 match2 match3
stop

***************
* Table 3
***************
eststo probit1: ///
	probit GSR rainfed upland soil1 ce_method gender education farm_exper hhsize ///
	if province=="Gaza"
eststo probit2: ///
	probit GSR rainfed upland soil1 soil2 ce_method gender education farm_exper hhsize ///
	if province=="Nampula"
eststo probit3: ///
	probit GSR rainfed upland soil1 soil2 ce_method gender education farm_exper hhsize ///
	if province=="Sofala"
** save 
label var GSR "GSR (Simão)"
esttab probit1 probit2 probit3 using $root/tex/table3_revised.tex, ///
	se label nogap b(%4.3f) nonotes nonumber mtitle("Gaza" "Nampula" "Sofala") ///
	order(rainfed upland soil1 soil2) ///
	s(chi2 r2_p N, fmt(%9.3f %9.3f %9.0g) ///
	labels("Log-likelihood $\chi^2$" "Pseudo $\R^2$" "Observations")) ///
	star(* 0.10 ** 0.05 *** 0.01) replace 
esttab probit1 probit2 probit3 using $root2/table3.csv, ///
	se label nogap b(%4.3f) nonotes nonumber mtitle("Gaza" "Nampula" "Sofala") ///
	order(rainfed upland soil1 soil2) ///
	s(chi2 r2_p N, fmt(%9.3f %9.3f %9.0g) ///
	labels("Log-likelihood chi^2" "Pseudo R^2" "Observations")) ///
	star(* 0.10 ** 0.05 *** 0.01) replace 
estimates drop probit1 probit2 probit3

***************
* Table 5A
***************
eststo col1: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Gaza == 1 & GSR==1
eststo col2: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 improved traditional if Gaza == 1 & GSR==0
eststo col3: qui estpost ttest yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Gaza == 1, by(GSR2) // Diff.
eststo col4: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance if Nampula == 1 & GSR==1
eststo col5: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 improved traditional if Nampula == 1 & GSR==0
eststo col6: qui estpost ttest yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Nampula == 1, by(GSR2) // Diff.
eststo col7: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Sofala == 1 & GSR==1
eststo col8: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 improved traditional if Sofala == 1 & GSR==0
eststo col9: qui estpost ttest yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Sofala == 1, by(GSR2) // Diff.
esttab col1 col2 col3 col4 col5 col6 col7 col8 col9 ///
using $root/tex/table5a_revised.tex, ///
label nogap nonotes cells("mean(pattern(1 1 0 1 1 0 1 1 0) fmt(3)) b(star pattern(0 0 1 0 0 1 0 0 1) fmt(3))" ) nomtitle nonumber ///
mgroups("Gaza province" "Nampula province" "Sofala province", ///
		pattern(1 0 0 1 0 0 1 0 0) ///
		span prefix(\multicolumn{@span}{c}{) suffix(}) ///
		erepeat(\cmidrule(lr){@span})) replace
esttab col1 col2 col3 col4 col5 col6 col7 col8 col9 ///
using $root2/table5A.csv, ///
label nogap nonotes cells("mean(pattern(1 1 0 1 1 0 1 1 0) fmt(3)) b(star pattern(0 0 1 0 0 1 0 0 1) fmt(3))" ) nomtitle nonumber ///
mgroups("Gaza province" "Nampula province" "Sofala province", ///
		pattern(1 0 0 1 0 0 1 0 0) span ) replace
for num 1/9: estimates drop colX

***************
* Table 5B
***************
eststo col1: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Gaza == 1 & GSR==1 & match==1
eststo col2: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 improved traditional if Gaza == 1 & GSR==0 & match==1
eststo col3: qui estpost ttest yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Gaza == 1 & match==1, by(GSR2) // Diff.
eststo col4: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Nampula == 1 & GSR==1 & match==1
eststo col5: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 improved traditional if Nampula == 1 & GSR==0 & match==1
eststo col6: qui estpost ttest yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Nampula == 1 & match==1, by(GSR2) // Diff.
eststo col7: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Sofala == 1 & GSR==1 & match==1
eststo col8: qui estpost su yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 improved traditional if Sofala == 1 & GSR==0 & match==1
eststo col9: qui estpost ttest yield cost seed_input fert_input paid_labor_cost irrigation rainfed upland ce_method soil1 soil2 soil3 gender education farm_exper hhsize distance distance2 if Sofala == 1 & match==1, by(GSR2) // Diff.
esttab col1 col2 col3 col4 col5 col6 col7 col8 col9 ///
	using $root/tex/table5b_revised.tex, ///
	label nogap nonotes nomtitle nonumber ///
	cells("mean(pattern(1 1 0 1 1 0 1 1 0) fmt(3)) b(star pattern(0 0 1 0 0 1 0 0 1) fmt(3))" ) ///
	mgroups("Gaza province" "Nampula province" "Sofala province", ///
		pattern(1 0 0 1 0 0 1 0 0) ///
		span prefix(\multicolumn{@span}{c}{) suffix(}) ///
		erepeat(\cmidrule(lr){@span})) replace
esttab col1 col2 col3 col4 col5 col6 col7 col8 col9 ///
	using $root2/table5B.csv, ///
	label nogap nonotes nomtitle nonumber ///
	cells("mean(pattern(1 1 0 1 1 0 1 1 0) fmt(3)) b(star pattern(0 0 1 0 0 1 0 0 1) fmt(3))" ) ///
	mgroups("Gaza province" "Nampula province" "Sofala province", ///
		pattern(1 0 0 1 0 0 1 0 0) span) replace
for num 1/9: estimates drop colX

********************
*Table 6A and 7A
********************
esttab ESR1 using $root/tex/table6a_revised2.tex, ///
	prehead("{\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi} \begin{tabular}{l*{3}{c}} \hline") ///
	posthead("Yield, log &	\multicolumn{2}{c}{Regime equation} & \multicolumn{1}{c}{Selection}	 \\ \cline{2-3} &\multicolumn{1}{c}{GSR(Simão)} & \multicolumn{1}{c}{Non-GSR} & \multicolumn{1}{c}{equation}\\ \hline") ///
	keep(GSR: ln_yield_1: ln_yield_0:) ///
	s(wald indp N, fmt(%9.3f %9.3f %9.0g) ///
	labels("Wald test" "LR test of indep. eqns." "Observations")) ///
	unstack se label nogap nonotes nomtitles nonumber b(%4.3f) replace
esttab ESR1 using $root2/table6A.csv, ///
	prehead("Regime equation") ///
	posthead("Regime equation") ///
	keep(GSR: ln_yield_1: ln_yield_0:) ///
	s(wald indp N, fmt(%9.3f %9.3f %9.0g) ///
	labels("Wald test" "LR test of indep. eqns." "Observations")) ///
	unstack se label nogap nonotes nomtitles nonumber b(%4.3f) replace


esttab ESR5 using $root/tex/table7a_revised2.tex, ///
	prehead("{\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi} \begin{tabular}{l*{3}{c}} \hline") ///
	posthead("Cost efficiency, log &	\multicolumn{2}{c}{Regime equation} & \multicolumn{1}{c}{Selection}	 \\ \cline{2-3} &\multicolumn{1}{c}{GSR(Simão)} & \multicolumn{1}{c}{Non-GSR} & \multicolumn{1}{c}{equation}\\ \hline") ///
	keep(GSR: ln_cost_1: ln_cost_0:) ///
	s(wald indp N, fmt(%9.3f %9.3f %9.0g) ///
	labels("Wald test" "LR test of indep. eqns." "Observations")) ///
	unstack se label nogap nonotes b(%4.3f) replace
esttab ESR5 using $root2/table7A.csv, ///
	prehead("Regime equation") ///
	posthead("Regime equation") ///
	keep(GSR: ln_cost_1: ln_cost_0:) ///
	s(wald indp N, fmt(%9.3f %9.3f %9.0g) ///
	labels("Wald test" "LR test of indep. eqns." "Observations")) ///
	unstack se label nogap nonotes nomtitles nonumber b(%4.3f) replace

********************
*Table A4
********************
eststo fal1: probit GSR ln_distance ln_distance2
estadd local wald = round(e(chi2), .001)
eststo fal2: reg ln_yield ln_distance ln_distance2 if GSR == 0 
estadd local wald = round(e(F), .001)
eststo fal3: reg ln_cost ln_distance ln_distance2 if GSR == 0
estadd local wald = round(e(F), .001)
eststo fal4: probit GSR ln_distance ln_distance2 if match==1
estadd local wald = round(e(chi2), .001)
eststo fal5: reg ln_yield ln_distance ln_distance2 if GSR == 0 & match == 1
estadd local wald = round(e(F), .001)
eststo fal6: reg ln_cost ln_distance ln_distance2 if GSR == 0 & match == 1
estadd local wald = round(e(F), .001)
** save 
esttab fal1 fal2 fal3 using $root/tex/table_A4_revised.tex, ///
	prehead("{\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi} \begin{tabular}{l*{3}{cc}} \hline\hline") ///
	posthead("\hline \multicolumn{7}{c}{\bf{Before matching}}\\") ///	
	fragment ///
	mgroups("GSR adoption" "Yield, log" "Cost efficiency, log", ///
	pattern(1 1 1) span prefix(\multicolumn{@span}{c}{) suffix(})) ///
	se label nogap b(%4.3f) nomtitles nonumbers nolines wide ///
	s(wald N, fmt(%9.3f %9.0g) labels("Wald test" "Observations")) ///
	star(* 0.10 ** 0.05 *** 0.01) prefoot("\hline") replace 
esttab fal4 fal5 fal6 using $root/tex/table_A4_revised.tex, ///
	posthead("\hline \multicolumn{7}{c}{\bf{After matching}}\\") ///
	fragment append ///
	se label nogap b(%4.3f) nomtitles nonumbers nolines wide ///
	s(wald N, fmt(%9.3f %9.0g) labels("Wald test" "Observations")) ///
	star(* 0.10 ** 0.05 *** 0.01) prefoot("\hline") ///
	postfoot("\hline\hline \end{tabular}}") 
for num 1/6: estimates drop falX

********************
*Table A5 and A6
********************
esttab ESR2 ESR3 ESR4 using $root/tex/table_A5_revised2.tex, ///
	prehead("{\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi} \begin{tabular}{l*{9}{c}}\hline\hline") ///
	posthead("&\multicolumn{3}{c}{\bf{Gaza Province}} & \multicolumn{3}{c}{\bf{Nampula Province}} & \multicolumn{3}{c}{\bf{Sofala Province}} \\\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10} Yield, log & \multicolumn{2}{c}{Regime equation} &	\multicolumn{1}{c}{Selection} & \multicolumn{2}{c}{Regime equation} & \multicolumn{1}{c}{Selection} & \multicolumn{2}{c}{Regime equation} & \multicolumn{1}{c}{Selection} \\\cline{2-3}\cline{5-6}\cline{8-9} & \multicolumn{1}{c}{GSR (Simão)} & \multicolumn{1}{c}{Non-GSR} & \multicolumn{1}{c}{equation} & \multicolumn{1}{c}{GSR (Simão)} & \multicolumn{1}{c}{Non-GSR} & \multicolumn{1}{c}{equation} & \multicolumn{1}{c}{GSR (Simão)} & \multicolumn{1}{c}{Non-GSR} & \multicolumn{1}{c}{equation} \\ \hline") ///
	keep(GSR: ln_yield_1: ln_yield_0:) ///
	order(ln_seed ln_fert ln_labor ln_area rainfed upland ce_method soil1 soil2) ///
	s(wald indp N, fmt(%9.3f %9.3f %9.0g) ///
	labels("Wald test" "LR test of indep. eqns." "Observations")) ///
	unstack se label nogap nonotes nomtitles nonumber b(%4.3f) replace
esttab ESR6 ESR7 ESR8 using $root/tex/table_A6_revised2.tex, ///
	prehead("{\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi} \begin{tabular}{l*{9}{c}}\hline\hline") ///
	posthead("&\multicolumn{3}{c}{\bf{Gaza Province}} & \multicolumn{3}{c}{\bf{Nampula Province}} & \multicolumn{3}{c}{\bf{Sofala Province}} \\\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10} Cost efficiency, log & \multicolumn{2}{c}{Regime equation} &	\multicolumn{1}{c}{Selection} & \multicolumn{2}{c}{Regime equation} & \multicolumn{1}{c}{Selection} & \multicolumn{2}{c}{Regime equation} & \multicolumn{1}{c}{Selection} \\\cline{2-3}\cline{5-6}\cline{8-9} & \multicolumn{1}{c}{GSR (Simão)} & \multicolumn{1}{c}{Non-GSR} & \multicolumn{1}{c}{equation} & \multicolumn{1}{c}{GSR (Simão)} & \multicolumn{1}{c}{Non-GSR} & \multicolumn{1}{c}{equation} & \multicolumn{1}{c}{GSR (Simão)} & \multicolumn{1}{c}{Non-GSR} & \multicolumn{1}{c}{equation} \\ \hline") ///
	keep(GSR: ln_cost_1: ln_cost_0:) ///
	order(ln_seed ln_fert ln_labor ln_area rainfed upland ce_method soil1 soil2) ///
	s(wald indp N, fmt(%9.3f %9.3f %9.0g) ///
	labels("Wald test" "LR test of indep. eqns." "Observations")) ///
	unstack se label nogap nonotes b(%4.3f) replace
for num 1/6: estimates drop ESRX

**************************************************************************
* List of figures 
**************************************************************************
***************
* Figure 3
***************
// histogram 
set scheme s1mono
twoway (hist pscore if GSR==1, width(0.025) color(orange%40)) ///        
       (hist pscore if GSR==0, width(0.025) color(midblue%50)) /// 
	   (kdensity pscore if GSR==1, color(dkorange)) ///
	   (kdensity pscore if GSR==0, color(midblue)), ///
       graphregion(color(white)) legend(order(1 "Adoptoer" 2 "Non-adoptoer")) ytitle(Frequency) xtitle("")
graph export $root/tex/fig3a.jpg, as(jpg) name("Graph") quality(100) replace
twoway (hist pscore if GSR==1 & match==1, width(0.025) color(orange%40)) ///        
       (hist pscore if GSR==0 & match==1, width(0.025) color(midblue%50)) ///   
	   (kdensity pscore if GSR==1 & match==1, color(dkorange)) ///
	   (kdensity pscore if GSR==0 & match==1, color(midblue)), ///
       graphregion(color(white)) legend(order(1 "Adoptoer" 2 "Non-adoptoer")) ytitle(Frequency) xtitle("")
graph export $root/tex/fig3b.jpg, as(jpg) name("Graph") quality(100) replace

//EPS format
twoway (hist pscore if GSR==0, width(0.025) fcolor(gs14) lcolor(blue)) /// 
	   (hist pscore if GSR==1, width(0.025) fcolor(none) lcolor(dkorange)) ///        
	   (kdensity pscore if GSR==1, color(dkorange)) ///
	   (kdensity pscore if GSR==0, color(blue)), ///
       graphregion(color(white)) legend(order(1 "Non-adopter" 2 "Adopter")) ytitle(Frequency) xtitle("")
graph export $root2/eps-figures/fig3a.eps, as(eps) name("Graph") replace
twoway (hist pscore if GSR==0 & match==1, width(0.025) color(gs14) lcolor(blue)) /// 
	   (hist pscore if GSR==1 & match==1, width(0.025) color(none) lcolor(dkorange)) ///        
	   (kdensity pscore if GSR==1 & match==1, color(dkorange)) ///
	   (kdensity pscore if GSR==0 & match==1, color(blue)), ///
       graphregion(color(white)) legend(order(1 "Non-adopter" 2 "Adopter")) ytitle(Frequency) xtitle("")
graph export $root2/eps-figures/fig3b.eps, as(eps) name("Graph") replace


***************
* Figure A1 
***************
gen cost = (tcost/yield)*1000
gen ln_yield = log(yield)
gen ln_cost = log(tcost/yield+1)
// histogram 
set scheme s1mono
twoway (hist ln_yield if GSR == 1, width(0.2) color(gray%50)) ///        
	   (kdensity ln_yield if GSR == 1, color(black%50)), ///
       graphregion(color(white)) legend(off) ytitle(Frequency) xtitle("")
graph export $root/tex/figA1a.jpg, as(jpg) name("Graph") quality(100) replace
twoway (hist ln_yield if GSR == 0, width(0.2) color(gray%50)) /// 
	   (kdensity ln_yield if GSR == 0, color(black%50)), ///
       graphregion(color(white)) legend(off) ytitle(Frequency) xtitle("")
graph export $root/tex/figA1b.jpg, as(jpg) name("Graph") quality(100) replace
twoway (hist ln_cost if GSR == 1, width(0.2) color(gray%50)) ///        
	   (kdensity ln_cost if GSR == 1, color(black%50)), ///
       graphregion(color(white)) legend(off) ytitle(Frequency) xtitle("")
graph export $root/tex/figA1c.jpg, as(jpg) name("Graph") quality(100) replace
twoway (hist ln_cost if GSR == 0, width(0.2) color(gray%50)) /// 
	   (kdensity ln_cost if GSR == 0, color(black%50)), ///
       graphregion(color(white)) legend(off) ytitle(Frequency) xtitle("")
graph export $root/tex/figA1d.jpg, as(jpg) name("Graph") quality(100) replace

//EPS format
set scheme s1mono
twoway (hist ln_yield if GSR == 1, width(0.2) color(gs14) lcolor(gs7)) ///        
	   (kdensity ln_yield if GSR == 1, color(gs1)), ///
       graphregion(color(white)) legend(off) ytitle(Frequency) xtitle("")
graph export $root2/eps-figures/figA1a.eps, as(eps) name("Graph") replace
twoway (hist ln_yield if GSR == 0, width(0.2) color(gs14) lcolor(gs7)) /// 
	   (kdensity ln_yield if GSR == 0, color(gs1)), ///
       graphregion(color(white)) legend(off) ytitle(Frequency) xtitle("")
graph export $root2/eps-figures/figA1b.eps, as(eps) name("Graph") replace
twoway (hist ln_cost if GSR == 1, width(0.2) color(gs14) lcolor(gs7)) ///        
	   (kdensity ln_cost if GSR == 1, color(gs1)), ///
       graphregion(color(white)) legend(off) ytitle(Frequency) xtitle("")
graph export $root2/eps-figures/figA1c.eps, as(eps) name("Graph") replace
twoway (hist ln_cost if GSR == 0, width(0.2) color(gs14) lcolor(gs7)) /// 
	   (kdensity ln_cost if GSR == 0, color(gs1)), ///
       graphregion(color(white)) legend(off) ytitle(Frequency) xtitle("")
graph export $root2/eps-figures/figA1d.eps, as(eps) name("Graph") replace

***************
* Figure 1
***************
use "https://raw.githubusercontent.com/Wataru-jp/GSR-Impact-assessment-in-Mozambique/main/data/MZ_Rice.dta", clear
set scheme s1mono
gen self_eff = production/consumption*100
graph twoway ///
	bar self_eff time, yaxis(2) ytitle("Kilograms/ Percents",axis(2)) bcolor(gs14) ///
	|| line production consumption Improt time, ///
		yaxis(1) ytitle("Thousand tons",axis(1)) lcolor(green cranberry ebblue) ///
	|| line consumption_capita time, yaxis(2) graphregion(color(gs16)) xlabel(1990(5)2021) ///
		xtitle("") lcolor(black) lpattern(dash) ///
	legend(order(4 "Consumption" 2 "Consumption per capita (kg)" 3 "Production" 1 "Self-sufficiency ratio (%)" 5 "Imports" ) size(medsmall))
graph export $root2/Kodama-et-al_eps-figure/fig1.eps, as(eps) name("Graph") replace
*graph export $root/tex/fig1.jpg, as(jpg) name("Graph") quality(100) replace

