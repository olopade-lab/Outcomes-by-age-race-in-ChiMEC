use "C:\Users\eterman\OneDrive - The University of Chicago\CRDW_Registry_IndexDx_2021Mar18_coded_newJuly8.dta", clear
*group based on ages: older >= 55 and younger <= 40
egen agegrp = cut(Age), at(19 41 55 100) icodes 


tab agegrp raceethnic
drop if sex =="1 Male"
*create 4 target groups: Young White, Young Black, Older White, Older Black
*variable ageethnic = 1, 2, 3, 4 respectively^

gen ageethnic = 1 if agegrp==0 & raceethnic == 1 & sex == "2 Female"


 replace ageethnic = 2 if agegrp == 0 & raceethnic == 2 & sex == "2 Female"

 replace ageethnic = 3 if agegrp == 2 & raceethnic == 1 & sex == "2 Female"

 replace ageethnic = 4 if agegrp == 2 & raceethnic == 2 & sex == "2 Female"

*recode so ageethnic displays descriptors:
label define ageethnic_label 1 "young White" 2 "young Black" 3 "older White" 4 "older Black"
label values ageethnic ageethnic_label
tab ageethnic

** select only invasive early breast cancer: we want to use patholgical stage for patients with adjuvant treatment 
keep if inlist(stage1, 1, 2, 3)

egen CCI3 = cut(CCI), at(0 1 2 10)
egen grade3 = cut(grade), at(1 3 5) icodes 

tab ageethnic

***CHECK PCR VARIABLE
tab pCR

tab ageethnic pCR, row 
table ageethnic pCR subtype 
table ageethnic subtype , c(mean pCR) format(%5.3f)


*LOG REGRESSION WITH PCR

*unadjusted:
logistic pCR ib3.ageethnic

 * young black vs young white 
 lincom 2.ageethnic - 1.ageethnic , or
 *young black vs old black 
 lincom 2.ageethnic - 4.ageethnic , or 

 *adjusted for subtype, stage:
 
logistic pCR ib3.ageethnic ib1.subtype i.stage1 
* young black vs young white 
  lincom 2.ageethnic - 1.ageethnic , or 
* young black vs old black 
  lincom 2.ageethnic - 4.ageethnic , or
  
 *adjusted for subtype, stage, grade:
logistic pCR ib3.ageethnic ib1.subtype i.stage1 i.grade3

* young black vs young white 
  lincom 2.ageethnic - 1.ageethnic , or
* young black vs old black 
  lincom 2.ageethnic - 4.ageethnic , or 

******************
*SET UP FOR SURVIVAL FUNCTIONS/KM GRAPHS 
*overall survival, event is death
*need to set up for a piecewise Cox model b/c proportional hazard assumption doesn't hold otherwise 
gen ID= _n
stset year_FU, failure (dead) id(ID)

** show # of patients included and # of events 
tab ageethnic dead
sts test (ageethnic)
sts list, at(0 5 10) by(ageethnic)

stsplit fu2, at(5 100)   
  stset
  gen exptime = _t - _t0
  tab fu2 _d
  stsum , by(fu2)

 stcox ib3.ageethnic##fu2 

 *to get older black vs. older white for >5 years
  lincom 4.ageethnic + 4.ageethnic#5.fu2, hr 
  
 *need to get young black vs. young white for >5 years
 lincom (2.ageethnic - 1.ageethnic) + (2.ageethnic - 1.ageethnic)#5.fu2 , hr 
 
 *young black vs. old black for >5 years
 
 
 * young black vs young white 
        lincom 2.ageethnic - 1.ageethnic , hr 
* young black vs old black 
        lincom 2.ageethnic - 4.ageethnic , hr
    estat phtest, detail
*test to see that cox proportional hazard model fits

*get median survival time 
stsum 

sts test (ageethnic)
stcox ib(frequent).ageethnic 
testparm i.ageethnic 
* young black vs young white 
        lincom 2.ageethnic - 1.ageethnic , hr 
* young black vs old black 
        lincom 2.ageethnic - 4.ageethnic , hr 

* multivariable analysis: different way of controlling molecular markers 
stcox ib3.ageethnic i.ER i.PR i.her2
testparm i.ageethnic

stcox ib3.ageethnic ib1.subtype
testparm i.ageethnic

stcox ib3.ageethnic ib1.subtype i.PR
testparm i.ageethnic

* multivariable analysis: final model I    
stcox ib3.ageethnic ib1.subtype i.grade i.stage1 i.CCI3
testparm i.ageethnic
* young black vs young white 
     lincom 2.ageethnic - 1.ageethnic , hr 
* young black vs old black 
    lincom 2.ageethnic - 4.ageethnic , hr 

	stcox c.Age##ib1.raceethnic if inlist(raceethnic, 1, 2)  
 
  lincom 2.raceethnic#c.Age + c.Age, hr
 
 generate OldYoung = 0 if agegrp == 0
replace OldYoung = 1 if agegrp == 2
  stcox i.OldYoung##ib1.raceethnic if inlist(raceethnic, 1, 2)  
 
  lincom 1.OldYoung#2.raceethnic + 1.OldYoung, hr
  
estat phtest, detail


*make KM graphs:

*make KM graphs: use tmax option to cut time display (but all data are used), do not use "if condition" as the later removse patients with short follow-up. 
.  sts graph, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Overall Survival) xtitle (Time from diagnosis (years)) ytitle (Probability of survival) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 

*by ageethnic categories:
 sts graph if ageethnic ==1, by(subtype) tmax(15) risktable(, order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC")) title(Overall survival in young white women) xtitle (Time from diagnosis (years)) ytitle (Probability of survival) legend(order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
 
 sts graph if ageethnic ==2, by(subtype) tmax(15) risktable(, order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC")) title(Overall survival in young black women) xtitle (Time from diagnosis (years)) ytitle (Probability of survival) legend(order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
 
 sts graph if ageethnic ==3, by(subtype) tmax(15) risktable(, order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC")) title(Overall survival in older white women) xtitle (Time from diagnosis (years)) ytitle (Probability of survival) legend(order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
 
  sts graph if ageethnic ==4, by(subtype) tmax(15) risktable(, order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC")) title(Overall survival in older black women) xtitle (Time from diagnosis (years)) ytitle (Probability of survival) legend(order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
  
 *by subtype 
 sts graph if subtype ==1, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Overall Survival HR+) xtitle (Time from diagnosis (years)) ytitle (Probability of survival) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
  
  sts graph if subtype ==2 , by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Overall Survival HER2+) xtitle (Time from diagnosis (years)) ytitle (Probability of survival) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
  
  sts graph if subtype ==3, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Overall Survival HR-/HER2+) xtitle (Time from diagnosis (years)) ytitle (Probability of survival) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
  
  sts graph if subtype ==4, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Overall Survival TNBC) xtitle (Time from diagnosis (years)) ytitle (Probability of survival) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
  
*********RECURRENCE FREE SURVIVAL event is recurrence
**Recurrence-free survival (disease-free survival)
*Recurrence-free survival, event includes both death and recurrence 

stset year_RFS, failure (RFS) id(ID)

** show # of patients included and # of events 
tab ageethnic RFS
sts test (ageethnic)
sts list, at(0 5 10) by(ageethnic)

stsplit fu2, at(5 100)   
  stset
  gen exptime = _t - _t0
  
  tab fu2 _d
  stsum , by(fu2)

 stcox ib3.ageethnic##fu2 

** show # of patients included and # of events 
  tab ageethnic RFS
  sts test (ageethnic)
  
*categorical age:


  stcox i.OldYoung##ib1.raceethnic if inlist(raceethnic, 1, 2)  
 
  lincom 1.OldYoung#2.raceethnic + 1.OldYoung, hr
  
Here, you can code OldYoung as 0 for young and 1 for old.

*continuous age:
 
stcox c.Age##ib1.raceethnic if inlist(raceethnic, 1, 2)  
 
  lincom 2.raceethnic#c.Age + c.Age*10, hr

stcox ib(frequent).ageethnic 
testparm i.ageethnic 
* young black vs young white 
       lincom 2.ageethnic - 1.ageethnic , hr
* young black vs old black 
       lincom 2.ageethnic - 4.ageethnic , hr 
estat phtest, detail

* multivariable analysis: different way of controlling molecular markers 
   stcox ib3.ageethnic i.ER i.PR i.her2  
testparm i.ageethnic 
stcox ib3.ageethnic ib1.subtype
testparm i.ageethnic 
stcox ib3.ageethnic ib1.subtype i.PR 
testparm i.ageethnic

* multivariable analysis: final model II
  stcox ib3.ageethnic ib1.subtype i.grade i.stage1 i.CCI3
testparm i.ageethnic  

 * young black vs young white 
     lincom 2.ageethnic - 1.ageethnic , hr
* young black vs old black 
       lincom 2.ageethnic - 4.ageethnic , hr 
estat phtest, detail 


 *make KM graphs: use tmax option to cut time display (but all data are used), do not use "if condition" as the later remove patients with short follow-up. 
.  sts graph, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Disease-free Survival) xtitle (Time from diagnosis (years)) ytitle (Disease-free survival prob.) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal)  ) xsize(7) ysize(6) 
  
*by ageethnic categories:
 sts graph if ageethnic ==1, by(subtype) tmax(15) risktable(, order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC")) title(Disease-free survival in young white women) xtitle (Time from diagnosis (years)) ytitle (Disease-free survival prob.) legend(order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
 
 sts graph if ageethnic ==2, by(subtype) tmax(15) risktable(, order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC")) title(Disease-free survival in young black women) xtitle (Time from diagnosis (years)) ytitle (Disease-free survival prob.) legend(order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
 
 sts graph if ageethnic ==3, by(subtype) tmax(15) risktable(, order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC")) title(Disease-free survival in older white women) xtitle (Time from diagnosis (years)) ytitle (Disease-free survival prob.) legend(order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
 
  sts graph if ageethnic ==4, by(subtype) tmax(15) risktable(, order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC")) title(Disease-free survival in older black women) xtitle (Time from diagnosis (years)) ytitle (Disease-free survival prob.) legend(order(1 "HR+/HER2-" 2 "HR+/HER2+" 3 "HR-/HER2+" 4 "TNBC") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
  
*by subtype 
 sts graph if subtype ==1, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Disease-free survival HR+) xtitle (Time from diagnosis (years)) ytitle (Disease-free survival prob.) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
  
  sts graph if subtype ==2 | subtype ==3, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Disease-free survival HER2+) xtitle (Time from diagnosis (years)) ytitle (Disease-free survival prob.) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
  
sts graph if subtype ==3, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Disease-free survival HR-/HER2+) xtitle (Time from diagnosis (years)) ytitle (Disease-free survival prob.) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
  
 sts graph if subtype ==4, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Disease-free survival TNBC) xtitle (Time from diagnosis (years)) ytitle (Disease-free survival prob.) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(7) ring(0)) ylabel(, angle(horizontal) ) xsize(7) ysize(6) 
 
**TIME TO RECURRENCE
*Event include only recurrence, but death without recurrence is censored 
   ** This does give different picture from overall or disease-free survival 
 ** Young patients are more likely to recur but less likely to die compared with old patients. 
   ** competing risk between dying from other diseases and breast cancer recurrence. 

stset year_recur, failure (recur)
** show # of patients included and # of events 
 tab ageethnic recur
 sts test (ageethnic)
 
*make KM graphs: use tmax option to cut time display (but all data are used), do not use "if condition" as the later remove patients with short follow-up. 
sts graph, by(ageethnic) tmax(15) risktable(, order(1 "Young White" 2 "Young Black" 3 "Older White" 4 "Older Black")) title(Time-to-recurrence) xtitle (Time from diagnosis (years)) ytitle (Probability of recurrence) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(1) position(11) ring(0)) ylabel(0(0.25)0.75, angle(horizontal) ) xsize(7) ysize(6) failure

stcox ib(frequent).ageethnic 
testparm i.ageethnic 
* young black vs young white 
       lincom 2.ageethnic - 1.ageethnic , hr
* young black vs old black 
       lincom 2.ageethnic - 4.ageethnic , hr 
estat phtest, detail

* multivariable analysis: different way of controlling molecular markers 
   stcox ib3.ageethnic i.ER i.PR i.her2  
testparm i.ageethnic 
stcox ib3.ageethnic ib1.subtype
testparm i.ageethnic 
stcox ib3.ageethnic ib1.subtype i.PR 
testparm i.ageethnic

* multivariable analysis: final model II     
  stcox ib3.ageethnic ib1.subtype i.grade i.stage1 i.CCI3
testparm i.ageethnic  

 * young black vs young white 
     lincom 2.ageethnic - 1.ageethnic , hr
* young black vs old black 
       lincom 2.ageethnic - 4.ageethnic , hr 
estat phtest, detail 


***[3D] TIME TO RECURRENCE/DEATH WITHOUT RECURRENCE, COMPETING RISK 
tab recur dead if ageethnic!=., m
tab recur dead if year_recur!=. & ageethnic!=., m 
tab typeRec1 ageethnic if ageethnic!=., m col 
tab recur dead if RFS !=. & year_RFS!=., m
gen outcome = 0 if RFS !=. & year_RFS!=.
replace outcome = 1 if RFS !=. & year_RFS!=. & recur==1
replace outcome = 2 if RFS !=. & year_RFS!=. & outcome!=1 & dead==1
label define outcome 0 "Event-free" 1 "Relapse" 2 "Dead w/o relapse"
label value outcome outcome 
tab outcome ageethnic, col nokey 
stset year_RFS, failure(outcome = 1) //  relapse/recurrence is the event of interest
stcrreg ib3.ageethnic , compete(outcome = 2) // Dead without relapse is the competing events

testparm i.ageethnic
* young black vs young white 
  lincom 2.ageethnic - 1.ageethnic , hr 
* young black vs old black 
    lincom 2.ageethnic - 4.ageethnic , hr 


stcurve, cif at(ageethnic = (1 2 3 4)) range(0 15) title(Recurrence) xtitle (Years from diagnosis) ytitle(Cumulative incidence) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(2) position(6) ring(1)) ylabel(0(0.1)0.4, angle(horizontal) ) xsize(4) ysize(5) name(cox_ci1, replace)

.  
stcrreg ib3.ageethnic ib1.subtype i.grade i.stage1 i.CCI3, compete(outcome = 2) // dead without relapse is the ecompeting event

testparm i.ageethnic 
testparm i.CCI3

testparm i.ageethnic  
  * young black vs young white 
    lincom 2.ageethnic - 1.ageethnic , hr 
 * young black vs old black 
    lincom 2.ageethnic - 4.ageethnic , hr 

stcrreg ib3.ageethnic ib1.subtype i.grade i.stage1, compete(outcome = 2) // relapse/recurrence is the event of interest, this is the model we used 

testparm i.ageethnic  
  * young black vs young white 
    lincom 2.ageethnic - 1.ageethnic , hr 
 * young black vs old black 
    lincom 2.ageethnic - 4.ageethnic , hr 

*cif is cumulative incidence function
stcurve, cif at(ageethnic = (1 2 3 4)) range(0 15) title(Recurrence) xtitle (Years from diagnosis) ytitle(Cumulative incidence)legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(2) position(6) ring(1)) ylabel(0(0.1)0.4, angle(horizontal) ) xsize(4) ysize(5) 

***not quite sure why reversing it is relevant***
stset year_RFS, failure(outcome = 2) //  Dead without relapse is the event of interest
stcrreg ib3.ageethnic , compete(outcome = 1) // relapse/recurrence is the competing event
testparm i.ageethnic
* young black vs young white 
.         lincom 2.ageethnic - 1.ageethnic , hr 
* young black vs old black 
.         lincom 2.ageethnic - 4.ageethnic , hr 
stcurve, cif at(ageethnic = (1 2 3 4)) range(0 15) title(Dead witout recurrence) xtitle (Years from diagnosis) ytitle(Cumulative incidence) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(2) position(6) ring(1)) ylabel(0(0.1)0.4, angle(horizontal) ) xsize(4) ysize(5) name(cox_ci2, replace)


graph combine cox_ci1 cox_ci2

stcrreg ib3.ageethnic ib1.subtype i.grade i.stage1 i.CCI3, compete(outcome = 1) // Dead without relapse is the competing event 
testparm i.ageethnic 
  * young black vs young white 
    lincom 2.ageethnic - 1.ageethnic , hr 
 * young black vs old black 
    lincom 2.ageethnic - 4.ageethnic , hr 

 stcurve, cif at(ageethnic = (1 2 3 4)) range(0 15) title(Dead witout recurrence) xtitle (Years from diagnosis) ytitle(Cumulative incidence) legend(order(1 "White <=40" 2 "Black <=40" 3 "White >=55" 4 "Black >=55") cols(2) position(6) ring(1)) ylabel(0(0.1)0.4, angle(horizontal) ) xsize(4) ysize(5) 
 

 
=


*****************
*TABLE 1
ssc install table1_mc
*label variables as I want them to be displayed in table
label variable clinicalstage1234 "Clinical Stage"
label define clinicalstage_labels 0 "DCIS"
label values clinicalstage1234 clinicalstage_labels 
label variable nodegrp "Lymph Node Status"
label variable subtype "Receptor Status"
label variable grade "Tumor Grade"
label variable stage1 "Stage"
*create table, variables listed are the rows 
table1_mc, by(ageethnic) vars(Age contn %4.1f \ tumor_size contn %4.1f \ nodegrp cat %4.1f \ subtype cat %4.1f \  stage1 cat %4.1f \ grade cat %4.0f) nospace onecol



