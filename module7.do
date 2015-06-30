* Module 7

* Longitudinal Analysis Review: GA babies data
use "GAbabies.dta", clear

* Exploratory Data Analysis
twoway (scatter bweight birthord, jitter(10) mfcolor(none) mlcolor(midblue)) ///
       (lowess bweight birthord, lc(dknavy) lw(2)), legend(off) ///
       xscale(range(0.75 5.25)) xtitle("Birth order") ytitle("Birth weight (grams)")

graph box bweight, by(birthord, cols(5)) ytitle("Birth weight (grams)")

* Birth weight individual lines
xtline bweight if momid<=2383, i(momid) t(birthord) overlay legend(off) ///
     xlab(1(1)5) ylab(0(1000)5000) xtitle("Birth order") ytitle("Birth weight (grams)") ///
     addplot(lowess bweight birthord, lw(2))

* Fitted lines
xi: regress bweight i.momid*birthord if momid<=2383
predict p
xtline p if momid<=2383, overlay t(birthord) i(momid) legend(off) xlab(1(1)5) ylab(0(1000)5000) ///
   xtitle("Birth order") ytitle("Birth weight (grams)") addplot(lowess bweight birthord, lw(2))



* GEE
* Declare the dataset to be "panel" data, grouped by momid
* with time variable birthord
xtset momid birthord

* Fit a linear model with independence correlation
xtgee bweight birthord initage, i(momid) corr(ind) robust

* Fit a linear model with exchangeable correlation
xtgee bweight birthord initage, i(momid) corr(exc) robust


* Mixed Models
* Declare the dataset to be "panel" data, grouped by momid
* with time variable birthord
xtset momid birthord

* Fit a linear model with random intercepts
xtmixed bweight birthord initage || momid:

* Fit a linear model with random intercepts and slopes
xtmixed bweight birthord initage || momid: birthord






* Case Study: Postnatal Depression
use "depress.dta", clear

* Exploratory Data Analysis
sort group
by group: summarize dep0 dep1 dep2 dep3 dep4 dep5 dep6

graph matrix dep0 dep1 dep2 dep3 dep4 dep5 dep6, half

sort group
by group: ttest dep0==dep6

corr dep0 dep1 dep2 dep3 dep4 dep5 dep6

* change data to long form
reshape long dep, i(subj) j(visit)

graph twoway (lowess dep visit, mean) (scatter dep visit, jitter(10))
graph box dep, over(visit)


xtline dep, i(subj) t(visit) overlay legend(off)
     xlab(0(1)6) xtitle("Visit") ytitle("Depression Score")


* GEE
xtset subj visit
xtgee dep visit group, i(subj) corr(ind) robust

xtgee dep visit group, i(subj) corr(ar1) robust
estat wcorr


xi: xtgee dep group i.visit, i(subj) corr(ar1) robust
test _Ivisit_1 _Ivisit_2 _Ivisit_3 _Ivisit_4 _Ivisit_5 _Ivisit_6

gen group_visit=group*visit
xtgee dep visit group group_visit, i(subj) corr(ar1) robust
lincom visit + group_visit



* Case Study: ICHS
* Exploratory data analysis

use ichs.dta, replace

list id age time infection xerop gender hfora cost sint

* Categorize by year of age
gen visit=floor((age+36)/12)

* Cross-tabulate vitamin A deficiency and respiratory infection
bysort xerop: tab visit infection


* Declare the dataset to be "panel" data, grouped by id
* with time variable time 
xtset id time

* Individual trajectories
xtline infection if id<121141, legend(off)
xtline xerop if id<121141, legend(off)

* Probabilities of infection and vitamin A deficiency
bysort xerop: summarize infection if time==1
bysort xerop: summarize infection if time==2
bysort xerop: summarize infection if time==3
bysort xerop: summarize infection if time==4
bysort xerop: summarize infection if time==5
bysort xerop: summarize infection if time==6

gen ageyear=floor((age+36)/12)
bysort xerop: summarize infection if ageyear==0
bysort xerop: summarize infection if ageyear==1
bysort xerop: summarize infection if ageyear==2
bysort xerop: summarize infection if ageyear==3
bysort xerop: summarize infection if ageyear==4
bysort xerop: summarize infection if ageyear==5
bysort xerop: summarize infection if ageyear==6

* Covariance and correlation
drop age ageyear cost sint hfora
reshape wide infection xerop, i(id) j(time)
corr infection1 infection2 infection3 infection4 infection5 infection6
corr infection1 infection2 infection3 infection4 infection5 infection6, cov

* Generalized estimating equations

use ichs.dta, replace

* Declare the dataset to be "panel" data, grouped by id
* with time variable time
xtset id time

* Fit a GEE model with exchangeable working correlation
xtgee infection xerop age gender hfora cost sint, family(binomial) link(logit) corr(exch) robust
estat wcorr

* Generalized linear mixed-effects models

* Fit a model with random intercepts
help melogit
melogit infection i.xerop age i.gender hfora cost sint || id:

* Obtain predicted probabilities of infection, 
* setting the random effects to 0
margins i.xerop, predict(mu fixed)




* Case Study: Carpal Tunnel Syndrome
use "cts.dta", clear

*** generate baseline variable to be used when data are in long format
gen ctsaqfbase = ctsaqf0

***
*** some exploratory data analysis
***
  
summarize age
tab gender

tab idgroup 	
tab treatassign
tab treatassign surgical 
tab treatassign surgical if ctsaqf4!=.

bysort treatassign: summarize ctsaqf0 ctsaqf1 ctsaqf2 ctsaqf3 ctsaqf4 
bysort treatassign: cor ctsaqf0 ctsaqf1 ctsaqf2 ctsaqf3 ctsaqf4 

***
*** plot means by treatment group
***

*** change data from wide format to long format
reshape long surgreported ctsaqs ctsaqf, i(ID) j(visit)
save "cts_long.dta", replace



***
*** plot some individual series - one plot for treatassign=0 and one for treatassign==1
***

use "cts_long.dta", clear

graph twoway (lowess ctsaqf visit) (scatter ctsaqf visit) ///
      if (ID<=13062 & ID!=13009 & treatassign==0), by(ID)
list ID surgical if (ID<=13062 & ID!=13009 & treatassign==0)

graph twoway (lowess ctsaqf visit) (scatter ctsaqf visit) ///
      if (ID<=13101 & ID!=13009 & treatassign==1), by(ID)
list ID surgical if (ID<=13101 & ID!=13009 & treatassign==1)

***
*** plot some fitted lines
***

graph twoway (lfit ctsaqf visit) (scatter ctsaqf visit) ///
      if (ID<=13062 & ID!=13009 & treatassign==0), by(ID)

graph twoway (lfit ctsaqf visit) (scatter ctsaqf visit) ///
      if (ID<=13101 & ID!=13009 & treatassign==1), by(ID)


***
*** plot average outcome over time by treatment group
***


collapse (mean) ctsaqf, by(visit treatassign)

graph twoway (scatter ctsaqf visit if treatassign==0, msize(large) mfcolor(black)) (line ctsaqf visit if treatassign==0, lcolor(black) lpattern(solid)) ///
	(scatter ctsaqf visit if treatassign==1, msize(large)  mfcolor(white)) (line ctsaqf visit if treatassign==1, lcolor(black) lpattern(dash)), ///
	yscale(range(1 5)) ylabel(1(1)5) legend(order(2 "non-surgical" 4 "surgical"))



***
*** derive change scores
***

use "cts.dta", clear

gen change1 = ctsaqf1 - ctsaqf0 
gen change2 = ctsaqf2 - ctsaqf0 
gen change3 = ctsaqf3 - ctsaqf0 
gen change4 = ctsaqf4 - ctsaqf0 

***
*** [1] Analysis of post only
***
ttest ctsaqf1, by(treatassign) unequal   
ttest ctsaqf2, by(treatassign) unequal
ttest ctsaqf3, by(treatassign) unequal
ttest ctsaqf4, by(treatassign) unequal

***
*** [2] Analysis of change
***
ttest change1, by(treatassign) unequal    
ttest change2, by(treatassign) unequal    
ttest change3, by(treatassign) unequal    
ttest change4, by(treatassign) unequal    

***
*** [3] Analysis of post controlling for baseline
***
reg ctsaqf1 treatassign ctsaqf0
reg ctsaqf2 treatassign ctsaqf0
reg ctsaqf3 treatassign ctsaqf0
reg ctsaqf4 treatassign ctsaqf0

***
*** [4] Analysis of change controlling for baseline
***
reg change1 treatassign ctsaqf0
reg change2 treatassign ctsaqf0
reg change3 treatassign ctsaqf0
reg change4 treatassign ctsaqf0


***
*** PRIMARY ANALYSIS: post controlling for baseline, plus site
***
xi: reg ctsaqf4 treatassign ctsaqf0 i.idgroup




**************************************************************
* longitudinal regression analysis - ITT				 *
**************************************************************

*****
***** Vector of outcomes, with robust standard error estimates
*****  

xi: reg ctsaqf treatassign ctsaqfbase i.idgroup if visit!=0
xi: reg ctsaqf treatassign ctsaqfbase i.idgroup if visit!=0, cluster(ID)

*****
***** Generalized Estimating Equations
*****  

xtset ID visit

xi: xtgee ctsaqf treatassign ctsaqfbase i.idgroup if visit!=0, corr(independent) robust
xi: xtgee ctsaqf treatassign ctsaqfbase i.idgroup if visit!=0, corr(independent) robust rgf

xi: xtgee ctsaqf treatassign ctsaqfbase i.idgroup if visit!=0, corr(exchangeable) robust
xtcorr

*****
***** Linear Mixed Models Analysis
*****  

*** Random Intercepts, all post time points in outcome vector, adjusted for baseline

xi: xtmixed ctsaqf treatassign i.idgroup ctsaqfbase if visit!=0 ||  ID:  

 

 
**************************************************************
* as-treated analysis							 *
**************************************************************

*** generate as-treated surgery indicators
gen surgby3 = (surgical==1)
gen surgby6 = (surgical==1 | surgical==2)
gen surgby9 = (surgical==1 | surgical==2 | surgical==3)

save "cts_long.dta", replace

collapse (mean) surgby3 surgby6 surgby9 treatassign, by(ID)
tab surgby3 treatassign
tab surgby6 treatassign
tab surgby9 treatassign

use "cts_long.dta", clear


collapse (mean) ctsaqf, by(visit surgby3)

graph twoway (scatter ctsaqf visit if surgby3==0, msize(large) mfcolor(black)) (line ctsaqf visit if surgby3==0, lcolor(black) lpattern(solid)) ///
	(scatter ctsaqf visit if surgby3==1, msize(large)  mfcolor(white)) (line ctsaqf visit if surgby3==1, lcolor(black) lpattern(dash)), ///
	yscale(range(1 5)) ylabel(1(1)5) legend(order(2 "non-surgical" 4 "surgical"))


use "cts_long.dta", clear

collapse (mean) ctsaqf, by(visit surgby9)

graph twoway (scatter ctsaqf visit if surgby9==0, msize(large) mfcolor(black)) (line ctsaqf visit if surgby9==0, lcolor(black) lpattern(solid)) ///
	(scatter ctsaqf visit if surgby9==1, msize(large)  mfcolor(white)) (line ctsaqf visit if surgby9==1, lcolor(black) lpattern(dash)), ///
	yscale(range(1 5)) ylabel(1(1)5) legend(order(2 "non-surgical" 4 "surgical"))



***
*** SECONDARY ANALYSIS: LMM controlling for baseline, site, age
***
use "cts_long.dta", clear

xi: xtmixed ctsaqf surgby3 i.idgroup ctsaqfbase age if visit!=0 ||  ID:  
xi: xtgee ctsaqf surgby3 ctsaqfbase i.idgroup age if visit!=0, corr(independent) robust

xi: xtmixed ctsaqf surgby9 i.idgroup ctsaqfbase age if visit!=0 ||  ID:  
xi: xtgee ctsaqf surgby9 ctsaqfbase i.idgroup age if visit!=0, corr(independent) robust





* Case Study: Guatemalan Immunization Campaign
use "guatemala.dta", clear

* Exploratory Data Analysis
bysort cluster: gen comm_size=_N 
summ comm_size
bysort cluster mom: gen fam_size=_N 
summ fam_size

tab immun kid2p

* Mixed Effects Models

xtmelogit immun kid2p indNoSpa  indSpa momEdPri momEdSec husEdPri ///
	husEdSec husEdDK rural pcInd81 || cluster: || mom:
lincom kid2p, eform

xtmelogit immun kid2p rural pcInd81 || cluster: kid2p, cov(unstructured) || mom:
estimates store ri_rs
xtmelogit, or

xtmelogit immun kid2p rural pcInd81 || cluster: || mom:, or
estimates store ri
lrtest ri ri_rs




