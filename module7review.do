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
xtmixed bweight birthord initage || momid: birthord, cov(independent)




