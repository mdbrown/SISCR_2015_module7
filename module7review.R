### Module 7

library(lme4)
?glmer

library(geepack)
?geeglm

summ.fxn <- function(x) {
	nobs = length(x[!is.na(x)])
	mean = mean(x[!is.na(x)])
	sd = sd(x[!is.na(x)])
	min = min(x[!is.na(x)])
	max = max(x[!is.na(x)])

	return(c(nobs=round(nobs,0), mean=mean, sd=sd, min=round(min,0),max=round(max,0)))
}


### Motivating Example: birthweight plots

load('GAbabies.RData')
attach(data)

par(bty="n", mar=c(5, 5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
plot(jitter(birthord), bweight, col="lightsteelblue3", xlab="Birth order", ylab="Birth weight (grams)", ylim=c(0,5000))
lines(smooth.spline(birthord, bweight, df=4), lwd=5, col="dodgerblue4")

par(bty="n", mar=c(5, 5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
boxplot(split(bweight, birthord), lwd=2, pch=19, xlab="Birth order", ylab="Birth weight (grams)", ylim=c(0,5000))

## Birth weight individual lines

samp <- sample(unique(momid), 30)

par(bty="n", mar=c(5, 5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
plot(birthord, bweight, pch=19, col="lightsteelblue3", xlab="Birth order", ylab="Birth weight (grams)", ylim=c(0,5000), type="n", main="Observed data")
for(i in samp){
	lines(birthord[momid==i], bweight[momid==i], lwd=2, col="dodgerblue4")}
lines(lowess(birthord, bweight), lwd=5)

plot(birthord, bweight, pch=19, col="lightsteelblue3", xlab="Birth order", ylab="Birth weight (grams)", ylim=c(0,5000), type="n", main="Individual fits")
for(i in samp){
	m <- lm(bweight ~ birthord, subset=momid==i)
	lines(birthord[momid==i], fitted(m), lwd=2, col="dodgerblue4")}
lines(lowess(birthord, bweight), lwd=5)

## GEE
# Fit a linear model with independence correlation
summary(geeglm(bweight ~ birthord + initage, id=momid, data=data, corstr="independence"))

# Fit a linear model with exchangeable correlation
summary(geeglm(bweight ~ birthord + initage, id=momid, data=data, corstr="exchangeable"))

## Mixed Models
# Fit a linear model with random intercepts
summary(lmer(bweight ~ (1 | momid) + birthord + initage, data=data))

# Fit a linear model with random intercepts and slopes
summary(lmer(bweight ~ (birthord | momid) + birthord + initage, data=data))
