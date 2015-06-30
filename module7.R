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

	return(c(nobs=round(nobs,0), mean=mean, sd=sd, min=min, max=max))
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






### Case Study: Postnatal Depression
rm(list=ls())
load('depress.RData')

attach(depress)

summ.fxn(dep0[group=="placebo"])
summ.fxn(dep1[group=="placebo"])
summ.fxn(dep2[group=="placebo"])
summ.fxn(dep3[group=="placebo"])
summ.fxn(dep4[group=="placebo"])
summ.fxn(dep5[group=="placebo"])
summ.fxn(dep6[group=="placebo"])

summ.fxn(dep0[group=="estrogen"])
summ.fxn(dep1[group=="estrogen"])
summ.fxn(dep2[group=="estrogen"])
summ.fxn(dep3[group=="estrogen"])
summ.fxn(dep4[group=="estrogen"])
summ.fxn(dep5[group=="estrogen"])
summ.fxn(dep6[group=="estrogen"])

# pairwise correlation plots:
pairs(~dep0+dep1+dep2+dep3+dep4+dep5+dep6)

subp=(!is.na(dep6) & !is.na(dep0) & group=="placebo")
t.test(dep0[subp],dep6[subp], paired=T)

sube=(!is.na(dep6) & !is.na(dep0) & group=="estrogen")
t.test(dep0[sube],dep6[sube], paired=T)

## calculate correlation across measurements (also var-cov matrix)
## note: Stata "corr" command uses participants measured at all times; 
##	    code below only requires pairwise-complete observations
ntime <- 7
rmat <- with(depress, cbind(dep0, dep1, dep2, dep3, dep4, dep5, dep6))
	
vmat <- nmat <- covmatj <- covmatk <- matrix(0, nrow=ntime, ncol=ntime)
rownames(vmat) <- colnames(vmat) <- 1:ntime
rownames(nmat) <- colnames(nmat) <- 1:ntime

for(j in 1:ntime){
    for(k in 1:ntime){
        nmat[j,k] <- sum(!is.na(rmat[,j]*rmat[,k]))
 	  vmat[j,k] <- cov(rmat[,j],rmat[,k],use="pairwise.complete.obs")
	 if(nmat[j,k]>0) covmatj[j,k] <- cov(rmat[,j][!is.na(rmat[,j]*rmat[,k])],rmat[,j][!is.na(rmat[,j]*rmat[,k])],use="pairwise.complete.obs")
	 if(nmat[j,k]>0) covmatk[j,k] <- cov(rmat[,k][!is.na(rmat[,j]*rmat[,k])],rmat[,k][!is.na(rmat[,j]*rmat[,k])],use="pairwise.complete.obs")
	}}

cmat <- vmat/sqrt(covmatj)/sqrt(covmatk)

round(vmat,3)
round(cmat, 2)
nmat
		


# change data to long form
detach(depress)
depresslong <- reshape(data=depress, varying = c("dep0", "dep1", "dep2", "dep3", "dep4", "dep5", "dep6"), 
	v.names="dep", timevar="visit", idvar="subj", times = 0:6, 
	new.row.names = 1:427, direction="long")
depresslsort <- depresslong[order(depresslong$subj),]

attach(depresslsort)

# plot depression scores vs visit, by treatment group
par(bty="n", mar=c(5, 5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
plot(jitter(visit[group=="placebo"]), dep[group=="placebo"], 
	pch=21, col="blue", bg="blue", xlab="Visit", ylab="Depression Score", ylim=c(0,30))
points(jitter(visit[group=="estrogen"]), dep[group=="estrogen"], 
	pch=21, col="darkgreen", bg="darkgreen", xlab="Visit", ylab="Depression Score", ylim=c(0,30))
lines(lowess(visit[!is.na(dep) & group=="placebo"], 
	dep[!is.na(dep) & group=="placebo"]), lwd=5, col="blue")
lines(lowess(visit[!is.na(dep) & group=="estrogen"], 
	dep[!is.na(dep) & group=="estrogen"]), lwd=5, col="darkgreen")
legend(3,30, legend=c("Placebo", "Estrogen"), pch=c(21,21), col=c("blue", "darkgreen"), pt.bg=c("blue", "darkgreen"), lty=c(1,1), bty="n", cex=1.5, lwd=2)


boxplot(dep~visit, data=depresslsort)


# plot depression scores vs visit, by individual
par(bty="n", mar=c(5, 5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
plot(visit, dep, xlab="Visit", ylab="Depression Score", ylim=c(0,30), type="n", main="")
for(i in 1:61){
	lines(visit[subj==i], dep[subj==i], lwd=2, col="dodgerblue4")}


# GEE
m1 <- geeglm(dep~group+visit, id=subj, data=depresslsort, corstr="independence")
m2 <- geeglm(dep~group+visit, id=subj, data=depresslsort, corstr="ar1")

summary(m1)
summary(m2, correlation=T)

genCor <- function(maxClustSize, alpha, corstr){
	if(corstr=="independence"){#indep
		return(as.matrix(Diagonal(x=rep(1,maxClustSize))))	
	}
	if(corstr=="exchangeable"){#exch
			return(alpha^(!outer(1:maxClustSize,1:maxClustSize,"==")))
	}
	if(corstr=="ar1"){#ar1
			return(alpha^abs(outer(1:maxClustSize,1:maxClustSize,"-")))
	}
	if (corstr == "unstructured") {
		m <- as.matrix(Diagonal(x=rep(1,maxClustSize)))
		m[lower.tri(m)] <- alpha
		m <- t(m)
		m[lower.tri(m)] <- alpha
	}
	stop("Unrecognized correlation structure, should be one of 'independence', 'exchangeable', 'unstructured', or 'ar1'")
}

genCor(maxClustSize=max(summary(m2)$clusz),alpha=summary(m2)$corr[1,1], corstr="ar1")


# categorical visit term
m3 <- geeglm(dep~group+factor(visit), id=subj, data=depresslsort, corstr="ar1")
m0 <- geeglm(dep~group, id=subj, data=depresslsort, corstr="ar1")
summary(m3)
anova(m0,m3)


# interaction between visit and treatment group
m4 <- geeglm(dep~group*visit, id=subj, data=depresslsort, corstr="ar1")
summary(m4)

# calculate linear combination
m4$coef["visit"]+m4$coef["groupestrogen:visit"]
# and its standard error
sqrt(summary(m4)$cov.scaled[3,3] + summary(m4)$cov.scaled[4,4] + 2*summary(m4)$cov.scaled[3,4]) 

### alternate method using function in 'car' package:
library(car)
depresslsort$groupvisit=with(depresslsort,(as.numeric(group)-1)*as.numeric(visit))
m4=geeglm(dep~group+visit+groupvisit, id=subj, data=depresslsort, corstr="ar1")
deltaMethod(coef(m4),"visit+groupvisit",vcov.=summary(m4)$cov.scaled)

### evaluate whether dropout depends on depression levels
library(plyr)
library(ggplot2)

depresslong <- ddply(depresslong, .(subj), mutate, 
              maxtime = max(visit[!is.na(dep)]))

ggplot(aes(x=visit, y=dep), data=depresslong) +
  geom_point(aes(x=visit, y=dep)) +
  geom_smooth(aes(x=visit, y=dep, group=maxtime), se=FALSE) +
  theme_bw()

### fix mixed models which have weaker assumptions about missing data
summary(lmer(dep ~ (1|subj) + visit + factor(group), 
             data=depresslong))
summary(lmer(dep ~ (visit | subj) + visit*factor(group), 
             data=depresslong))




### Case Study: ICHS

rm(list=ls())
load('ichs.RData')

## Exploratory data analysis

# Categorize by year of age
ichs <- within(ichs, visit <- floor((age+36)/12))

# Cross-tabulate vitamin A deficiency and respiratory infection
with(ichs, table(visit,xerop,infection))
    		
# Selected individual trajectories
ichsIDs <- c(121013,121113,121114,121140,121414,122214,122315,152150,141517)

par(bty="n", mar=c(4.5, 4.5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, mfrow=c(3,3))
for(i in 1:length(ichsIDs)) {
	with(ichs, plot(age+36, infection, type="n", xlab="Age (years)", ylab="", ylim=c(0, 3), axes=FALSE))
	title(main=paste("subject",ichsIDs[i],sep=""))
	axis(1, at=c(12,24,36,48,60,72,84), labels=1:7)
	axis(2, at=c(0,1,2,3), labels=c("N", "Y", "N", "Y"), las=1)	
	mtext("Xerop",side=2, las=1, at=0.5)
	mtext("Infection",side=2, las=1, at=2.5)
	with(ichs, lines(age[id==ichsIDs[i]]+36, infection[id==ichsIDs[i]]+2, col="firebrick4", lwd=3))
	with(ichs, lines(age[id==ichsIDs[i]]+36, xerop[id==ichsIDs[i]], col="skyblue4", lwd=3))
}

# Probabilities of infection and vitamin A deficiency by age
infectionbyyear.xerop <- with(ichs, split(infection[xerop==1],as.factor(age[xerop==1])))
meaninfbyyear.xerop <- lapply(infectionbyyear.xerop,mean)
time.xerop <- with(ichs, sort(unique(age[xerop==1])))

infectionbyyear.noxerop <- with(ichs, split(infection[xerop==0],as.factor(age[xerop==0])))
meaninfbyyear.noxerop <- lapply(infectionbyyear.noxerop,mean)
time.noxerop <- with(ichs, sort(unique(age[xerop==0])))

par(bty="n", mar=c(4.5, 4.5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5)
plot(time.noxerop+36, meaninfbyyear.noxerop, type="n", xlab="Age (years)", ylab="Percent with Respiratory Infection",  xlim=c(0,86),ylim=c(0, 1), axes=FALSE)
axis(1, at=c(0,12,24,36,48,60,72,84), labels=0:7)
axis(2, c(0,0.2,0.4,0.6,0.8,1),las=1)
points(time.noxerop+36, meaninfbyyear.noxerop, col="firebrick4", pch=1, cex=2)	
points(time.xerop+36, meaninfbyyear.xerop, col="skyblue4", pch=16,cex=2)
legend(48,0.75, legend=c("No Xerop", "Xerop"), pch=c(1,16), col=c("firebrick4", "skyblue4"), bty="n", cex=1.5)

infectionbyvisit.xerop <- with(ichs, split(infection[xerop==1],as.factor(visit[xerop==1])))
meaninfbyvisit.xerop <- lapply(infectionbyvisit.xerop,mean)
visit.xerop <- with(ichs, sort(unique(visit[xerop==1])))

infectionbyvisit.noxerop <- with(ichs, split(infection[xerop==0],as.factor(visit[xerop==0])))
meaninfbyvisit.noxerop <- lapply(infectionbyvisit.noxerop,mean)
visit.noxerop <- with(ichs, sort(unique(visit[xerop==0])))

par(bty="n", mar=c(4.5, 4.5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5)
plot(visit.noxerop, meaninfbyvisit.noxerop, type="n", xlab="Age (years)", ylab="Percent with Respiratory Infection",  ylim=c(0, 1), axes=FALSE)
axis(1, at=0:7, labels=0:7)
axis(2, c(0,0.2,0.4,0.6,0.8,1),las=1)
points(visit.noxerop, meaninfbyvisit.noxerop, col="firebrick4", pch=1, cex=2)	
points(visit.xerop, meaninfbyvisit.xerop, col="skyblue4", pch=16,cex=2)
legend(4,0.75, legend=c("No Xerop", "Xerop"), pch=c(1,16), col=c("firebrick4", "skyblue4"), bty="n", cex=1.5)

# Logistic regression model (ignoring correlation)
model <- glm(infection ~ xerop + age + gender + hfora + cost + sint, 
                data=ichs, family="binomial")
summary(model)
exp(model$coef[2])
exp(model$coef[2]+qnorm(0.975)*summary(model)$coef[2,2])
exp(model$coef[2]-qnorm(0.975)*summary(model)$coef[2,2])

# Variance, covariance, correlations
ichswide <- reshape(data=with(ichs, data.frame(infection, id, time)), v.names="infection", timevar="time", idvar="id", direction="wide")

ntime <- with(ichs, length(unique(time)))

rmat <- with(ichswide, cbind(infection.1, infection.2, infection.3, 
					infection.4, infection.5, infection.6))
	
vmat <- nmat <- covmatj <- covmatk <- matrix(0, nrow=ntime, ncol=ntime)
rownames(vmat) <- colnames(vmat) <- 1:ntime
rownames(nmat) <- colnames(nmat) <- 1:ntime

for(j in 1:ntime){
    for(k in 1:ntime){
        nmat[j,k] <- sum(!is.na(rmat[,j]*rmat[,k]))
 	  vmat[j,k] <- cov(rmat[,j],rmat[,k],use="pairwise.complete.obs")
	 if(nmat[j,k]>0) covmatj[j,k] <- cov(rmat[,j][!is.na(rmat[,j]*rmat[,k])],rmat[,j][!is.na(rmat[,j]*rmat[,k])],use="pairwise.complete.obs")
	 if(nmat[j,k]>0) covmatk[j,k] <- cov(rmat[,k][!is.na(rmat[,j]*rmat[,k])],rmat[,k][!is.na(rmat[,j]*rmat[,k])],use="pairwise.complete.obs")
	}}

cmat <- vmat/sqrt(covmatj)/sqrt(covmatk)

round(vmat,3)
round(cmat, 2)
nmat
		
round(colMeans(rmat,na.rm=T),2) 

## Generalized estimating equations

rm(list=ls())
load('ichs.RData')

m1 <- geeglm(infection ~ xerop + age + gender + hfora + cost + sint, 
             id=id, data=ichs, family="binomial", corstr="independence")
             
m2 <- geeglm(infection ~ xerop + age + gender + hfora + cost + sint, 
             id=id, data=ichs, family="binomial", corstr="exchangeable")
             
m3 <- geeglm(infection ~ xerop + age + gender + hfora + cost + sint,
              id=id, data=ichs, family="binomial", corstr="ar1")

exp(m1$coef[2])
exp(m1$coef[2]+qnorm(0.975)*summary(m1)$coef[2,2])
exp(m1$coef[2]-qnorm(0.975)*summary(m1)$coef[2,2])

exp(m2$coef[2])
exp(m2$coef[2]+qnorm(0.975)*summary(m2)$coef[2,2])
exp(m2$coef[2]-qnorm(0.975)*summary(m2)$coef[2,2])

summary(m2)$geese$correlation[1]

exp(m3$coef[2])
exp(m3$coef[2]+qnorm(0.975)*summary(m3)$coef[2,2])
exp(m3$coef[2]-qnorm(0.975)*summary(m3)$coef[2,2])

summary(m3)$geese$correlation[1]
summary(m3)$geese$correlation[1]^2
summary(m3)$geese$correlation[1]^3
summary(m3)$geese$correlation[1]^4
summary(m3)$geese$correlation[1]^5

## Generalized linear mixed-effects models

rm(list=ls())
load('ichs.RData')


m_ri <- glmer(infection ~ (1 | id) + factor(xerop) + age + factor(gender) + hfora + cost + sint,
	 family=binomial, data=ichs)

expit <- function(x){exp(x)/(1+exp(x))}

expit(fixef(m_ri)[1])
expit(fixef(m_ri)[1]-1.96*VarCorr(m_ri)$id)
expit(fixef(m_ri)[1]+1.96*VarCorr(m_ri)$id)



### Case Study: Carpal Tunnel Syndrome
rm(list=ls())
load('cts.RData')

attach(cts)

table(treatassign,surgical)
table(treatassign[!is.na(ctsaqf4)],surgical[!is.na(ctsaqf4)])

summ.fxn(ctsaqf0[treatassign==0])
summ.fxn(ctsaqf1[treatassign==0])
summ.fxn(ctsaqf2[treatassign==0])
summ.fxn(ctsaqf3[treatassign==0])
summ.fxn(ctsaqf4[treatassign==0])

summ.fxn(ctsaqf0[treatassign==1])
summ.fxn(ctsaqf1[treatassign==1])
summ.fxn(ctsaqf2[treatassign==1])
summ.fxn(ctsaqf3[treatassign==1])
summ.fxn(ctsaqf4[treatassign==1])


# change to multiple records per person ("long" format)
detach(cts)
cts <- within(cts, ctsaqfbase<-ctsaqf0)
ctslong <- reshape(data=cts, varying = list(c("ctsaqf0", "ctsaqf1", "ctsaqf2", "ctsaqf3", "ctsaqf4"), 
	c("ctsaqs0", "ctsaqs1", "ctsaqs2", "ctsaqs3", "ctsaqs4"), 
	c("surgreported0","surgreported1","surgreported2","surgreported3","surgreported4")),
	v.names=c("ctsaqf","ctsaqs","surgreported"), timevar="visit", idvar="ID", times = 0:4, 
	new.row.names = 1:580, direction="long")
ctslsort <- ctslong[order(ctslong$ID),]

attach(ctslsort)

# plot outcomes vs visit, for a sample of individuals
ctsIDs0 <- unique(ID[ID<=13062 & ID!=13009 & treatassign==0])
ctsIDs1 <- unique(ID[ID<=13101 & ID!=13009 & treatassign==1])

par(bty="n", mar=c(4.5, 4.5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, mfrow=c(4,5))
for(i in 1:length(ctsIDs0)) {
	with(ctslsort, plot(visit, ctsaqf, type="n", xlab="Visit", ylab="", ylim=c(1,5), main=""))
	title(main=paste("ID",ctsIDs0[i],sep=""))
	with(ctslsort, lines(visit[ID==ctsIDs0[i]], ctsaqf[ID==ctsIDs0[i]], lwd=2))
	with(ctslsort, points(visit[ID==ctsIDs0[i]], surgreported[ID==ctsIDs0[i]], pch=21, col="red", bg="red"))
}

par(bty="n", mar=c(4.5, 4.5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, mfrow=c(4,5))
for(i in 1:length(ctsIDs1)) {
	with(ctslsort, plot(visit, ctsaqf, type="n", xlab="Visit", ylab="", ylim=c(1,5), main=""))
	title(main=paste("ID",ctsIDs1[i],sep=""))
	with(ctslsort, lines(visit[ID==ctsIDs1[i]], ctsaqf[ID==ctsIDs1[i]], lwd=2))
	with(ctslsort, points(visit[ID==ctsIDs1[i]], surgreported[ID==ctsIDs1[i]], pch=21, col="red", bg="red"))
}


# plot outcomes vs visit, by treatment group
ctsaqfmean0 <- c(mean(ctsaqf[treatassign==0 & visit==0],na.rm=T),
			mean(ctsaqf[treatassign==0 & visit==1],na.rm=T),
			mean(ctsaqf[treatassign==0 & visit==2],na.rm=T),
			mean(ctsaqf[treatassign==0 & visit==3],na.rm=T),
			mean(ctsaqf[treatassign==0 & visit==4],na.rm=T))

ctsaqfmean1 <- c(mean(ctsaqf[treatassign==1 & visit==0],na.rm=T),
			mean(ctsaqf[treatassign==1 & visit==1],na.rm=T),
			mean(ctsaqf[treatassign==1 & visit==2],na.rm=T),
			mean(ctsaqf[treatassign==1 & visit==3],na.rm=T),
			mean(ctsaqf[treatassign==1 & visit==4],na.rm=T))

par(bty="n", mar=c(5, 5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
plot(jitter(visit[treatassign==0]), ctsaqf[treatassign==0], 
	pch=21, col="blue", bg="blue", xlab="Visit", ylab="CTSAQF", ylim=c(1,5))
points(jitter(visit[treatassign==1]), ctsaqf[treatassign==1], 
	pch=21, col="darkgreen", bg="darkgreen", xlab="Visit", ylab="CTSAQF", ylim=c(1,5))
lines(unique(visit), ctsaqfmean0, lwd=5, col="blue")
lines(unique(visit), ctsaqfmean1, lwd=5, col="darkgreen")
legend(0,5, legend=c("Non-surgical", "Surgical"), pch=c(21,21), col=c("blue", "darkgreen"), pt.bg=c("blue", "darkgreen"), lty=c(1,1), bty="n", cex=1.2, lwd=2)

detach(ctslsort)

# back to wide data
# compute correlations by treatment group
ntime <- 5
sub0 <- (treatassign==0 & !is.na(ctsaqf0) & !is.na(ctsaqf1) & !is.na(ctsaqf2) 
		& !is.na(ctsaqf3) & !is.na(ctsaqf4))
rmat0 <- with(cts, cbind(ctsaqf0[sub0], ctsaqf1[sub0], ctsaqf2[sub0], ctsaqf3[sub0], ctsaqf4[sub0]))

sub1 <- (treatassign==1 & !is.na(ctsaqf0) & !is.na(ctsaqf1) & !is.na(ctsaqf2) 
		& !is.na(ctsaqf3) & !is.na(ctsaqf4))
rmat1 <- with(cts, cbind(ctsaqf0[sub1], ctsaqf1[sub1], ctsaqf2[sub1], ctsaqf3[sub1], ctsaqf4[sub1]))

# non-surgical group	
rmat=rmat0
vmat <- nmat <- covmatj <- covmatk <- matrix(0, nrow=ntime, ncol=ntime)
rownames(vmat) <- colnames(vmat) <- 1:ntime
rownames(nmat) <- colnames(nmat) <- 1:ntime

for(j in 1:ntime){
    for(k in 1:ntime){
        nmat[j,k] <- sum(!is.na(rmat[,j]*rmat[,k]))
 	  vmat[j,k] <- cov(rmat[,j],rmat[,k],use="pairwise.complete.obs")
	 if(nmat[j,k]>0) covmatj[j,k] <- cov(rmat[,j][!is.na(rmat[,j]*rmat[,k])],rmat[,j][!is.na(rmat[,j]*rmat[,k])],use="pairwise.complete.obs")
	 if(nmat[j,k]>0) covmatk[j,k] <- cov(rmat[,k][!is.na(rmat[,j]*rmat[,k])],rmat[,k][!is.na(rmat[,j]*rmat[,k])],use="pairwise.complete.obs")
	}}

cmat <- vmat/sqrt(covmatj)/sqrt(covmatk)

round(vmat,3)
round(cmat, 2)
nmat
		

# surgical group
rmat=rmat1
vmat <- nmat <- covmatj <- covmatk <- matrix(0, nrow=ntime, ncol=ntime)
rownames(vmat) <- colnames(vmat) <- 1:ntime
rownames(nmat) <- colnames(nmat) <- 1:ntime

for(j in 1:ntime){
    for(k in 1:ntime){
        nmat[j,k] <- sum(!is.na(rmat[,j]*rmat[,k]))
 	  vmat[j,k] <- cov(rmat[,j],rmat[,k],use="pairwise.complete.obs")
	 if(nmat[j,k]>0) covmatj[j,k] <- cov(rmat[,j][!is.na(rmat[,j]*rmat[,k])],rmat[,j][!is.na(rmat[,j]*rmat[,k])],use="pairwise.complete.obs")
	 if(nmat[j,k]>0) covmatk[j,k] <- cov(rmat[,k][!is.na(rmat[,j]*rmat[,k])],rmat[,k][!is.na(rmat[,j]*rmat[,k])],use="pairwise.complete.obs")
	}}

cmat <- vmat/sqrt(covmatj)/sqrt(covmatk)

round(vmat,3)
round(cmat, 2)
nmat


# generate change variables
attach(cts)
cts <- within(cts, change1 <- ctsaqf1-ctsaqf0)
cts <- within(cts, change2 <- ctsaqf2-ctsaqf0)
cts <- within(cts, change3 <- ctsaqf3-ctsaqf0)
cts <- within(cts, change4 <- ctsaqf4-ctsaqf0)

# post only
with(cts, t.test(ctsaqf1~treatassign))
with(cts, t.test(ctsaqf2~treatassign))
with(cts, t.test(ctsaqf3~treatassign))
with(cts, t.test(ctsaqf4~treatassign))

# change only
with(cts, t.test(change1~treatassign))
with(cts, t.test(change2~treatassign))
with(cts, t.test(change3~treatassign))
with(cts, t.test(change4~treatassign))

# post adjusted for baseline
summary(lm(ctsaqf1~treatassign+ctsaqf0, data=cts))
summary(lm(ctsaqf2~treatassign+ctsaqf0, data=cts))
summary(lm(ctsaqf3~treatassign+ctsaqf0, data=cts))
summary(lm(ctsaqf4~treatassign+ctsaqf0, data=cts))

# primary analysis, 12-mo timepoint
summary(lm(ctsaqf4~treatassign+ctsaqf0+factor(idgroup), data=cts))

# longitudinal ITT -- naive
ctslsortpost <- ctslsort[ctslsort$visit!=0,]
summary(lm(ctsaqf~treatassign+ctsaqfbase+factor(idgroup), data=ctslsortpost))

# GEE
summary(geeglm(ctsaqf~treatassign+ctsaqfbase+factor(idgroup),
              id=ID, data=ctslsortpost, corstr="independence"))

# Mixed Models
summary(lmer(ctsaqf~(1|ID)+treatassign+ctsaqfbase+factor(idgroup), data=ctslsortpost))


# As-treated Analyses
# generate variables
cts <- within(cts, surgby3 <- (surgical==1))
cts <- within(cts, surgby6 <- (surgical==1 | surgical==2))
cts <- within(cts, surgby9 <- (surgical==1 | surgical==2 | surgical==3))

with(cts, table(surgby3,treatassign))
with(cts, table(surgby9,treatassign))

ctslsort <- within(ctslsort, surgby3 <- (surgical==1))
ctslsort <- within(ctslsort, surgby6 <- (surgical==1 | surgical==2))
ctslsort <- within(ctslsort, surgby9 <- (surgical==1 | surgical==2 | surgical==3))

attach(ctslsort) 

# plot outcomes vs visit, by treatment group
ctsaqfmean0surg3 <- c(mean(ctsaqf[surgby3==0 & visit==0],na.rm=T),
			mean(ctsaqf[surgby3==0 & visit==1],na.rm=T),
			mean(ctsaqf[surgby3==0 & visit==2],na.rm=T),
			mean(ctsaqf[surgby3==0 & visit==3],na.rm=T),
			mean(ctsaqf[surgby3==0 & visit==4],na.rm=T))

ctsaqfmean1surg3 <- c(mean(ctsaqf[surgby3==1 & visit==0],na.rm=T),
			mean(ctsaqf[surgby3==1 & visit==1],na.rm=T),
			mean(ctsaqf[surgby3==1 & visit==2],na.rm=T),
			mean(ctsaqf[surgby3==1 & visit==3],na.rm=T),
			mean(ctsaqf[surgby3==1 & visit==4],na.rm=T))

par(bty="n", mar=c(5, 5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
plot(jitter(visit[surgby3==0]), ctsaqf[surgby3==0], 
	pch=21, col="blue", bg="blue", xlab="Visit", ylab="CTSAQF", ylim=c(1,5))
points(jitter(visit[surgby3==1]), ctsaqf[surgby3==1], 
	pch=21, col="darkgreen", bg="darkgreen", xlab="Visit", ylab="CTSAQF", ylim=c(1,5))
lines(unique(visit), ctsaqfmean0surg3, lwd=5, col="blue")
lines(unique(visit), ctsaqfmean1surg3, lwd=5, col="darkgreen")
legend(0,5, legend=c("No surgery by 3 months", "Surgery by 3 months"), pch=c(21,21), col=c("blue", "darkgreen"), pt.bg=c("blue", "darkgreen"), lty=c(1,1), bty="n", cex=1.2, lwd=2)


ctsaqfmean0surg9 <- c(mean(ctsaqf[surgby9==0 & visit==0],na.rm=T),
			mean(ctsaqf[surgby9==0 & visit==1],na.rm=T),
			mean(ctsaqf[surgby9==0 & visit==2],na.rm=T),
			mean(ctsaqf[surgby9==0 & visit==3],na.rm=T),
			mean(ctsaqf[surgby9==0 & visit==4],na.rm=T))

ctsaqfmean1surg9 <- c(mean(ctsaqf[surgby9==1 & visit==0],na.rm=T),
			mean(ctsaqf[surgby9==1 & visit==1],na.rm=T),
			mean(ctsaqf[surgby9==1 & visit==2],na.rm=T),
			mean(ctsaqf[surgby9==1 & visit==3],na.rm=T),
			mean(ctsaqf[surgby9==1 & visit==4],na.rm=T))

par(bty="n", mar=c(5, 5, 2, 2) + 0.1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
plot(jitter(visit[surgby9==0]), ctsaqf[surgby9==0], 
	pch=21, col="blue", bg="blue", xlab="Visit", ylab="CTSAQF", ylim=c(1,5))
points(jitter(visit[surgby9==1]), ctsaqf[surgby9==1], 
	pch=21, col="darkgreen", bg="darkgreen", xlab="Visit", ylab="CTSAQF", ylim=c(1,5))
lines(unique(visit), ctsaqfmean0surg9, lwd=5, col="blue")
lines(unique(visit), ctsaqfmean1surg9, lwd=5, col="darkgreen")
legend(0,5, legend=c("No surgery by 9 months", "Surgery by 9 months"), pch=c(21,21), col=c("blue", "darkgreen"), pt.bg=c("blue", "darkgreen"), lty=c(1,1), bty="n", cex=1.2, lwd=2)

detach(ctslsort)


# Mixed Models
ctslsortpost <- ctslsort[ctslsort$visit!=0,]
summary(lmer(ctsaqf~(1|ID)+surgby3+ctsaqfbase+factor(idgroup)+age, data=ctslsortpost))
summary(lmer(ctsaqf~(1|ID)+surgby9+ctsaqfbase+factor(idgroup)+age, data=ctslsortpost))







### Case Study: Guatemala Immunization Campaign
rm(list=ls())
load('guatemala.RData')

attach(guatemala)

# Exploratory Data Analysis
comm_size=NULL
for(i in 1:length(unique(cluster))) {
	comm_size[cluster==unique(cluster)[i]] <- length(cluster[cluster==unique(cluster)[i]])
}
summ.fxn(comm_size)

fam_size=NULL
for(i in 1:length(unique(mom))) {
	fam_size[mom==unique(mom)[i]] <- length(mom[mom==unique(mom)[i]])
}
summ.fxn(fam_size)

table(immun,kid2p)

* Mixed Effects Models
m1 <- glmer(immun ~ (1 | cluster) + (1 | mom) + kid2p + indNoSpa + indSpa + momEdPri + momEdSec
		+ husEdPri + husEdSec + husEdDK + rural + pcInd81, data=guatemala, family=binomial)
summary(m1)
exp(summary(m1)$coef[2,1])
exp(summary(m1)$coef[2,1]+qnorm(.025)*summary(m1)$coef[2,2])	
exp(summary(m1)$coef[2,1]+qnorm(.975)*summary(m1)$coef[2,2])

G3 <- as.numeric(VarCorr(m1)$cluster)
G2 <- as.numeric(VarCorr(m1)$mom)
# ICC for random intercepts model:
as.numeric((G3)/(G3+G2+pi^2/3))
as.numeric(G3+G2)/(G3+G2+pi^2/3)


m2 <- glmer(immun ~ (kid2p | cluster) + (1 | mom) + kid2p + rural + pcInd81, data=guatemala, family=binomial)
summary(m2)

# fit reduced random intercepts model and test need for random slopes with LR test
m0 <- glmer(immun ~ (1 | cluster) + (1 | mom) + kid2p + rural + pcInd81, data=guatemala, family=binomial)
anova(m0,m2)






