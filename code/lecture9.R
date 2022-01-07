library(lattice)
library(drc) 
library(ggplot2)

#### I. Modeling Pr(Success) as a function of Dose ####
d <- c(25, 27.5, 30, 32.5, 35)
n <- c(2, 5, 17, 25, 5)
x <- c(1, 3, 13, 24, 5)
p <- x/n
tox <- as.data.frame(cbind(d, p))

resp <- c(); dose <- c()
for (i in 1:length(d)){
	resp <- c(resp, rep(1, x[i]), rep(0,n[i]-x[i]))
	dose <- c(dose, rep(d[i],n[i]))
}

data.tox <- as.data.frame(cbind(resp, dose))

## Logistic regression model: 
logistic <- glm(resp ~ dose, data=data.tox, family=binomial())
l.curve <- predict(logistic, data.frame(dose=seq(25, 35, by=.25)), type="response")
pred <- as.data.frame(cbind(l.curve, seq(25,35, by=.25)))

## 95% CIs for alpha and beta 
confint(logistic) ## for alpha and beta

## 95% CIs for P(d) 
vc <- vcov(logistic)
CI <- matrix(0,nrow=length(d),ncol=4)
for (i in 1:length(d)){
	est <- logistic$coef[1] + d[i]*logistic$coef[2]
	V <- vc[1,1] + (d[i])^2*vc[2,2] + 2*d[i]*vc[1,2]
	lwr <- est - 1.96*sqrt(V)
	upr <- est + 1.96*sqrt(V)
	CI[i,] <- c(d[i], 1/(1+exp(-lwr)), 1/(1+exp(-upr)), 1/(1+exp(-est)))
}
as.data.frame(CI)

## 95% CI for LD50
alpha.est <- logistic$coef[1]
beta.est <- logistic$coef[2]
LD50.est <- -alpha.est/beta.est 
V50 <- 1/beta.est^2*(vc[1,1] + LD50.est^2*vc[2,2] + 2*LD50.est*vc[1,2])
lwr <- LD50.est - 1.96*sqrt(V50)
upr <- LD50.est + 1.96*sqrt(V50)
LD50 <- c(LD50.est, lwr, upr)

## 95% CI for LD90
g <- log(.9/.1) 
LD90.est <- (g - alpha.est)/beta.est
V90 <- 1/beta.est^2*(vc[1,1] + LD90.est^2*vc[2,2] + 2*LD90.est*vc[1,2])
lwr <- LD90.est - 1.96*sqrt(V90)
upr <- LD90.est + 1.96*sqrt(V90)
LD90 <- c(LD90.est, lwr, upr)

plot(p ~ d, data=data.tox, cex=1.25, pch=16, col="black", ylim=c(0,1), xlab="Dose", ylab="Probability of success")
points(l.curve ~ V2, data=pred, type='l', lty=1, lwd=4, col="darkgray") 
text(30,.4,"Observed proportion",adj=0, cex=1.25)
points(29.5,.4, cex=1.25, pch=16, col="black")
text(30,.35,"Fitted logistic model",adj=0, cex=1.25)
segments(29, .35, 29.5, .35, lty=1, lwd=4, col="darkgray")
text(30,.3,"Estimated MEV90",adj=0, cex=1.25)
points(29.5, .3, cex=1.25, pch=15, col="darkgreen")
text(30,.25,"95% CI for MEV90",adj=0, cex=1.25)
segments(29, .25, 29.5, .25, lty=1, lwd=4, col="darkgreen")
abline(h=.9, lty=2, col="darkgreen")
points(LD90[1],.9,pch=15,cex=1.25,col="darkgreen")
segments(LD90[2], .9, LD90[3], .9, lty=1, lwd=4, col="darkgreen")


## 4-parameter Emax model
d <- seq(0, 100, by=1)

emax <- function(d, e0, emax, ed50, r){
	f <- c()
	for (i in 1:length(d)) f[i] <- e0 + emax*(d[i])^r/((ed50)^r + (d[i])^r)
	return(f)
}

f0 <- emax(d, e0=0, emax=1, ed50=20, r=.5)
f1 <- emax(d, e0=0, emax=1, ed50=20, r=1)
f2 <- emax(d, e0=0, emax=1, ed50=20, r=2.5)
f3 <- emax(d, e0=0, emax=1, ed50=20, r=5)

plot(f1 ~ d, type='l', col="darkgray", lwd=4, lty=1, xlab="d", ylab="Effect", ylim=c(0,1))
points(f2 ~ d, type='l', lty=1, lwd=4, col="blue")
points(f3 ~ d, type='l', lty=1, lwd=4, col="brown")
points(f0 ~ d, type='l', lty=1, lwd=4, col="darkgreen")
title(expression(paste("4-parameter Emax model:  ",E[0]," + ",E[max]*d^r/(ED[50]^r+d^r))))
text(65, .4, "r=0.5", adj=0)
segments(58, .4, 63, .4, lty=1, lwd=4, col="darkgreen")
text(65, .35, "r=1", adj=0)
segments(58, .35, 63, .35, lty=1, lwd=4, col="darkgray")
text(65, .3, "r=2.5", adj=0)
segments(58, .3, 63, .3, lty=1, lwd=4, col="blue")
text(65, .25, "r=5", adj=0)
segments(58, .25, 63, .25, lty=1, lwd=4, col="brown")
text(65, .15, expression(paste(E[0],"=0; ",E[max],"=1; ",ED[50],"=20")), adj=0)


