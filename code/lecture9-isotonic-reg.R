## Isotonic regression - examples

library(cir)
library(ggplot2)
library(plyr)

## Plot of non-monotone proportions and isotonic regresion
d <- seq(1:6)
p <- c(1/10, 2/10, 19/100, 3/20, 12/16, 8/10)
p
p1 <- oldPAVA(p)

plot(p ~ d, cex=1.25, pch=16, col="black", ylim=c(0,1), xlab="Dose", ylab="Probability of success")
points(p1 ~ d, type='o', cex=1.25, pch=4, lty=2, lwd=2, col='darkgray')
text(1.5,.9,"Observed proportion",adj=0, cex=1.25)
points(1.25,.9, cex=1.25, pch=16, col="black")
text(1.5,.8,"Isotonic regression",adj=0, cex=1.25)
points(1.25,.8, cex=1.25, pch=4, lwd=2, col="darkgray")
segments(1, .8, 1.5, .8, lty=2, lwd=2, col="darkgray")


## Bar plot of the number of dose allocations in the case study
trial <- data.frame(outcome=rep(c("#Failures", "#Successes"), each=5),
			dose=rep(c("25", "27.5", "30", "32.5", "35"), 2),
			obs=c(1, 2, 4, 1, 0, 1, 3, 13, 24, 5)
)

trial_sorted <- arrange(trial, desc(outcome), dose)
trial_cumsum <- ddply(trial_sorted, "dose", transform, label_obs=cumsum(obs))

ggplot(data=trial_cumsum, aes(x=dose, y=obs, fill=outcome)) + 
geom_bar(stat="identity") + 
geom_text(aes(y=label_obs, label=obs), vjust=1.6, color="white", size=3.5) +
scale_fill_brewer(palette="Paired") +
scale_y_continuous(name="Number assigned", limits=c(0, 25)) +
ggtitle("Distribution of dose assignments") +
theme_minimal()


## Plot of proportions and isotonic regression in the case study
d <- c(25, 27.5, 30, 32.5, 35)
p <- c(1/2, 3/5, 13/17, 24/25, 5/5)
w <- c(2, 5, 17, 25, 5)
p1 <- oldPAVA(p, wt=w)

plot(p ~ d, cex=1.25, pch=16, col="black", xaxt='n', ylim=c(0,1), xlab="Dose", ylab="Probability of success")
points(p1 ~ d, type='o', cex=1.25, pch=4, lty=2, lwd=2, col='darkgray')
axis(1, at=d, labels=d)
text(25.5,.9,"Observed proportion",adj=0, cex=1.25)
points(25,.9, cex=1.25, pch=16, col="black")
text(25.5,.8,"Isotonic regression",adj=0, cex=1.25)
points(25,.8, cex=1.25, pch=4, lwd=2, col="darkgray")
segments(24.5, .8, 25.5, .8, lty=2, lwd=2, col="darkgray")


## 95% CI for MEV90 using bootstrap
d <- c(25, 27.5, 30, 32.5, 35)
n <- c(2, 5, 17, 25, 5)
p <- c(1/2, 3/5, 13/17, 24/25, 5/5)
B <- 10000
MEV90.est <- c()
isoton.reg <- matrix(0, ncol=5, nrow=B)
MEV90.est[1] <- (.9 - 13/17)/(24/25 - 13/17)*(32.5 - 30) + 30
isoton.reg[1,] <- p
for (i in 2:B){ 
	x <- c()
	for (j in 1:length(n)) x[j] <- rbinom(1, n[j], p[j])
	isoton.reg[i,] <- oldPAVA(x/n, wt=n)

	if (((isoton.reg[i,1] <= .9) & (isoton.reg[i,2] > .9))|(isoton.reg[i,1] > .9)) { mu.est <- 25 }
	else {mu.est <- max(d[isoton.reg[i,] <= 0.9])}

	index <- which(d==mu.est)

	if (isoton.reg[i,index+1] - isoton.reg[i,index]==0) {MEV90.est[i] <- d[index+1]}
	else {MEV90.est[i] <- (.9 - isoton.reg[i,index])/(isoton.reg[i,index+1] - isoton.reg[i,index])*(d[index+1] - mu.est) + mu.est}
}


j <- round(.025*B)
k <- round(.975*B)
est.sorted <- sort(MEV90.est)
c(est.sorted[j], est.sorted[k])

bootstrap <- as.data.frame(MEV90.est)
bootstrap$MEV90.est

ggplot(data=bootstrap, aes(bootstrap$MEV90.est)) +
geom_histogram(breaks=seq(25, 37.5, by = .5), col="black", fill="blue") +
geom_vline(xintercept=est.sorted[j], lty=2) + 
geom_vline(xintercept=est.sorted[k], lty=2) + 
geom_text(aes(x=est.sorted[j], label="2.5th percentile", y=2000), angle=90, vjust=-1, text=element_text(size=5)) +
geom_text(aes(x=est.sorted[k], label="97.5th percentile", y=2000), angle=90, vjust=-1, text=element_text(size=5)) +
xlab("Dose")  +
ggtitle("Bootstrap distribution of estimated MEV90") +
theme_minimal()


