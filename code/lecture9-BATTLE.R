## Bar plots visualizing BATTLE trial data

library(ggplot2)

df2 <- data.frame(drug=rep(c("Erlotinib", "Vandetanib", "Erlotinib + Bexarotene", "Sorafenib"), each=5),
			BM=rep(c("EGFR", "KRAS/BRAF", "VEGF", "RXR/CyclinD1", "None"), 4),
			n=c(17, 7, 25, 1, 8, 27, 3, 16, 0, 6, 20, 3, 3, 1, 9, 23, 14, 39, 4, 18),
			resp=c(6, 1, 10, 0, 3, 11, 0, 6, 0, 0, 11, 1, 0, 1, 5, 9, 11, 25, 1, 11)
)

df2$prop <- round(df2$resp/df2$n, digits=2)
df3 <- df2[df2$BM!="RXR/CyclinD1",]

ggplot(data=df2, aes(x=drug, y=n, fill=BM)) + geom_bar(stat="identity") + 
scale_fill_brewer(palette="Paired")+
ggtitle("Treatment allocation numbers in the BATTLE trial") +
theme_minimal()

ggplot(data=df3, aes(x=drug, y=prop, fill=BM)) + 
  geom_bar(stat="identity", , position=position_dodge()) +
scale_fill_brewer(palette="Paired")+
geom_text(aes(label=prop), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5) +
ggtitle("Observed 8-week DCR's in the BATTLE trial") +
theme_minimal()

m1 <- c(58, 52, 36, 98)
m1/sum(m1)

m2 <- c(87, 27, 83, 6, 41)
m2/sum(m2)