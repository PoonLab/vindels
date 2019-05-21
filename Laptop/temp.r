
genetic.dists <- read.csv("~/vindels/Laptop/genetic-dists.csv", row.names=1)

# GENETIC DISTANCES PLOTS --------------

# qcquire the order of means 
new.df <- split(genetic.dists, genetic.dists$subtype)
means <- sapply(new.df, function (x) median(x$GD) )
means <- means[order(means)]
ordered <- names(means)


gd.order <- genetic.dists
# this applies the properly ordered factor to the subtype column 
gd.order$subtype <- factor(gd.order$subtype, levels=ordered)
# this sorts the entire data frame by the properly ordered 'levels' in the subtype column
gd.order <- gd.order[order(gd.order$subtype),]

new.order <- split(gd.order, gd.order$subtype)

# A) Standard plot 
par(mar=c(5,5,1,1))
plot(gd.order$subtype, gd.order$GD, xlab="Group M Clade", ylab="Genetic Distance (Cherries)", cex.axis=1.2, cex.lab=1.5)

# perform a wilcoxon test on this 
for (idx in 1:6){
  test <- wilcox.test(new.order[[idx]]$GD, new.order[[idx+1]]$GD)
  if (test$p.value < 0.05){
    arrows(idx, 0.15, idx+1, 0.15, length=0, lwd=3)
    text(idx+0.5, 0.16, labels="*", cex=1.5)
  }
}

# A.2) Log plot 
gd.order2 <- genetic.dists[genetic.dists$GD != 0,]
gd.order2$logged <- log10(gd.order2$GD)
gd.order2 <- gd.order2[gd.order2$logged>-4,]

# create the order of medians
new.df <- split(gd.order2, gd.order2$subtype)
means <- sapply(new.df, function (x) median(x$logged) )
means <- means[order(means)]
ordered <- names(means)

gd.order2$subtype <- factor(gd.order2$subtype, levels=ordered)
gd.order2 <- gd.order2[order(gd.order2$subtype),]

par(mar=c(5,5,1,1))
plot(gd.order2$subtype, gd.order2$logged, xlab="Group M Clade", ylab="log(Genetic Distance)", cex.axis=1.2, cex.lab=1.5)


# perform a wilcoxon test on this 
for (idx in 1:6){
  test <- wilcox.test(new.df[[idx]]$logged, new.df[[idx+1]]$logged)
  if (test$p.value < (0.05/6)){
    arrows(idx+0.1, -3.4, idx+0.9, -3.4, length=0, lwd=2.5)
    text(idx+0.5, -3.5, labels="*", cex=2)
  }
  print(test)
}

# B) Only GDs under 0.05 
gd.05 <- genetic.dists[which(genetic.dists$GD <= 0.05),]
par(mar=c(5,5,1,1))
plot(gd.05$subtype, gd.05$GD, xlab="Group M Clade", ylab="Genetic Distance (Cherries)", cex.axis=1.2, cex.lab=1.5)


# C) Staggered density plot of the seven subtypes  
require(ggplot2)
require(ggridges)
gd.15 <- genetic.dists[which(genetic.dists$GD <= 0.15),]
ggplot(gd.15, aes(x=GD, y=subtype, group=subtype)) + 
  geom_density_ridges(colour="white", fill="blue3", scale=1, bandwidth=0.002) + 
  labs(x="Genetic Distance", y="Subtype") + 
  theme(axis.title.x=element_text(size=15,margin=margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=15,margin=margin(t = 15, r = 0, b = 0, l = 0)),
        #strip.text.x = element_text(size=16),
        axis.text=element_text(size=13))


big.df <- big.df[which(big.df$dates>1960),]
new.df <- split(big.df[,2:3],big.df[,1])