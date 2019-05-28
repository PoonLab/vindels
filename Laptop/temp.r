# CHERRY GENETIC DISTANCES

# this script is to plot different breakdowns of cherry genetic distances 


genetic.dists <- read.csv("~/vindels/Pipeline_2_Within/genetic-dists.csv", row.names=1)
filtered.indels <- read.csv("~/vindels/Pipeline_2_Within/filtered-indels.csv", row.names = 1)

genetic.dists <- cbind(genetic.dists, filtered.indels[match(genetic.dists$tip1.label, filtered.indels$tip1.label),])

# GENETIC DISTANCES PLOTS --------------

# qcquire the order of means 
new.df <- split(genetic.dists, genetic.dists$subtype)
means <- sapply(new.df, function (x) median(x$total.length) )
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
plot(gd.order$subtype, gd.order$total.length, xlab="Group M Clade", ylab="Genetic Distance (Cherries)", cex.axis=1.2, cex.lab=1.5)

# perform a wilcoxon test on this 
for (idx in 1:6){
  test <- wilcox.test(new.order[[idx]]$total.length, new.order[[idx+1]]$total.length)
  if (test$p.value < 0.05){
    arrows(idx, 0.15, idx+1, 0.15, length=0, lwd=3)
    text(idx+0.5, 0.16, labels="*", cex=1.5)
  }
}

# A.2) Log plot 
gd.order2 <- genetic.dists[genetic.dists$total.length != 0,]
gd.order2$logged <- log10(gd.order2$total.length)
gd.order2 <- gd.order2[gd.order2$logged>-4,]

# create the order of medians
new.df <- split(gd.order2, gd.order2$subtype)
means <- sapply(new.df, function (x) median(x$logged) )
means <- means[order(means)]
ordered <- names(means)

gd.order2$subtype <- factor(gd.order2$subtype, levels=ordered)
gd.order2 <- gd.order2[order(gd.order2$subtype),]

par(mar=c(5,5,1,1))
plot(gd.order2$subtype, gd.order2$logged, xlab="Group M Clade", ylab="log10(Cherry Genetic Distances)", cex.axis=1.2, cex.lab=1.5, varwidth=T)


# perform a wilcoxon test on this 
for (idx in 1:6){
  test <- wilcox.test(new.df[[idx]]$logged, new.df[[idx+1]]$logged)
  if (test$p.value < (0.05/6)){
    arrows(idx+0.1, -2.9, idx+0.9, -2.9, length=0, lwd=2.5)
    text(idx+0.5, -2.95, labels="*", cex=2)
  }
  print(test)
}

# B) Only GDs under 0.05 
gd.05 <- genetic.dists[which(genetic.dists$total.length <= 0.05),]
par(mar=c(5,5,1,1))
plot(gd.05$subtype, gd.05$total.length, xlab="Group M Clade", ylab="Genetic Distance (Cherries)", cex.axis=1.2, cex.lab=1.5)


# C) Staggered density plot of the seven subtypes  
require(ggplot2)
require(ggridges)
gd.15 <- genetic.dists[which(genetic.dists$total.length <= 0.15),]
ggplot(gd.15, aes(x=total.length, y=subtype, group=subtype)) + 
  geom_density_ridges(colour="white", fill="blue3", scale=1, bandwidth=0.002) + 
  labs(x="Genetic Distance", y="Subtype") + 
  theme(axis.title.x=element_text(size=15,margin=margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=15,margin=margin(t = 15, r = 0, b = 0, l = 0)),
        #strip.text.x = element_text(size=16),
        axis.text=element_text(size=13))


big.df <- big.df[which(big.df$dates>1960),]
new.df <- split(big.df[,2:3],big.df[,1])