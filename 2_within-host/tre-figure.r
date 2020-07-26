require(ape)


tre <- read.tree("~/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/rooted_trees/16362.tree")

par(mar=c(2,2,2,2))
plot(tre, show.tip.label=F, edge.width=1.8)
add.scale.bar(x=0.08,y=-0.5,lwd=2.0,cex=1.5)


xcr <- c(0.062,0.031,0.025,0.023, 0.037,0.049,0.025, 0.013) 
ycr <- c(27.5,26,37.3,79.5,  14,33,57.5, 101)
col <- c(rep("green",4), rep("red",4))
ty <- c(rep(25,4),rep(24,4))

points(xcr,ycr, pch=ty, bg=col,cex=1.5)
