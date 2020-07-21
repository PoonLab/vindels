require(ape)


tre <- read.tree("~/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/rooted_trees/16362.tree")

par(mar=c(2,2,2,2))
plot(tre, show.tip.label=F)
add.scale.bar(x=0.08,y=-0.5,lwd=1.5)
