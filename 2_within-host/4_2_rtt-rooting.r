require(ape)

trefolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/4_5_Raxml/RAxML_bestTree*")


for (file in trefolder){
  tre <- read.tree(file)
  filename <- strsplit(basename(file), "bestTree.") [[1]][2]
  print(filename)
  
  tip.dates <- sapply(tre$tip.label, function(x){strsplit(x, "_")[[1]][2]})
  
  rtd <- rtt(tre, as.numeric(tip.dates))
  plot(rtd)
  rtd$edge.length <- rep(1, times=length(rtd$edge.length))
  
  plot(rtd)
  #write.tree(rtd, file=paste0("~/PycharmProjects/hiv-withinhost/4_5_Raxml/guide_trees/", filename))
  break
}


tre2 <- read.tree("~/PycharmProjects/hiv-withinhost/4_5_Raxml/RAxML_bestTree.51322.tree")
tip.dates2 <- sapply(tre2$tip.label, function(x){strsplit(x, "_")[[1]][2]})

rtd2 <- rtt(tre2, as.numeric(tip.dates2))
plot(rtd2)
rtd2$edge.length <- rep(1, times=length(rtd2$edge.length))
plot(rtd2)