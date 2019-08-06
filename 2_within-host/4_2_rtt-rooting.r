require(ape)

trefolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/4_5_Raxml/RAxML_bestTree*")


for (file in trefolder){
  tre <- read.tree(file)
  filename <- strsplit(basename(file), "bestTree.") [[1]][2]
  print(filename)
  
  tip.dates <- sapply(tre$tip.label, function(x){strsplit(x, "_")[[1]][2]})
  
  rtd <- rtt(tre, as.numeric(tip.dates))
  
  write.tree(rtd, file=paste0("~/PycharmProjects/hiv-withinhost/4_5_Raxml/rooted_trees/", filename))
}

rtdfolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/4_5_Raxml/rooted_trees/*.tree")


for (file in rtdfolder){
  filename <- basename(file)
  rescaled <- read.tree(file)
  rescaled$edge.length <- rep(1, times=length(rescaled$edge.length))
  write.tree(rescaled, file=paste0("~/PycharmProjects/hiv-withinhost/4_5_Raxml/guide_trees/", filename))
}

