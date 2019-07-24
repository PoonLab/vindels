require(ape)

trefolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/4_5_Raxml/RAxML_bestTree*")


for (file in trefolder){
  tre <- read.tree(file)
  filename <- strsplit(basename(file), "bestTree.") [[1]][2]
  print(filename)
  
  if (filename == "VN_Data.tree"){
    tip.dates <- sapply(tre$tip.label, function(x){strsplit(x, "_")[[1]][3]})
  }else{
    tip.dates <- sapply(tre$tip.label, function(x){strsplit(x, "_")[[1]][2]})
  }
  
  rtd <- rtt(tre, as.numeric(tip.dates))
  write.tree(rtd, file=paste0("~/PycharmProjects/hiv-withinhost/4_5_Raxml/rooted_trees/", filename))
}


