
require(ape)
tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/6Trees", full.names=TRUE)

dists <- data.frame()
for (infile in tfolder){
  tre <- read.tree(infile)
  
  name <- strsplit(basename(infile), "\\.")[[1]][1]
  tiplens <- tre$edge.length[which(tre$edge[,2] <= Ntip(tre))]
  
  dists <- rbind(dists, data.frame(subtype=rep(name, length(tiplens)), GD=tiplens))
}

boxplot(dists$GD,dists$subtype)
