
require(ape)
tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/6Trees", full.names=TRUE)

dists <- data.frame()
for (infile in tfolder){
  tre <- read.tree(infile)
  
  tre <- multi2di(tre)
  
  name <- strsplit(basename(infile), "\\.")[[1]][1]
  
  ntips <- Ntip(tre)
  numtips <- tabulate(tre$edge[which(tre$edge[,2] <= ntips),1])
  is.cherry <- which(numtips == 2)
  
  for (node in is.cherry){
    idx.cherry <- which(tre$edge[,1] == node)
    
  }
  
  tiplens <- tre$edge.length[which(tre$edge[,2] <= Ntip(tre))]
  
  dists <- rbind(dists, data.frame(subtype=rep(name, length(tiplens)), GD=tiplens))
}

plot(dists$GD,dists$subtype, )
