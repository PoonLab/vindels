tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/8_Dated_Trees/", full.names=FALSE)
setwd("~/PycharmProjects/hiv-evolution-master/8_Dated_Trees")




for (i in 1:length(tfolder)){
  tre <- read.tree(paste("~/PycharmProjects/hiv-evolution-master/8_Dated_Trees",tfolder[5],sep="/"))
  
  n <- length(tre$tip.label)
  
  numtips <- tabulate(tre$edge[,1][tre$edge[,2] <= n])
  
  #determines which nodes contain cherries (returns vector with their integer positions)
  vect = c()
  for (x in 1:length(numtips)){
    if (numtips[x] == 2){
      vect <- c(vect, x)
    }
  }
  
  
  #makes a vector containing the positions of cherry tips
  lens <- c()
  labels <- c()
  for (i in vect){
    print(i)
    for (n in 1:length(tre$edge[,1])){
      if (i == tre$edge[,1][n]){
        lens <- c(lens, tre$edge.length[n])
        labels <- c(labels, tre$tip.label[tre$edge[n,2]])
      }
    }
  }
  
  

  
  
  
  print(numtips)
  
}