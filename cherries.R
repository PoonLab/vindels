
setwd("~/PycharmProjects/hiv-evolution-master/8_Dated_Trees")
tfolder <- list.files(path=".", full.names=TRUE)


tre <- read.tree(text='(((A:1.1,B:1.2):0.1,(C:1.3,D:1.4):0.2):0.3,E:1.5):0.4;')


for (i in 1:length(tfolder)){
  #tre <- read.tree(paste("~/PycharmProjects/hiv-evolution-master/8_Dated_Trees",tfolder[5],sep="/"))
  tre <- read.tree(tfolder[5])
  
  n <- length(tre$tip.label)  # Ntip(tre)
  
  # number of tips per internal node
  numtips <- tabulate(tre$edge[,1][tre$edge[,2] <= n])
  
  #determines which nodes contain cherries (returns vector with their integer positions)
  is.cherry <- sapply(numtips, function(x) x==2)
  #vect = c()
  #for (x in 1:length(numtips)){
  #  if (numtips[x] == 2){
  #    vect <- c(vect, x)
  #  }
  #}
  
  # get row indices for edges from each cherry
  
  
  # construct data frame where each row corresponds to a cherry
  m <- sapply(which(is.cherry), function(i) { 
    edge.idx <- tre$edge[,1]==i  # select rows in edge matrix with parent node index
    c(i,   # index of the parent node
      which(edge.idx),
      t(     # transpose so we can concatenate this vector to i
        tre$edge[edge.idx, 2]    # column of tip indices
        )
      )
    })
  df <- data.frame(node.index=m[1,], edge1=m[2,], edge2=m[3,], tip1=m[4,], tip2=m[5,])

  # to add additional columns, use '$' operator
  df$tip1.label <- tre$tip.label[df$tip1]
  df$tip2.label <- tre$tip.label[df$tip2]
  df$tip1.len <- tre$edge.length[df$edge1]
  df$tip2.len <- tre$edge.length[df$edge2]
  
  
  
  
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