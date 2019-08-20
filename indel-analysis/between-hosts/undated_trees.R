require(ape)
require(stringr)
require(ggplot2)

tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/6Trees", full.names=TRUE)




# This code reads phylogenetic trees and variable loop sequences in csv format. 
branch.lengths <- data.frame()

for (i in 1:7){
  tre <- read.tree(tfolder[i])
 
  
  #for output 
  name <- strsplit(tfolder[i], "/")[[1]][7]
  subtype <- strsplit(name, "_CR")[[1]][1]
  filename <- paste0(subtype,"+.csv" )
  
  tre <- multi2di(tre) 
  #naming the csv 
  #names(csv) <- c('accno', 'VR1', 'VR2','VR3', 'VR4', 'VR5')
  
  #counting tips
  n <- Ntip(tre)
  
  # number of tips per internal node
  # count the number of instances that first column (node) corresponds to a second column number which is <= n (meaning it is a tip)
  numtips <- tabulate(tre$edge[,1][tre$edge[,2] <= n])
  
  #determines which nodes contain cherries (returns vector with their integer positions)
  is.cherry <- sapply(numtips, function(d) d==2)
  
  
  # construct data frame where each row corresponds to a cherry
  m <- sapply(which(is.cherry), function(a) { #will input the numbers of nodes containing 2 tips?
    edge.idx <- tre$edge[,1]==a  # FINDS THE EDGES (row #) corresponding with the parent node ; ap: select rows in edge matrix with parent node index
    
    c(a,   # index of the parent node
      which(edge.idx),
      t(     # transpose so we can concatenate this vector to i
        tre$edge[edge.idx, 2]    # column of tip indices
      )
    )
  })
  df <- data.frame(node.index=m[1,], edge1=m[2,], edge2=m[3,], tip1=m[4,], tip2=m[5,])
  
  
  df$tip1.label <- tre$tip.label[df$tip1]
  df$tip2.label <- tre$tip.label[df$tip2]
  df$tip1.len <- tre$edge.length[df$edge1]
  df$tip2.len <- tre$edge.length[df$edge2]
  
  
  indels <- df[,c(6:9)]
  indels$total.length <- indels$tip1.len + indels$tip2.len
  
  branch.lengths <- rbind(branch.lengths, data.frame(subtype=rep(subtype,nrow(indels)), length=indels$total.length))
  
}

branch.len2 <- split(branch.lengths$length, branch.lengths$subtype)
par()
boxplot(branch.len2)
