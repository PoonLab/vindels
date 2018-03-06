require(ape)
setwd("~/PycharmProjects/hiv-evolution-master/9_indels")
tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/8_1_Trees", full.names=TRUE)
vfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/3RegionSequences/VRegions2", full.names=TRUE)


for (i in 1:length(tfolder)){
  tre <- read.tree(tfolder[i])
  csv <- read.csv(vfolder[i], header=FALSE)
  
  names(csv) <- c('accno', 'VR1', 'VR2','VR3', 'VR4', 'VR5')
  
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

  # to add additional columns, use '$' operator
  df$tip1.label <- tre$tip.label[df$tip1]
  df$tip2.label <- tre$tip.label[df$tip2]
  df$tip1.len <- tre$edge.length[df$edge1]
  df$tip2.len <- tre$edge.length[df$edge2]
  
  
  indels <- df[,c(6:9)]

  for (x in 1:nrow(indels)){
    idxA <- match(indels$tip1.label[x], csv$accno)
    idxB <- match(indels$tip2.label[x], csv$accno)
    
    for (t in 2:ncol(csv)){
      Avr <- as.character(csv[idxA,t])
      Bvr <- as.character(csv[idxB,t])
      
      diff <- paste0("VR",as.character(t-1),".indel")
      nt <- paste0("VR",as.character(t-1),".nt")
      
      Alength <- nchar(Avr)
      Blength <- nchar(Bvr)
      
      indels[x,diff] <- Alength != Blength
      indels[x,nt] <- abs(Alength - Blength)
      
    }

  }
  name <- strsplit(tfolder[i], "/")[[1]][7]
  filename <- paste0(strsplit(name, "\\.")[[1]][1], ".csv" )
  
  write.csv(indels, filename)

}