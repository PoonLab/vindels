require(ape)
setwd("~/PycharmProjects/hiv-evolution-master/9_1_indels")
tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/8_1_Trees_cut", full.names=TRUE)
vfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/3RegionSequences/VRegions3", full.names=TRUE)

zeros <- data.frame()
len.diff <- data.frame()

for (i in 1:length(tfolder)){
  tre <- read.tree(tfolder[i])
  csv <- read.csv(vfolder[i], header=FALSE)
  
  #for output 
  name <- strsplit(tfolder[i], "/")[[1]][7]
  filename <- paste0(strsplit(name, "\\.")[[1]][1], ".csv" )
  
  #naming the csv 
  names(csv) <- c('accno', 'VR1', 'VR2','VR3', 'VR4', 'VR5')
  
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

  # to add additional columns, use '$' operator
  df$tip1.label <- tre$tip.label[df$tip1]
  df$tip2.label <- tre$tip.label[df$tip2]
  df$tip1.len <- tre$edge.length[df$edge1]
  df$tip2.len <- tre$edge.length[df$edge2]
  
  
  indels <- df[,c(6:9)]
  indels$total.length <- indels$tip1.len + indels$tip2.len
  
  
  filtered.indels <- indels
  # for (y in 1:nrow(indels)){
  #   value <- indels$total.length[y]
  #   if (value != 0){
  #      filtered.indels <- rbind(filtered.indels, data.frame(indels[y,]))
  #   }
  # }

  
  
  for (x in 1:nrow(filtered.indels)){
    idxA <- match(filtered.indels$tip1.label[x], csv$accno)
    idxB <- match(filtered.indels$tip2.label[x], csv$accno)
    
    for (t in 2:ncol(csv)){
      Avr <- as.character(csv[idxA,t])
      Bvr <- as.character(csv[idxB,t])
      
      bln <- paste0("VR",as.character(t-1),".indel")
      num <- paste0("VR",as.character(t-1),".nt")
      
      Alength <- nchar(Avr)
      Blength <- nchar(Bvr)
      
      # a.count <- str_count(q.data$string, "a")
      # t.count <- str_count(q.data$string, "t")
      # c.count <- str_count(q.data$string, "c")
      # g.count <- str_count(q.data$string, "g")
      *************************
      diff <- abs(Alength - Blength)
      # if (diff != 0){
      #   len.diff[,paste0(filename, ".",as.character(t-1))] <- 
      # }
      
      filtered.indels[x,bln] <- Alength == Blength
      filtered.indels[x,num] <- diff
      
    }

  }
  
  for (y in 1:nrow(filtered.indels)){
    value <- filtered.indels$total.length[y]
    if (value == 0){
      zeros <- rbind(zeros, data.frame(filtered.indels[y,]))
    }
  }
  
  

  #write.csv(filtered.indels, filename)

}