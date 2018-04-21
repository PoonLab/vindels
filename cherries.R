require(ape)
require(stringr)
setwd("~/PycharmProjects/hiv-evolution-master/9_3_indels")
tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/8_4_Filtered_Run2", full.names=TRUE)
vfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/3RegionSequences/VRegions3", full.names=TRUE)

zeros <- data.frame()
len.diff <- list() #data.frame(subtype=character(),stringsAsFactors = F)
nt.prop <- data.frame(stringsAsFactors = F)

for (i in 1:length(tfolder)){
  tre <- read.tree(tfolder[i])
  csv <- read.csv(vfolder[i], header=FALSE)
  
  #for output 
  name <- strsplit(tfolder[i], "/")[[1]][7]
  subtype <- strsplit(name, "\\.")[[1]][1]
  filename <- paste0(subtype, ".csv" )
  
  
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

  
  
  temp <- data.frame(stringsAsFactors = F)
  for (x in 1:nrow(filtered.indels)){
    idxA <- match(filtered.indels$tip1.label[x], csv$accno)
    idxB <- match(filtered.indels$tip2.label[x], csv$accno)
    
    for (t in 2:ncol(csv)){
      Avr <- as.character(csv[idxA,t])
      Bvr <- as.character(csv[idxB,t])
      
      Alength <- nchar(Avr)
      Blength <- nchar(Bvr)
      bln <- Alength == Blength
      
      A.B <- paste0(Avr,Bvr)
      temp[x,paste0("VR",(t-1),".A")] <- str_count(A.B, "A")/nchar(A.B)
      temp[x,paste0("VR",(t-1),".C")] <- str_count(A.B, "C")/nchar(A.B)
      temp[x,paste0("VR",(t-1),".T")] <- str_count(A.B, "T")/nchar(A.B)
      temp[x,paste0("VR",(t-1),".G")] <- str_count(A.B, "G")/nchar(A.B)
      temp[x,paste0("VR",(t-1),".indel")] <- bln
      
      
      name.bln <- paste0("VR",as.character(t-1),".indel")
      name.nt <- paste0("VR",as.character(t-1),".nt")

      diff <- abs(Alength - Blength)
      filtered.indels[x,name.bln] <- bln
      filtered.indels[x,name.nt] <- diff
    }
  
  }
  new.df <- filtered.indels[filtered.indels$total.length != 0,]
  
  write.csv(new.df, filename)
  
  #load the nt proportions data frame from the temp
  # for (g in 1:5){
  #   no <- which(temp[,paste0("VR",g,'.indel')])
  #   yes <- which(!temp[,paste0("VR",g,'.indel')])
  # 
  #   nt.prop[i,paste0('VR',g,'.A.no')] <- mean(temp[,paste0('VR',g,'.A')][no])
  #   nt.prop[i,paste0('VR',g,'.G.no')] <- mean(temp[,paste0('VR',g,'.G')][no])
  #   nt.prop[i,paste0('VR',g,'.T.no')] <- mean(temp[,paste0('VR',g,'.T')][no])
  #   nt.prop[i,paste0('VR',g,'.C.no')] <- mean(temp[,paste0('VR',g,'.C')][no])
  #   if (length(yes) == 0){
  #     nt.prop[i,paste0('VR',g,'.A.yes')] <- 0
  #     nt.prop[i,paste0('VR',g,'.G.yes')] <- 0
  #     nt.prop[i,paste0('VR',g,'.T.yes')] <- 0
  #     nt.prop[i,paste0('VR',g,'.C.yes')] <- 0
  #   }else{
  #     nt.prop[i,paste0('VR',g,'.A.yes')] <- mean(temp[,paste0('VR',g,'.A')][yes])
  #     nt.prop[i,paste0('VR',g,'.G.yes')] <- mean(temp[,paste0('VR',g,'.G')][yes])
  #     nt.prop[i,paste0('VR',g,'.T.yes')] <- mean(temp[,paste0('VR',g,'.T')][yes])
  #     nt.prop[i,paste0('VR',g,'.C.yes')] <- mean(temp[,paste0('VR',g,'.C')][yes])
  #   }
  # }
  # #Used to build the len.diff data frame needed for 3/6+ comparison
  # for (j in 1:5){
  #   len.diff[[paste0(filename,".VR",j,".three")]] <- filtered.indels[,paste0("VR",j,".nt")] <= 3
  #   len.diff[[paste0(filename, ".VR",j,".!three")]] <-  filtered.indels[,paste0("VR",j,".nt")] > 3
  # }
  # if (filename == "B.csv"){
  #   break
  # }
  
}
  
  
  #Used to isolate only the 0 length cherries
  # for (y in 1:nrow(filtered.indels)){
  #   value <- filtered.indels$total.length[y]
  #   if (value == 0){
  #     zeros <- rbind(zeros, data.frame(filtered.indels[y,]))
  #   }
  # }
  #write.csv(filtered.indels, filename)


#Used to load the indel.sizes data frame containing 3/6+ indel frequencies
indel.sizes <- data.frame(stringsAsFactors = FALSE)
count <- 1
for (z in 1:length(len.diff)){
  subtype <- strsplit(names(len.diff)[[z]], "\\.")[[1]][1]
  vregion <- strsplit(names(len.diff)[[z]], "\\.")[[1]][3]
  size <- strsplit(names(len.diff)[[z]], "\\.")[[1]][4]
  
  indel.sizes[z,"subtype"] <- subtype #vregion
  indel.sizes[z,"count"] <- sum(len.diff[[z]])
  indel.sizes[z,"vregion"] <- vregion
  if (size == "three"){
    indel.sizes[z,"three"] <- "<= 3 nt"
  }else{
    indel.sizes[z,"three"] <- "> 3 nt"
  }
  
  #old version
  #indel.sizes[count,paste0(vregion, "_", size)] <- sum(len.diff[[z]])
  # if ((z - 1) %% 10 == 9){
  #   indel.sizes[count,1] <- subtype
  #   count <- count + 1
  # }
}



# for the generation of figures used in poster
df3 <- data.frame(VRegions=rep(indel.sizes$vregion, indel.sizes$count), Indel_Sizes=rep(indel.sizes$three,indel.sizes$count))
table2 <- table(rep(indel.sizes$subtype, indel.sizes$count), rep(indel.sizes$three, indel.sizes$count))
mosaic(~Variable_Region+Indel_Sizes,data=df3, shade=T, gp_labels=(gpar(fontsize=11)))

par(ps = 29, cex.lab = 1.1, cex.axis = 0.5, cex.sub=0.5, las=0, xpd=T, mar=c(5,4, 0.2,2))
par(ps = 29, cex.lab = 1.1, cex.axis = 0.5, cex.sub=0.5, las=0, xpd=T, mar=c(5,4, 2,2))
mosaicplot(~VRegions + Indel_Sizes, data=df3, xlab = "Variable Region", ylab = "Indel Sizes",
           shade=T, main=NULL)
mosaicplot(~Subtype + Indel_Sizes, data=df3, xlab = "Subtype", ylab = "Indel Sizes",
           shade=T, main=NULL)
