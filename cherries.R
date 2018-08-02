require(ape)
require(stringr)
require(ggplot2)

tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/8_Dated_Trees", full.names=TRUE)
vfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/3RegionSequences/VRegions_edit", full.names=TRUE)

len.diff <- list() #data.frame(subtype=character(),stringsAsFactors = F)


# This code reads phylogenetic trees and variable loop sequences in csv format. 
branch.lengths <- data.frame()

for (i in 1:length(tfolder)){
  tre <- read.tree(tfolder[i])
  csv <- read.csv(vfolder[i], header=FALSE, stringsAsFactors = F)
  
  #for output 
  name <- strsplit(tfolder[i], "/")[[1]][7]
  subtype <- strsplit(name, "\\.")[[1]][1]
  filename <- paste0(subtype,"+.csv" )
  
  
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


  df$tip1.label <- tre$tip.label[df$tip1]
  df$tip2.label <- tre$tip.label[df$tip2]
  df$tip1.len <- tre$edge.length[df$edge1]
  df$tip2.len <- tre$edge.length[df$edge2]
  

  indels <- df[,c(6:9)]
  indels$total.length <- indels$tip1.len + indels$tip2.len
  
  
  filtered.indels <- indels[indels$total.length != 0,]
  
  filtered.indels2 <- data.frame(stringsAsFactors = F)
  temp3 <- data.frame(stringsAsFactors = F)
  lens <- c(0,78,120,108,102,33)
  
  count = 0
  for (x in 1:nrow(filtered.indels)){
    idxA <- match(filtered.indels$tip1.label[x], csv$accno)
    idxB <- match(filtered.indels$tip2.label[x], csv$accno)
    
    
    temp2 <- data.frame(stringsAsFactors = F)
    
    #BY VARIABLE LOOPS / COLUMNS 
    for (t in 2:ncol(csv)){

      Avr <- as.character(csv[idxA,t])
      Bvr <- as.character(csv[idxB,t])
      
      Alength <- nchar(Avr)
      Blength <- nchar(Bvr)
      bln <- Alength == Blength
      
      A.B <- paste0(Avr,Bvr)
      
      name.bln <- paste0("VR",as.character(t-1),".indel")
      name.nt <- paste0("VR",as.character(t-1),".nt")
      
      
      
      #FILTER ----------------------------
      if (Alength < lens[t]/2 || Blength < lens[t]/2 || 
          (str_count(Avr, "\\?")/Alength) > 0.15 || 
          (str_count(Bvr, "\\?")/Blength) > 0.15 ||
          (t !=6 && (substr(paste(trans(as.DNAbin.DNAString(Avr)),collapse=""),1,1) != "C"|| 
                     substr(paste(trans(as.DNAbin.DNAString(Bvr)),collapse=""),1,1) != "C")))
      {
        count = count + 1
        print(Avr)
        print(Bvr)
        print("")
        filtered.indels[x,name.bln] <- NA
        filtered.indels[x,name.nt] <- NA
      }else{
        diff <- abs(Alength - Blength)
        filtered.indels[x,name.bln] <- bln
        filtered.indels[x,name.nt] <- diff
      }
      
      #for generating 10_Cherries, only indel-containing sequences
      if (!is.na(bln) && !bln){
        temp2 <- rbind(temp2, data.frame(accno1=as.character(csv[idxA,1]), seq1=Avr, accno2=as.character(csv[idxB,1]), seq2=Bvr, Vr=(t-1)))
      }
      
    } #COLUMNS END
    temp3 <- rbind(temp3, temp2)
    filtered.indels2 <- rbind(filtered.indels2, filtered.indels[x,])
    
    
  } #ROWS END 
  print(count)
  
  setwd("~/PycharmProjects/hiv-evolution-master/10_Cherries2/")
  #temp3 = output for cherry sequences in csv format (accno1,seq1,accno2,seq2,vregion) FOLDER: 10_Cherries
  write.csv(temp3,filename)
  
  
  setwd("~/PycharmProjects/hiv-evolution-master/9_2_indels/")
  #filtered indels = true/false outcomes along with the indel sizes  FOLDER: 9_indels
  write.csv(filtered.indels2, filename)
  
  branch.lengths <- rbind(branch.lengths, data.frame(subtype=rep(subtype,nrow(filtered.indels)), length=filtered.indels$total.length))
  for (j in 1:5){
    len.diff[[paste0(filename,".VR",j,".three")]] <- filtered.indels[which(!is.na(filtered.indels[paste0("VR",j,".indel")])),paste0("VR",j,".nt")] <= 3
    len.diff[[paste0(filename, ".VR",j,".!three")]] <-  filtered.indels[which(!is.na(filtered.indels[paste0("VR",j,".indel")])),paste0("VR",j,".nt")] > 3
  }

}
  
  
#Used to load the indel.sizes data frame containing 3/6+ indel frequencies
indel.sizes <- data.frame(stringsAsFactors = FALSE)

for (z in 1:length(len.diff)){
  subtype <- strsplit(names(len.diff)[[z]], "\\+.")[[1]][1]
  vregion <- strsplit(strsplit(names(len.diff)[[z]], "\\.")[[1]][3],"VR")[[1]][2]
  size <- strsplit(names(len.diff)[[z]], "\\.")[[1]][4]
  
  indel.sizes[z,"subtype"] <- subtype #vregion
  indel.sizes[z,"count"] <- sum(len.diff[[z]])
  indel.sizes[z,"vregion"] <- vregion
  if (size == "three"){
    indel.sizes[z,"three"] <- "<= 3 nt"
  }else{
    indel.sizes[z,"three"] <- "> 3 nt"
  }
  
}


# MOSAIC PLOTS - Figure 3
df3 <- data.frame(VRegions=rep(indel.sizes$vregion, indel.sizes$count), Indel_Sizes=rep(indel.sizes$three,indel.sizes$count))
df4 <- data.frame(Subtype=rep(indel.sizes$subtype, indel.sizes$count), Indel_Sizes=rep(indel.sizes$three,indel.sizes$count))
table2 <- table(rep(indel.sizes$subtype, indel.sizes$count), rep(indel.sizes$three, indel.sizes$count))

require(vcd)
par(ps = 50, cex.lab = 0.7, cex.axis = 0.5, cex.sub=0.5, las=0, xpd=T, mar=c(5,4, 2,2), mfrow=c(2,2))
mosaic(~VRegions + Indel_Sizes, data=df3,
       shade=T, main=NULL, direction="v",spacing=spacing_equal(sp = unit(0, "lines")),
       residuals_type=NULL,legend=T)



          
par(ps = 27, cex.lab = 0.7, cex.axis = 0.5, cex.sub=0.1, las=0, xpd=T, mar=c(5,4, 2,2), xaxt='n')
mosaic(~Subtype + Indel_Sizes, data=df4, xlab = "Subtype", ylab = "Indel Sizes",
       shade=T, main=NULL, direction="v",spacing=spacing_equal(sp = unit(0, "lines")),
       residuals_type=NULL,legend=T, set_varnames = c(Sex="Gender", Survived="survival"))


#PHYLOGENETIC PLOT - Figure 1
tre <- read.tree(tfolder[7])



#-------------------------
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
#Used to build the len.diff data frame needed for 3/6+ comparison

