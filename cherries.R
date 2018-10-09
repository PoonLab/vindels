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
  test <- strsplit(subtype,"_")[[1]]
  if (length(test) == 2){
    filename <- paste0(subtype,"+.csv" )
    subtype <- strsplit(subtype,"_")[[1]][2]
  }else if(subtype == "F1"){
    filename <- paste0(subtype,"+.csv" )
    subtype <- "F" 
  }else{
    filename <- paste0(subtype,"+.csv" )
  }
  
  
  
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
  indel2 <- data.frame(stringsAsFactors = F)
  nonindel2 <- data.frame(stringsAsFactors = F)
  lens <- c(0,78,120,108,102,33)
  
  count = 0
  for (x in 1:nrow(filtered.indels)){
    idxA <- match(filtered.indels$tip1.label[x], csv$accno)
    idxB <- match(filtered.indels$tip2.label[x], csv$accno)
    
    
    indel <- data.frame(stringsAsFactors = F)
    nonindel <- data.frame(stringsAsFactors = F)
    
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
        indel <- rbind(indel, data.frame(accno1=as.character(csv[idxA,1]), seq1=Avr, accno2=as.character(csv[idxB,1]), seq2=Bvr, Vr=(t-1)))
      #for generating 9_0_nonindel, containing only nonindel sequences
      }else if (!is.na(bln) && isTRUE(bln)){
        nonindel <- rbind(nonindel, data.frame(accno1=as.character(csv[idxA,1]), seq1=Avr, accno2=as.character(csv[idxB,1]), seq2=Bvr, Vr=(t-1)))
      }
      
    } #COLUMNS END
    indel2 <- rbind(indel2, indel)
    nonindel2 <- rbind(nonindel2, nonindel)
    filtered.indels2 <- rbind(filtered.indels2, filtered.indels[x,])
    
    
  } #ROWS END 
  print(count)
  
  #indel2 = output for cherry sequences in csv format (accno1,seq1,accno2,seq2,vregion) FOLDER: 10_Cherries
  #only cherries with AT LEAST ONE INDEL are directed here
  setwd("~/PycharmProjects/hiv-evolution-master/10_Cherries/")
  write.csv(indel2,filename)
  
  
  setwd("~/PycharmProjects/hiv-evolution-master/9_0_nonindel/")
  write.csv(nonindel2,paste0(filename,"+"))
  
  #filtered indels = true/false outcomes + indel sizes --- used for MLE analysis  FOLDER: 9_indels
  setwd("~/PycharmProjects/hiv-evolution-master/9_2_indels/")
  write.csv(filtered.indels, filename)
  
  branch.lengths <- rbind(branch.lengths, data.frame(subtype=rep(subtype,nrow(filtered.indels)), length=filtered.indels$total.length))
  for (j in 1:5){
    name <- paste0("VR",j,".indel")
    values <- filtered.indels[which(!is.na(filtered.indels[name]) & !filtered.indels[name]),paste0("VR",j,".nt")]
    len.diff[[paste0(filename,".VR",j,".three")]] <- values == 3
    len.diff[[paste0(filename, ".VR",j,".six")]] <-  values == 6
    len.diff[[paste0(filename, ".VR",j,".nine")]] <-  values >= 9
  }

}
  
  
#Used to load the indel.sizes data frame containing 3/6+ indel frequencies
indel.sizes <- data.frame(stringsAsFactors = FALSE)

for (z in 1:length(len.diff)){
  subtype <- strsplit(names(len.diff)[[z]], "\\+.")[[1]][1]
  vregion <- strsplit(strsplit(names(len.diff)[[z]], "\\.")[[1]][3],"VR")[[1]][2]
  size <- strsplit(names(len.diff)[[z]], "\\.")[[1]][4]
  
  if (subtype == "01_AE"){
    subtype <- "AE"
  }else if(subtype == "02_AG"){
    subtype <- "AG"
  }else if(subtype == "F1"){
    subtype <- "F"
  }
  
  indel.sizes[z,"subtype"] <- subtype #vregion
  indel.sizes[z,"count"] <- sum(len.diff[[z]])
  indel.sizes[z,"vregion"] <- paste0("V",vregion)
  if (size == "three"){
    indel.sizes[z,"size"] <- "3"
  }else if (size == "six"){
    indel.sizes[z,"size"] <- "6"
  }else{
    indel.sizes[z,"size"] <- "9+"
  }
  
}


# MOSAIC PLOTS - Figure 3
df3 <- data.frame(variable.loop=rep(indel.sizes$vregion, indel.sizes$count), indel.size=rep(indel.sizes$size,indel.sizes$count),stringsAsFactors = F)
df4 <- data.frame(subtype=rep(indel.sizes$subtype, indel.sizes$count), indel.size=rep(indel.sizes$size,indel.sizes$count),stringsAsFactors = F)

df3$indel.size <- factor(df3$indel.size,levels=c("9+","6","3"))
df3 <- df3[order(df3$indel.size),]

df4$indel.size <- factor(df4$indel.size,levels=c("9+","6","3"))
df4$subtype <- factor(df4$subtype, levels=c("AE", "AG", "A1","B","C","D","F"))
df4 <- df4[order(df4$indel.size),]

require(vcd)
library(gridGraphics)
library(gridExtra)
#par(ps = 50, cex.lab = 0.7, cex.axis = 0.5, cex.sub=0.5, las=0, xpd=T, mar=c(5,4, 2,2), mfrow=c(2,2))


mosaic(~variable.loop + indel.size, data=df3,
       shade=T, main=NULL, direction="v",
       spacing=spacing_equal(sp = unit(0.7, "lines")),
       residuals_type="Pearson",
       margins=c(2,2,4,2),
       labeling_args = list(tl_labels = c(F,T), 
                            tl_varnames=c(F,T),
                            gp_labels=gpar(fontsize=19),
                            gp_varnames=gpar(fontsize=24),
                            set_varnames = c(variable.loop="Variable Loop", 
                                             indel.size="Indel Length (nt)"),
                            offset_labels=c(0,0,0,0),rot_labels=c(0,0,0,0), just_labels=c("center","center","center","center")),
       legend=legend_resbased(fontsize = 16, fontfamily = "",
                       x = unit(0.2, "lines"), y = unit(3,"lines"),
                       height = unit(0.8, "npc"),
                       width = unit(1, "lines"), range=c(-10,10),pvalue=F),
       set_labels=list(Variable.Loop=c("V1","V2","V3","V4","V5")))
a <- grid.grab()

#par(ps = 27, cex.lab = 0.7, cex.axis = 0.5, cex.sub=0.1, las=0, xpd=T, mar=c(5,4, 2,2), xaxt='n')
mosaic(~subtype + indel.size, data=df4,
       shade=T, main=NULL, direction="v",
       spacing=spacing_equal(sp = unit(0.7, "lines")),
       residuals_type="Pearson",legend=F, 
       margins=c(1,4,4,4),
       labeling_args = list(tl_labels = c(F,T), 
                            tl_varnames=c(F,T), 
                            gp_labels=gpar(fontsize=19),
                            gp_varnames=gpar(fontsize=24),
                            set_varnames = c(subtype="Subtype", 
                                             indel.size="Indel Length (nt)"),
                            offset_labels=c(0,0,0,0),rot_labels=c(0,0,35,0), 
                            just_labels=c("center","center","center","center")))
b <- grid.grab()
grid.newpage()
grid.arrange(b,a,ncol=2)

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

