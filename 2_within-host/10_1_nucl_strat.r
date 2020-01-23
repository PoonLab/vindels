#
nucleotides <- c("A","C","G","T")
path <- "~/Lio/"
path <- "~/PycharmProjects/hiv-withinhost/"
total.ins <- read.csv(paste0(path,"10_nucleotide/total-ins.csv"), row.names = 1, stringsAsFactors = F)
total.del <- read.csv(paste0(path,"10_nucleotide/total-del.csv"), row.names = 1, stringsAsFactors = F)


# RANDOMIZATION TEST 
# -----------------------------------------
sampleString <- function(len, vloop){
  len <- len-1
  idx <- sample(1:(nchar(vloop)-len),100, replace=TRUE)
  strings <- sapply(idx, function(x){substr(vloop, x, x+len)})
  a <- unname(sapply(strings, function(x){str_count(x, "A")/nchar(x)}))
  c <- unname(sapply(strings, function(x){str_count(x, "C")/nchar(x)}))
  g <- unname(sapply(strings, function(x){str_count(x, "G")/nchar(x)}))
  t <- unname(sapply(strings, function(x){str_count(x, "T")/nchar(x)}))
  list(a,c,g,t)
}
ins3 <- total.ins[nchar(total.ins$Seq) %% 3 ==0,]
insnon3 <- total.ins[nchar(total.ins$Seq) %% 3 !=0,]
del3 <- total.del[nchar(total.del$Seq) %% 3 ==0,]
delnon3 <- total.del[nchar(total.del$Seq) %% 3 !=0,]

ins.3.t <- c(sum(unname(sapply(ins3[,'Seq'], nchar))), sum(unname(sapply(ins3[,"Vseq"], nchar))))
ins.non3.t <- c(sum(unname(sapply(insnon3[,'Seq'], nchar))), sum(unname(sapply(insnon3[,"Vseq"], nchar))))
del.3.t <- c(sum(unname(sapply(del3[,"Seq"], nchar))), sum(unname(sapply(del3[,"Vseq"], nchar))))
del.non3.t <- c(sum(unname(sapply(delnon3[,"Seq"], nchar))), sum(unname(sapply(delnon3[,"Vseq"], nchar))))

iProps <- c()
dProps <- c()
iVProps <- c()
dVProps <- c()
counts <- data.frame()

df1 <- del3
df2 <- delnon3
total1 <- del.3.t
total2 <- del.non3.t

for (nuc in nucleotides){
  icount <- sum(str_count(df1$Seq, nuc))
  dcount <- sum(str_count(df2$Seq, nuc))
  counts <- rbind(counts, data.frame(nucl=nuc, ins=icount, del=dcount))
  
  iProps <- c(iProps, icount / total1[1])
  dProps <- c(dProps, dcount / total2[1])
  
  iVProps <- c(iVProps, sum(str_count(df1$Vseq, nuc)) / total1[2])
  dVProps <- c(dVProps, sum(str_count(df2$Vseq, nuc)) / total2[2])
}
require(reshape)
counts <- melt(counts)

ins.nt <- data.frame(nt=nucleotides,props=iProps,vprops=iVProps)
del.nt <- data.frame(nt=nucleotides,props=dProps,vprops=dVProps)
indel.nt <- rbind(ins.nt, del.nt)
indel.nt$indel <- c(rep(0,4),rep(3,4))
indel.nt$counts <- counts$value



# BOOTSTRAP SUPPORT  
# ------------------------------------
# involves taking random samples WITH replacement 100s or 1000s of times 

bs.props <- list()

for (n in 1:1000){
  
  # create the bootstrap sample 
  sam1 <- sample(nrow(df1), nrow(df1), replace=T)
  sam2 <- sample(nrow(df2), nrow(df2), replace=T)
  
  df1.bs <- df1[sam1,c('Seq','Vseq')]
  df2.bs <- df2[sam2,c('Seq','Vseq')]
  
  total1 <- c(sum(unname(sapply(df1.bs[,"Seq"], nchar))), sum(unname(sapply(df1.bs[,"Vseq"], nchar))))
  total2 <- c(sum(unname(sapply(df2.bs[,"Seq"], nchar))), sum(unname(sapply(df2.bs[,"Vseq"], nchar))))
  
  props <- sapply(nucleotides, function(nuc){
    icount <- sum(str_count(df1.bs$Seq, nuc))
    dcount <- sum(str_count(df2.bs$Seq, nuc))
    
    iProps <-  icount / total1[1]
    dProps <-  dcount / total2[1]
    
    iVProps <-  sum(str_count(df1$Vseq, nuc)) / total1[2]
    dVProps <-  sum(str_count(df2$Vseq, nuc)) / total2[2]
    
    return(c(iProps, dProps, iVProps, dVProps))
  })
  letter <- "d"
  for (n in nucleotides){
    bs.props[[paste0(letter,"3-",n)]] <- c(bs.props[[paste0(letter,"3-",n)]], unname(props[1,n]))
    bs.props[[paste0(letter,"n3-",n)]] <- c(bs.props[[paste0(letter,"n3-",n)]], unname(props[2,n]))
    bs.props[[paste0(letter,"v3-",n)]] <- c(bs.props[[paste0(letter,"v3-",n)]], unname(props[3,n]))
    bs.props[[paste0(letter,"vn3-",n)]] <- c(bs.props[[paste0(letter,"vn3-",n)]], unname(props[4,n]))
  }
}

con.int <- lapply(bs.props, function(x){quantile(x, c(0.025,0.975))})


#RANDOMIZATION TEST 
# ---------------------------------
iSample <- list(c(),c(),c(),c())
dSample <- list(c(),c(),c(),c())
# generates the randomly sampled substrings for each indel
for (row in 1:nrow(df1)){
  itemp <- sampleString(df1[row,"len"], df1[row,"Vseq"])
  for (i in 1:4){
    iSample[[i]] <- c(iSample[[i]], itemp[[i]])
  }
}
for (row in 1:nrow(df2)){
  dtemp <- sampleString(df2[row,"len"], df2[row,"Vseq"])
  for (i in 1:4){
    dSample[[i]] <- c(dSample[[i]], dtemp[[i]])
  }
}

# compares the observed proportion to the overall distribution of each nucleotide 
isign <- c()
dsign <- c()
for (i in 1:4){
  idist <- iSample[[i]]
  ddist <- dSample[[i]]
  
  iQT <- quantile(idist, probs=c(0.025,0.975))
  dQT <- quantile(ddist, probs=c(0.025,0.975))
  
  ins.p <- ins.nt[i,2]
  del.p <- del.nt[i,2]
  
  # highlight significant differences 
  if (ins.p < iQT[[1]]){
    isign <- c(isign, "lower")
  }else if(ins.p > iQT[[2]]){
    isign <- c(isign, "higher")
  }else{
    isign <- c(isign, "")
  }
  
  # highlight significant differences 
  if (del.p < dQT[[1]]){
    dsign <- c(dsign, "lower")
  }else if(del.p > dQT[[2]]){
    dsign <- c(dsign, "higher")
  }else{
    dsign <- c(dsign, "")
  }
  
}

indel.nt$sign <- c(isign,dsign)


# NT ALL INDELS PLOT 
# broken down by variable loop and nucleotide
# -------------------------------------
require(RColorBrewer)
colors <- brewer.pal(4, "Set1")

cex=1
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0.1,0.60)
plot(indel.nt[,c(3,2)], pch=indel.nt[,4]+21, bg=indel.nt[,1],xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.2,cex.main=1.8, ylab='', xlab='',cex=3, main="Deletions - Nt Proportions")
title(ylab="Proportion Inside Indels", line=3,cex.lab=1.3)
title(xlab="Proportion in Variable Loops", line=3,cex.lab=1.3)
legend(0.53,0.24,legend=nucleotides, pch=22,cex=1.3, pt.bg=indel.nt[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.10,0.58,legend=c("3", "Non-3"), pch=c(21,24),cex=1.3, pt.bg="black",x.intersp = 1.0,y.intersp=1.3, pt.cex=3)
par(xpd=F)
abline(0,1)



# NT PROPORTIONS BY VARIABLE LOOP 
# --------------------------------------------
vloops <- c(1, 2,3,4,5)
vloops2 <- c("V1","V2","V4","V5")
ins.props <- data.frame()
del.props <- data.frame()

for (i in c(1,2,4,5)){
  iTemp <- df1[df1$Vloop==i,]
  dTemp <- df2[df2$Vloop==i,]
  
  if (nrow(iTemp) != 0){
    iProps <- c()
    iVProps <- c()
    # a vector of two totals
    # iTotals[1] = total number of nucleotides in insertion sequences
    # iTotals[2] = total nuimebr of nucleotides in the vloops associated with insertion sequences
    iTotals <- c(sum(unname(sapply(iTemp$Seq, nchar))), sum(unname(sapply(iTemp$Vseq, nchar))))
    for (nuc in nucleotides){
      iProps <- c(iProps, sum(str_count(iTemp$Seq, nuc)) / iTotals[1])
      iVProps <- c(iVProps, sum(str_count(iTemp$Vseq, nuc)) / iTotals[2])
    }
    ins.props <- rbind(ins.props, data.frame(nt=nucleotides, iprops=iProps, vprops=iVProps, vloop=rep(vloops[i],4)))
  }
  
  if (nrow(dTemp) != 0){
    dProps <- c()
    dVProps <- c()
    dTotals <- c(sum(unname(sapply(dTemp$Seq, nchar))),sum(unname(sapply(dTemp$Vseq, nchar))))
    
    for (nuc in nucleotides){
      dProps <- c(dProps, sum(str_count(dTemp$Seq, nuc)) / dTotals[1])
      dVProps <- c(dVProps, sum(str_count(dTemp$Vseq, nuc)) / dTotals[2])
    }
    del.props <- rbind(del.props, data.frame(nt=nucleotides, dprops=dProps, vprops=dVProps, vloop=rep(vloops[i],4)))
  }
  

}


# NT PROP INSERTIONS PLOT 
# broken down by variable loop and nucleotide
# -------------------------------------
require(RColorBrewer)
colors <- brewer.pal(4, "Set1")

cex=2
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0.1,0.5)
plot(ins.props[,c(3,2)], pch=ins.props[,4]+20, bg=ins.props[,1],xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.8, main="Deletions")
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Inside Deletions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.43,0.22,legend=nucleotides, pch=21,cex=1.5, pt.bg=ins.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.33,0.22,legend=vloops2, pch=c(21,22,24,25),cex=1.5, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=F)
abline(0,1)


# NT PROP DELETIONS PLOT 
# broken down by variable loop and nucleotide
# -------------------------------------
require(RColorBrewer)
colors <- brewer.pal(4, "Set1")

#A
cex=2
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0.0,0.70)
plot(del.props[,c(3,2)], pch=del.props[,4]+20, bg=del.props[,1],xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Deletions - Non 3")
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Inside Deletions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.44,0.24,legend=nucleotides, pch=21,cex=1.5, pt.bg=del.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.35,0.24,legend=vloops2, pch=c(21,22,24,25),cex=1.5, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=F)
abline(0,1)

