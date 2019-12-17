# DINUCLEOTIDE PROPORTIONS 
# -------------------------------------------
nucleotides <- c("A","C","G","T")


total.ins <- read.csv("~/Lio/10_nucleotide/total-ins.csv", stringsAsFactors = F, row.names=1)
total.del <- read.csv("~/Lio/10_nucleotide/total-del.csv", stringsAsFactors = F, row.names = 1)


dinucleotide <- function(seq){
  nucleotides <- c("A","C","G","T")
  ditotal <- nchar(seq)-1
  dinucl <- matrix(nrow=4,ncol=4)
  dinucl[is.na(dinucl)] <- 0
  colnames(dinucl) <- nucleotides
  rownames(dinucl) <- nucleotides
  if (nchar(seq) > 1){
    for (n in 1:ditotal){
      di <- substr(seq,n,n+1)
      pos1 <- substr(di,1,1)
      pos2 <- substr(di,2,2)
      print(pos1)
      print(pos2)
      dinucl[pos1,pos2] <- dinucl[pos1,pos2] + 1
    }
    
  }
  return(dinucl)
}



di.i <- t(sapply(total.ins$Seq, dinucleotide))
di.iv <- t(sapply(total.ins$Vseq, dinucleotide))
di.d <- t(sapply(total.del$Seq, dinucleotide))
di.dv <- t(sapply(total.del$Vseq, dinucleotide))

rownames(di.i) <- NULL
rownames(di.iv) <- NULL
rownames(di.d) <- NULL
rownames(di.dv) <- NULL

di.i <- as.data.frame(di.i)
di.iv <- as.data.frame(di.iv)
di.d <- as.data.frame(di.d)
di.dv <- as.data.frame(di.dv)

nts <- c()
for (x in 1:4){for (y in 1:4){nts<- c(nts, paste0(nucleotides[x],nucleotides[y]))}}

names(di.i) <- nts
names(di.iv) <- nts
names(di.d) <- nts
names(di.dv) <- nts


di.i$sum <- rowSums(di.i)
di.iv$sum <- rowSums(di.iv)
di.d$sum <- rowSums(di.d)
di.dv$sum <- rowSums(di.dv)

di.ifinal <- data.frame(seq=colSums(di.i)[1:16]/colSums(di.i)[['sum']], vloop=colSums(di.iv)[1:16]/colSums(di.iv)[['sum']])
di.dfinal <- data.frame(seq=colSums(di.d)[1:16]/colSums(di.d)[['sum']], vloop=colSums(di.dv)[1:16]/colSums(di.dv)[['sum']])


# RESIDUAL ANALYSIS






# RANDOMIZATION TEST 
# -----------------------------------------
sampleString <- function(len, vloop){
  nucl <- c("A","C", "G" ,"T")
  len <- len-1
  idx <- sample(1:(nchar(vloop)-len),100, replace=TRUE)
  strings <- sapply(idx, function(x){substr(vloop, x, x+len)})
  result <- list()
  for (a in 1:4){
    for (b in 1:4){
      dinucl <- paste0(nucl[a], nucl[b])
      result[[dinucl]] <- unname(sapply(strings, function(x){str_count(x, dinucl)/nchar(x)}))
    }
  }
  result
}

iSample <- list()
dSample <- list()


# generates the randomly sampled substrings for each indel
for (row in 1:nrow(total.ins)){
  itemp <- sampleString(total.ins[row,"len"], total.ins[row,"Vseq"])
  dtemp <- sampleString(total.del[row,"len"], total.del[row,"Vseq"])
  for (i in 1:16){
    iSample[[nts[i]]] <- c(iSample[[nts[i]]], itemp[[nts[i]]])
    dSample[[nts[i]]] <- c(dSample[[nts[i]]], dtemp[[nts[i]]])
  }
}

# compares the observed proportion to the overall distribution of each nucleotide 
isign <- c()
dsign <- c()
for (i in 1:16){
  idist <- iSample[[nts[i]]]
  ddist <- dSample[[nts[i]]]
  
  iQT <- quantile(idist, probs=c(0.025,0.975))
  dQT <- quantile(ddist, probs=c(0.025,0.975))
  
  ins.p <- di.ifinal[i,1]
  del.p <- di.dfinal[i,1]
  
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
di.ifinal$sign <- isign
di.dfinal$sign <- dsign


# DINUCLEOTIDE FIGURES 
# -----------------------------------------
cex=2
par(pty="s", mar=c(6,8,4,1),las=0)

lim = c(0.0,0.15)
plot(NA, xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Insertions")
text(x=di.ifinal$vloop, y=di.ifinal$seq, labels=rownames(di.ifinal), cex=1.5)
title(ylab="Proportion In Indels", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
#legend(0.38,0.24,legend=nucleotides, pch=21,cex=1.9, pt.bg=ins.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
#legend(0.10,0.45,legend=c("Insertions", "Deletions"), pch=c(22,23),cex=1.9, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=F)
abline(0,1)


ex <- rep(1.5,16)
ex[7] <- 2.5

cex=2
par(pty="s", mar=c(6,8,4,1),las=0)

lim = c(0.0,0.15)
plot(NA, xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Deletions")
text(x=di.dfinal$vloop, y=di.dfinal$seq, labels=rownames(di.dfinal), cex=ex)
title(ylab="Proportion In Indels", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
#legend(0.38,0.24,legend=nucleotides, pch=21,cex=1.9, pt.bg=ins.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
#legend(0.10,0.45,legend=c("Insertions", "Deletions"), pch=c(22,23),cex=1.9, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=F)
abline(0,1)