# DINUCLEOTIDE PROPORTIONS 
# -------------------------------------------
nucleotides <- c("A","C","G","T")


total.ins <- read.csv("~/PycharmProjects/hiv-withinhost/10_nucleotide/all/dinucl-ins-current.csv", stringsAsFactors = F, row.names=1)
total.del <- read.csv("~/PycharmProjects/hiv-withinhost/10_nucleotide/all/dinucl-del-current.csv", stringsAsFactors = F, row.names = 1)

total.ins$anc <- gsub("-","",total.ins$anc)
total.del$anc <- gsub("-","",total.del$anc)

total.ins = total.ins[!grepl("[^ACGT]", total.ins$indel),]
total.del = total.del[!grepl("[^ACGT]", total.del$indel),]

dinucleotide <- function(seq){
  #initialize the matrix 
  nucleotides <- c("A","C","G","T")
  ditotal <- nchar(seq)-1
  dinucl <- matrix(nrow=4,ncol=4)
  dinucl[is.na(dinucl)] <- 0
  colnames(dinucl) <- nucleotides
  rownames(dinucl) <- nucleotides
  #print(seq)
  nts <- c()
  for (x in 1:4){for (y in 1:4){nts<- c(nts, paste0(nucleotides[x],nucleotides[y]))}}
  
  if (nchar(seq) > 1){
    for (n in 1:ditotal){
      di <- substr(seq,n,n+1)
      pos1 <- substr(di,1,1)
      pos2 <- substr(di,2,2)
      #print(pos1)
      #print(pos2)
      dinucl[pos2,pos1] <- dinucl[pos2,pos1] + 1
    }
    
  }
  dinucl <- c(dinucl)
  names(dinucl) <- nts
  dinucl
}



di.i <- lapply(total.ins$indel, dinucleotide)
di.iv <- lapply(total.ins$anc, dinucleotide)
di.d <- lapply(total.del$indel, dinucleotide)
di.dv <- lapply(total.del$anc, dinucleotide)

di.i <- Reduce("+", di.i)
di.iv <- Reduce("+", di.iv)
di.d <- Reduce("+", di.d)
di.dv <- Reduce("+", di.dv)

di.ins <- data.frame(indel=di.i, vloop=di.iv)
di.del <- data.frame(indel=di.d, vloop=di.dv)

di.ins$iprop <- sapply(di.ins[,1], function(x){x/colSums(di.ins)[[1]]})
di.ins$vprop <- sapply(di.ins[,2], function(x){x/colSums(di.ins)[[2]]})
di.del$dprop <- sapply(di.del[,1], function(x){x/colSums(di.del)[[1]]})
di.del$vprop <- sapply(di.del[,2], function(x){x/colSums(di.del)[[2]]})

#di.ifinal <- data.frame(seq=colSums(di.ins)[1:16]/colSums(di.i)[['sum']], vloop=colSums(di.iv)[1:16]/colSums(di.iv)[['sum']])
#di.dfinal <- data.frame(seq=colSums(di.d)[1:16]/colSums(di.d)[['sum']], vloop=colSums(di.dv)[1:16]/colSums(di.dv)[['sum']])


# RESIDUAL ANALYSIS




nt <- c("A","C","G","T")
nts <- c()
for(i in 1:4){for(j in 1:4){nts <- c(nts, paste0(nt[i],nt[j]))}}

# RANDOMIZATION TEST 
# -----------------------------------------
sampleString <- function(len, vloop){
  nucl <- c("A","C", "G" ,"T")
  len <- len-1
  idx <- sample(1:(nchar(vloop)-len),50, replace=TRUE)
  strings <- sapply(idx, function(x){substr(vloop, x, x+len)})
  result <- list()
  sapply(nts, function(x){
    result[[x]] <<- unname(sapply(strings, function(y){str_count(y, x)/nchar(y)}))
  })
  result
}

iSample <- list()
dSample <- list()
for (i in 1:16) {
  iSample[[nts[i]]] <- c()
  dSample[[nts[i]]] <- c()
}
# generates the randomly sampled substrings for each indel


size <- 50
null <- sapply(1:nrow(total.ins), function(i){
  temp <- sampleString(nchar(total.ins[i,'indel']), total.ins[i,"anc"])
  for (n in 1:16){
    iSample[[nts[n]]][((i-1)*size+1):(i*size)] <<- temp[[nts[n]]]
  }
})
size <- 50
null <- sapply(1:nrow(total.del), function(i){
  temp <- sampleString(nchar(total.del[i,'indel']), total.del[i,"anc"])
  for (n in 1:16){
    dSample[[nts[n]]][((i-1)*size+1):(i*size)] <<- temp[[nts[n]]]
  }
})
#   itemp <- sampleString(nchar(total.ins[row,'indel']), total.ins[row,"anc"])
#   dtemp <- sampleString(nchar(total.del[row,'indel']), total.del[row,"anc"])
#   for (i in 1:16){
#     iSample[[nts[i]]] <- c(iSample[[nts[i]]], itemp[[nts[i]]])
#     dSample[[nts[i]]] <- c(dSample[[nts[i]]], dtemp[[nts[i]]])
#   }
# }

# compares the observed proportion to the overall distribution of each nucleotide 
isign <- c()
dsign <- c()
for (i in 1:16){
  idist <- iSample[[nts[i]]]
  ddist <- dSample[[nts[i]]]
  
  iQT <- quantile(idist, probs=c(0.025,0.975))
  dQT <- quantile(ddist, probs=c(0.025,0.975))
  
  ins.p <- di.ins[i,3]
  del.p <- di.del[i,3]
  
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
di.ins$sign <- isign
di.del$sign <- dsign


# DINUCLEOTIDE FIGURES 
# -----------------------------------------
cex=2
par(pty="s", mar=c(6,8,2,1),las=0)
par(xpd=F)
lim = c(0.0,0.15)
plot(NA, xlim=lim,ylim=lim,
     cex.lab=1.75, 
     cex.axis=1.4,
     ylab='', 
     xlab='',
     cex=3.5,
     las=1)
text(x=di.ins[,"vprop"], y=di.ins[,"iprop"], labels=rownames(di.ins), cex=1.5)
#text(x=0.02,y=0.15, cex=1.5, labels="Insertions")
title(ylab="Proportion In Insertions", line=4.2,cex.lab=1.9)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.9)
abline(0,1)
par(xpd=NA)
text(-0.040, 0.15, labels="a)", cex=2)
#legend(0.38,0.24,legend=nucleotides, pch=21,cex=1.9, pt.bg=ins.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
#legend(0.10,0.45,legend=c("Insertions", "Deletions"), pch=c(22,23),cex=1.9, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)

ex <- c(6,10,11,14 )
pnts <- di.del

cex=2
par(pty="s", mar=c(6,8,2,1),las=0)
par(xpd=F)
lim = c(0.0,0.15)
plot(NA, xlim=lim,ylim=lim,
     cex.lab=1.75, 
     cex.axis=1.4,
     ylab='', 
     xlab='',
     cex=3.5,
     las=1)
text(x=di.del[,"vprop"], y=di.del[,"dprop"], labels=rownames(di.del), cex=1.5)
#text(x=0.02,y=0.15, cex=1.5, labels="Insertions")
title(ylab="Proportion In Deletions", line=4.2,cex.lab=1.9)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.9)
# points(x=pnts[,"vprop"], y=pnts[,"dprop"], pch=21, cex=1.3, bg="black")
# xcr <- c(0.017, 0.0307, 0.02,0.0293)
# ycr <- c(0.0293, 0.043, 0.039, 0.017)
# #arrows(pnts[,"vprop"], pnts[,"dprop"], xcr, ycr, length=0)
# text(xcr+c(-0.005,0,-0.003,0),ycr+c(0,0.003,0.003,-0.004), labels=rownames(pnts), cex=1.2)
abline(0,1)
par(xpd=NA)
text(-0.040, 0.15, labels="b)", cex=2)
