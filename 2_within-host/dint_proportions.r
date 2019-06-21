
# DINUCLEOTIDE PROPORTIONS 
# -------------------------------------------
nucleotides <- c("A","C","G","T")


total.ins <- read.csv("~/PycharmProjects/hiv-withinhost/12_dinucleotide/total-ins.csv", stringsAsFactors = F)
total.del <- read.csv("~/PycharmProjects/hiv-withinhost/12_dinucleotide/total-del.csv", stringsAsFactors = F)


dinucleotide <- function(seq){
  nucleotides <- c("A","C","G","T")
  ditotal <- nchar(seq)-1
  dinucl <- matrix(nrow=4,ncol=4)
  dinucl[is.na(dinucl)] <- 0
  colnames(dinucl) <- nucleotides
  rownames(dinucl) <- nucleotides
  for (n in 1:ditotal){
    di <- substr(seq,n,n+1)
    pos1 <- substr(di,1,1)
    pos2 <- substr(di,2,2)
    dinucl[pos2,pos1] <- dinucl[pos2,pos1] + 1
  }
  c(dinucl)
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


itot <- colSums(di.i)[1:16]/colSums(di.i)[['sum']]
ivtot <- colSums(di.iv)[1:16]/colSums(di.iv)[['sum']]
dtot <- colSums(di.d)[1:16]/colSums(di.d)[['sum']]
dvtot <- colSums(di.dv)[1:16]/colSums(di.dv)[['sum']]

di.ifinal <- data.frame(seq=itot, vloop=ivtot)
di.dfinal <- data.frame(seq=dtot, vloop=dvtot)

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



cex=2
par(pty="s", mar=c(6,8,4,1),las=0)

lim = c(0.0,0.15)
plot(NA, xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Deletions")
text(x=di.dfinal$vloop, y=di.dfinal$seq, labels=rownames(di.dfinal), cex=1.5)
title(ylab="Proportion In Indels", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
#legend(0.38,0.24,legend=nucleotides, pch=21,cex=1.9, pt.bg=ins.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
#legend(0.10,0.45,legend=c("Insertions", "Deletions"), pch=c(22,23),cex=1.9, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=F)
abline(0,1)