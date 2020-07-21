# flanking insertion sequences 
# trying to develop a model for interstrand jumping of polymerase 
require(ape)
require(stringr)
source("~/vindels/2_within-host/utils.r")

subs <- function(seq1, seq2){
  # if (nchar(seq1) != nchar(seq2)){
  #   return(NA)
  # }
  c1 <- str_split(seq1, "")[[1]]
  c2 <- str_split(seq2, "")[[1]]

  paste(which(unname(mapply(function(x,y){x!=y}, x=c1, y=c2))),collapse=",")
}

fixHeader <- function(header, pat){
  rep <- strsplit(pat, "-")[[1]][2]
  paste0(header,"_",rep)
}

slips <- function(vseq, pos, len){
  pos <- as.numeric(pos)+1
  len <- as.numeric(len)
  
  chars <- strsplit(vseq,"")[[1]]
  
}


# 14_flanking 
flankCheck <- function(indel,pos,vseq,wobble=1/12, offset=0){
  len <- nchar(indel)
  pos <- as.numeric(pos)
  
  beforeBool <- F
  afterBool <- F
  
  beforeIdx <- NaN
  afterIdx <- NaN
  
  beforeDiff <- NaN
  afterDiff <- NaN
  
  beforeSeq <- ""
  afterSeq <- ""
  
  # BEFORE (5')
  # find the best matching position from 0 to offset 
  lowest <- 1000
  bestIdx <- -1
  
  for (idx in 0:offset){
    # subtract length to get the start of the sequence 
    # needs to be enough nucleotides to check
    if ((pos - len - idx) >= 0){
      before <- substr(vseq, pos-len-idx, pos-idx-1)
      #print(before)
      diffs <- checkDiff(indel, before)
      if (length(diffs) < lowest){
        lowest <- length(diffs)
        bestIdx <- idx
      }
    }
  }
  
  if (bestIdx != -1 & lowest <= round(wobble*len)){
    beforeBool <- T
    beforeIdx <- bestIdx
    beforeDiff <- lowest
    beforeSeq <- substr(vseq,pos-len-bestIdx, pos-bestIdx-1)
  }
  
  # AFTER 
  # find the best matching position from 0 to offset 
  lowest <- 1000
  bestIdx <- -1
  for (idx in 0:offset){
    if ((pos + len + idx) <= nchar(vseq)){
      # then the PRECEDING position can be checked
      after <- substr(vseq, pos+idx, pos+len+idx-1)
      #print(after)
      diffs <- checkDiff(indel, after)
      if (length(diffs) < lowest){
        lowest <- length(diffs)
        bestIdx <- idx
      }
    }
  }
  
  if (bestIdx != -1 & lowest <= round(wobble*len)){
    afterBool <- T
    afterIdx <- bestIdx
    afterDiff <- lowest
    afterSeq <- substr(vseq,pos+bestIdx, pos+len+bestIdx-1)
  }
  
  c(indel, vseq, as.logical(beforeBool),  as.numeric(beforeIdx),  as.numeric(beforeDiff), beforeSeq, as.logical(afterBool), as.numeric(afterIdx), as.numeric(afterDiff),afterSeq)
}

# 14_flanking
flankProps <- function(indel, pos, vseq){
  len <- nchar(indel)
  pos <- as.numeric(pos)
  
  if ((pos - len - idx) >= 0){
    before <- substr(vseq, pos-len-idx+1, pos-idx)
    #print(before)
    diffs <- checkDiff(indel, before)
    if (length(diffs) < lowest){
      lowest <- length(diffs)
      bestIdx <- idx
    }
  }
  
}

ins <- read.csv("~/PycharmProjects/hiv-withinhost/10_nucleotide/all/flanking/ins-sep.csv", stringsAsFactors = F, row.names = 1)
del <- read.csv("~/PycharmProjects/hiv-withinhost/10_nucleotide/all/flanking/del-sep.csv", stringsAsFactors = F, row.names = 1)
#del <- del[-c(which(nchar(del$indel) > 50)),]
#ins <- ins[-c(which(nchar(ins$indel) > 50)),]
lens <- nchar(ins$indel)

del <- del[-which(nchar(del$indel) > nchar(del$anc)),]
del <- del[-which(del$pos > nchar(del$anc)),]

# apply an adjust to the deletion locations to make them the same as insertions 
#del$pos <- as.numeric(del$pos) + nchar(del$indel)

# Insertions : fill in gaps found in the tip sequences 
ins$anc <- gsub("-","",ins$anc)
ins$pos <- ins$pos - nchar(ins$indel) + 1
# Deletions : fill in gaps found in the ancestral sequences 
# this is to include any insertions in the ancestor 

del$anc <- mapply(restoreInsAnc2, del$anc, del$indel, del$pos)
del$anc <- gsub("-","",del$anc)
del$pos <- del$pos - nchar(del$indel) + 1


df <- del

flanking <- unname(mapply(flankCheck, 
                          indel=df$indel,
                          pos=df$pos, 
                          vseq=df$anc, 
                          wobble=0, 
                          offset=10000))

# modify flanking data.frame 
flanking <- as.data.frame(t(flanking), stringsAsFactors = F)
flanking <- cbind(df[,c(2,7)], len=nchar(df$indel), flanking)
colnames(flanking) <- c("vloop", "pos", "len","indel","anc", "before.bool", "before.offset", "before.diff", "before.seq","after.bool", "after.offset", "after.diff",  "after.seq")
flanking[,"before.bool"] <- as.logical(flanking[,"before.bool"] )
flanking[,"after.bool"] <- as.logical(flanking[,"after.bool"] )
flanking[,"before.offset"] <- as.numeric(flanking[,"before.offset"] )
flanking[,"after.offset"] <- as.numeric(flanking[,"after.offset"] )
flanking[,"before.diff"] <- as.numeric(flanking[,"before.diff"] )
flanking[,"after.diff"] <- as.numeric(flanking[,"after.diff"] )

idx <- which(flanking$after.bool | flanking$before.bool)

dist2 <- sapply(1:nrow(flanking), function(i){
  b1 <- flanking$after.bool[i]
  b2 <- flanking$before.bool[i]
  
  if (b1 & b2){
    os1 <- flanking$after.offset[i]
    os2 <- flanking$before.offset[i]
    if (os1 < os2){
      return(os1)
    }else{
      return(os2)
    }
  }else{
    if(b1){
      return(flanking$after.offset[i])
    }else if(b2){
      return(flanking$before.offset[i])
    }else{
      return(100)
    }
  }
})
dist<-dist[dist<=100]
dist2 <- dist2[dist2<=100]
par(mar=c(5,5,5,2))
caxis=1.3
clab=1.4
cmain=1.5

par(mfrow=c(1,2), mar=c(6,5,4,0),las=1,xpd=F)
# distribution of how far you need to travel to find an EXACT MATCH (wobble = 0, offset=10000)
hist(dist, breaks=seq(-0.5,100.5), xlim=c(0,29),
     ylim=c(0,0.25),freq=F, col='red',cex.lab=clab, 
     ylab="",main="Insertion - Exact Match Distance",cex.axis=caxis, 
     cex.main=cmain, xlab="Distance from Insertion Site (nt)")
title(ylab="Density", line=3.5, cex.lab=clab)
par(xpd=NA)
text(-5,0.285, "a)", cex=1.7)
par(mar=c(6,5,4,2),xpd=F)
hist(dist2, breaks=seq(-0.5,100.5), xlim=c(0,29),ylim=c(0,0.25),freq=F, 
     col='red',cex.lab=clab, ylab="",main="Deletion - Exact Match Distance",
     cex.axis=caxis, cex.main=cmain, 
     xlab="Distance from Deletion Site (nt)")
title(ylab="Density", line=3.5, cex.lab=clab)
par(xpd=NA)
text(-5,0.285, "b)", cex=1.7)
par(xpd=F)

sum(flanking$before.bool | flanking$after.bool) / nrow(flanking)

sub <- flanking[flanking$len > 5,]

sum(sub$before.bool | sub$after.bool) / nrow(sub)


# used to remove problematic cases 
all <- all[-as.numeric(rownames(all[all$Vlength==0,])),]
flanking <- flanking[-171,]
all <- all[-13514, ]

# proportion of insertions ACROSS ALL VARIABLE LOOPS that contain a match with 1/6 wobble directly 
#   adjacent (EITHER 5' or 3', no offset) 
sum(flanking$before.bool | flanking$after.bool) / nrow(flanking)


# retrieves insertions that have at least one instance of flanking sequence 
tabs <- table(flanking[flanking$before.bool | flanking$after.bool, "header"])
all$count.flanking <- 0
all[all$Header %in% names(tabs),"count.flanking"] <- tabs

# creation of the new.count column 
all$new.count <- 0
all[all$count.flanking!=0, "new.count"] <- nchar(all[all$count.flanking!=0,"indel"])

# creation of a substitution list column
#all$subs <- unname(mapply(subs, all$anc, all$Anc))

write.csv(flanking, "~/PycharmProjects/hiv-withinhost/14_flanking/flanking.csv")
write.csv(all, "~/PycharmProjects/hiv-withinhost/14_flanking/flanking-all.csv")



# split 'flanking' into 5 different variable loops 
sub.v <- split(flanking, flanking$vloop)

# calculate the probability of slippage occurring in each of the V-loops 
require(R.utils)
p.slip <- c()
for (i in c(1,2,4,5)){
  p.slip[i] <-  nrow(sub.v[[i]][sub.v[[i]]$before.bool | sub.v[[i]]$after.bool,]) / nrow(sub.v[[i]])
}
p.slip[3] <- 0

# probability of a slippage event occurring per nucleotide per v-loop (flanking contains a single insertion event per row)
nt.slip <- c()
for (i in 1:4){
  nt.slip[i] <-nrow(sub.v[[i]][sub.v[[i]]$before.bool | sub.v[[i]]$after.bool,]) / sum(nchar(all[all$Vloop==i & all$indel == "","Vseq"]))
}

# creates a transitional probability matrix for EACH V-LOOP (1,2,4,5 ; not V3) 
all$Vseq <- mapply(addX, all$Vseq, all$pos)
all.v <- split(all, all$Vloop)

p.trans <- list()
nt <- c("A", "C", "G", "T", "X")
for (i in c(1,2,4,5)){
  p.trans[[i]] <- transitionCounts("")
  total <- c(0,0,0,0,0)
  for (j in 1:nrow(all.v[[i]])){
    vseq <- all.v[[i]][j,"Vseq"]
    p.trans[[i]] <- p.trans[[i]] + transitionCounts(vseq)
    for (a in 1:5){
      total[a] <- total[a] + str_count(substr(vseq,1,nchar(vseq)-1), nt[a])
    }
  }
  p.trans[[i]] <- p.trans[[i]] / total
}





par(mfrow=c(3,2))
# 3 x 2 HISTOGRAM PLOT 
# -------------------------------------
# set flanking to infinite offset, 0 wobble
# determine how far you need to go to find a match
flanking[is.nan(flanking$after.offset),"after.offset"] <- 100
flanking[is.nan(flanking$before.offset),"before.offset"] <- 100

a.dist <- flanking[ , "after.offset"]
b.dist <- flanking[, "before.offset"]
a.dist1 <- flanking[!is.na(flanking$after.offset), "after.offset"]
b.dist1 <- flanking[!is.na(flanking$before.offset), "before.offset"]
a.dist2 <- flanking[!is.na(flanking$after.offset), "after.offset"]
b.dist2 <- flanking[!is.na(flanking$before.offset), "before.offset"]

par(mar=c(5,5,5,2))
caxis=1.3
clab=1.4
cmain=1.5

par(mfrow=c(1,2), mar=c(6,6,4,2),las=1)
# distribution of how far you need to travel to find an EXACT MATCH (wobble = 0, offset=10000)
hist(b.dist, breaks=100, xlim=c(0,50),freq=F, col='red',cex.lab=clab, main="Distances to next exact match - 5'",cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")
hist(a.dist, breaks=100,  xlim=c(0,50),freq=F,col='red',cex.lab=clab, main="Distances to next exact match - 3'", cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")


# distribution of how far you need to travel to find a MATCH WITHIN 1 NT (wobble = 1, offset=10000)
hist(b.dist1, breaks=seq(-0.5,max(b.dist1)+0.5), col='red',cex.lab=clab, main="Distances to next match \n(1/12 mismatch) - 5'",cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")
hist(a.dist1, breaks=seq(-0.5,max(a.dist1)+0.5), col='red',cex.lab=clab, main="Distances to next match \n(1/12 mismatch) - 3'", cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")


# distribution of how far you need to travel to find a MATCH WITHIN 2 NT (wobble = 2, offset=10000)
hist(b.dist2, breaks=seq(-0.5,max(b.dist2)+0.5), col='red',cex.lab=clab, main="Distances to next match \n(1/9 mismatch) - 5'",cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")
hist(a.dist2, breaks=seq(-0.5,max(a.dist2)+0.5), col='red',cex.lab=clab, main="Distances to next match \n(1/9 mismatch) - 3'", cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")


lens <- nchar(flanking$indel)

par(mar=c(6,6,6,2))
caxis=1.3
clab=1.5
cmain=1.8
hist(lens, breaks=seq(-0.5,max(lens)+0.5), col='red',cex.lab=clab, main="Insertion Lengths", cex.axis=caxis, cex.main=cmain, xlab="Length (Nucleotides)")
# Calculations: (used with wobble = 0, offset = 10000)
# proportion of duplicates in 5' and 3' positions
sum(flanking$before.bool) / nrow(flanking)
sum(flanking$after.bool) / nrow(flanking)

# same calculation with a subset of insertion lengths (6-11 nt)
subset <- flanking[flanking$len >= 6 & flanking$len < 12, ]
sum(subset$before.bool) / nrow(subset)
sum(subset$after.bool) / nrow(subset)



# HISTOGRAMS OF NUCLEOTIDE DIFFERENCES STRATIFIED BY INSERTION LENGTH
# ---------------------------------------------
# NEEDS FIXING 
toTest <- data.frame(counts=c(ins$before,ins$after), len=rep(nchar(ins$indel),2))
caxis=1.1
clab=1.3
cmain=1.6
#main=
  
dev.off()
cex=2
par(mfrow=c(4,2), xpd=NA, mar=c(4,6,4,5),las=0)
#hist(toTest[toTest$len==1,1], col="red", ylim=c(0,4),main="Indel Length: 1", xlab="Number of Differences ",cex.lab=clab, cex.axis=caxis, cex.main=cmain)
hist(toTest[toTest$len==3,1], breaks=seq(-0.5,max(toTest[toTest$len==3,1], na.rm=T)+0.5),col="red", main="Ins Length: 3", xlab="Number of Differences ",cex.lab=clab, cex.axis=caxis, cex.main=cmain)
hist(toTest[toTest$len==4,1], breaks=seq(-0.5,max(toTest[toTest$len==4,1], na.rm=T)+0.5),col="red", ylim=c(0,2), main="Ins Length: 4", xlab="Number of Differences ",cex.lab=clab, cex.axis=caxis, cex.main=cmain)
hist(toTest[toTest$len==6,1], breaks=seq(-0.5,max(toTest[toTest$len==6,1], na.rm=T)+0.5),col="red",  main="Ins Length: 6", xlab="Number of Differences ",cex.lab=clab, cex.axis=caxis, cex.main=cmain)
hist(toTest[toTest$len==7,1], breaks=seq(-0.5,max(toTest[toTest$len==7,1], na.rm=T)+0.5),col="red", ylim=c(0,4), main="Ins Length: 7", xlab="Number of Differences ",cex.lab=clab, cex.axis=caxis, cex.main=cmain)
hist(toTest[toTest$len==9,1], breaks=seq(-0.5,max(toTest[toTest$len==9,1], na.rm=T)+0.5),col="red", main="Ins Length: 9", xlab="Number of Differences ",cex.lab=clab, cex.axis=caxis, cex.main=cmain)
hist(toTest[toTest$len==12,1], breaks=seq(-0.5,max(toTest[toTest$len==12,1], na.rm=T)+0.5),col="red", main="Ins Length: 12", xlab="Number of Differences ",cex.lab=clab, cex.axis=caxis, cex.main=cmain)
hist(toTest[toTest$len==15,1],breaks=seq(-0.5,max(toTest[toTest$len==15,1], na.rm=T)+0.5), col="red", main="Ins Length: 15", xlab="Number of Differences ",cex.lab=clab, cex.axis=caxis, cex.main=cmain)
hist(toTest[toTest$len>15,1], breaks=seq(-0.5,max(toTest[toTest$len>15,1], na.rm=T)+0.5),col="red", main="Ins Length: > 15", xlab="Number of Differences ",cex.lab=clab, cex.axis=caxis, cex.main=cmain)

# SET A CUTOFF 
# negates any results that are below a certain length threshold (I chose 6 nt as the minimum)
ins[which(nchar(ins$indel)< 6),"before"] <- NA
ins[which(nchar(ins$indel)< 6),"after"] <- NA

# proportion of insertions that have a match within x nucl before, and after the insertion
length(ins[!is.na(ins$before), 'before'])/nrow(ins[which(nchar(ins$indel) >=6),])
length(ins[!is.na(ins$after), 'after'])/nrow(ins[which(nchar(ins$indel) >=6),])
