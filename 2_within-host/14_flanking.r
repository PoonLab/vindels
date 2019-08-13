# flanking insertion sequences 
# trying to develop a model for interstrand jumping of polymerase 
require(ape)
require(stringr)

checkDiff <- function(seq1, seq2){
  if (seq1 == seq2){
    return(NULL)
  }
  
  seq1 <- str_split(seq1, "")[[1]]
  seq2 <- str_split(seq2, "")[[1]]
  
  chars <- rbind(seq1, seq2)
  which(chars[1,]!=chars[2,])
}

insCheck <- function(indel,pos,vseq,wobble, offset=0){
  len <- nchar(indel)
  pos <- as.numeric(pos)
  
  beforeBool <- F
  afterBool <- F
  
  beforeIdx <- NA
  afterIdx <- NA
  
  beforeDiff <- NA
  afterDiff <- NA
  
  beforeSeq <- ""
  afterSeq <- ""
  
  # BEFORE 
  # find the best matching position from 0 to offset 
  lowest <- 1000
  bestIdx <- -1

  for (idx in 0:offset){
    #print(pos)
    #print(len)
    #print(idx)
    if ((pos - len - idx) >= 0){
      # then the PRECEDING position can be checked
      before <- substr(vseq, pos-len-idx+1, pos-idx)
      diffs <- checkDiff(indel, before)
      if (length(diffs) < lowest){
        lowest <- length(diffs)
        bestIdx <- idx
      }
    }
  }
  
  if (bestIdx != -1 & lowest <= wobble){
    beforeBool <- T
    beforeIdx <- bestIdx
    beforeDiff <- lowest
    beforeSeq <- substr(vseq,pos-len-bestIdx+1, pos-bestIdx)
  }

  # AFTER 
  # find the best matching position from 0 to offset 
  lowest <- 1000
  bestIdx <- -1
  for (idx in 0:offset){
    if ((pos + len + idx) <= nchar(vseq)){
      # then the PRECEDING position can be checked
      after <- substr(vseq, pos+idx+1, pos+len+idx)
      #print(after)
      diffs <- checkDiff(indel, after)
      if (length(diffs) < lowest){
        lowest <- length(diffs)
        bestIdx <- idx
      }
    }
  }
  
  if (bestIdx != -1 & lowest <= wobble){
    afterBool <- T
    afterIdx <- bestIdx
    afterDiff <- lowest
    afterSeq <- substr(vseq,pos+bestIdx+1, pos+len+bestIdx)
  }

  c(indel, vseq, as.logical(beforeBool),  as.numeric(beforeIdx),  as.numeric(beforeDiff), beforeSeq, as.logical(afterBool), as.numeric(afterIdx), as.numeric(afterDiff),afterSeq)
}


path <- '~/PycharmProjects/hiv-withinhost/'
ins <- read.csv(paste0(path,"10_nucleotide/ins.csv"), stringsAsFactors = F, row.names = 1)

rownames(ins) <- 1:nrow(ins)




# randomly rearrange the same nucleotides found in the vseq 

null_dist <- c()
for (i in 1:nrow(ins)){
  seq <- str_split(ins$Seq[i], ",")[[1]]
  pos <- str_split(ins$Pos[i], ",")[[1]]
  vseq <- ins$ins.unchanged[i]
  letters <- str_split(vseq, "")[[1]]
  
  for (s in 1:length(seq)){
    if (nchar(seq[s]) > 3){
    px <- as.numeric(pos[s]) - nchar(seq[s]) +1
    for (j in 1:1000){
      rseq <- paste(letters[sample(1:nchar(vseq))], collapse = "")
      res <- gregexpr(seq[s],rseq)[[1]]
      if (res != -1){
        res <- res - px
        null_dist <- c(null_dist, res)
      }
    }
    }
  }
}

test_dist <- c()
findMatch <- function(indel, pos, vseq){
  seq <- str_split(indel, ",")[[1]]
  pos <- str_split(pos, ",")[[1]]
  
  res <- c()
  for (s in 1:length(seq)){
    if (nchar(seq[s]) >= 3){
      px <- as.numeric(pos[s]) - nchar(seq[s]) +1
      #print(seq[s])
      #print(vseq)
      match <- gregexpr(seq[s],vseq)[[1]]
      if (length(match) > 1 | match != -1){
        match <- match - px
        res <- c(res, match)
      }
    }
  }
  res
}

testDist <- unname(unlist(mapply(findMatch, ins$Seq, ins$Pos, ins$Vseq)))

hist(testDist, breaks=seq(min(testDist)-0.5,max(testDist)+0.5), col='red', xlab="Length (nt)")

# find the locations of matches with INDEL in this random vseq

# repeat the process 1000 times and get a null distribution describing the distribution of matches expected by chance
 
# compare this distribution to the 


# TO DO 
# incorporate the 


ins <- read.csv(paste0(path,"10_nucleotide/ins-sep.csv"), stringsAsFactors = F, row.names = 1)

flanking <- unname(mapply(insCheck, indel=ins$Seq, pos=ins$Pos, vseq=ins$Vseq, wobble=1, offset=3))
flanking <- as.data.frame(t(flanking), stringsAsFactors = F)
flanking <- cbind( ins[,c(1,7)], len=nchar(ins$Seq), flanking)
colnames(flanking) <- c("accno","pos", "len", "indel", "vseq","before.bool", "before.offset", "before.diff", "before.seq","after.bool", "after.offset", "after.diff",  "after.seq")
flanking[,"before.bool"] <- as.logical(flanking[,"before.bool"] )
flanking[,"after.bool"] <- as.logical(flanking[,"after.bool"] )
flanking[,"before.offset"] <- as.numeric(flanking[,"before.offset"] )
flanking[,"after.offset"] <- as.numeric(flanking[,"after.offset"] )
flanking[,"before.diff"] <- as.numeric(flanking[,"before.diff"] )
flanking[,"after.diff"] <- as.numeric(flanking[,"after.diff"] )

# ANALYSES 
# infinite offset, 0 wobble
# determine how far you need to go to find a match
a.dist <- flanking[!is.na(flanking$after.offset), "after.offset"]
b.dist <- flanking[!is.na(flanking$before.offset), "before.offset"]
a.dist1 <- flanking[!is.na(flanking$after.offset), "after.offset"]
b.dist1 <- flanking[!is.na(flanking$before.offset), "before.offset"]
a.dist2 <- flanking[!is.na(flanking$after.offset), "after.offset"]
b.dist2 <- flanking[!is.na(flanking$before.offset), "before.offset"]

par(mar=c(5,5,5,2))
caxis=1.3
clab=1.5
cmain=1.8

par(mfrow=c(3,2), xpd=NA, mar=c(4,6,4,5),las=0)
# distribution of how far you need to travel to find an EXACT MATCH (wobble = 0, offset=10000)
hist(a.dist, breaks=seq(-0.5,max(a.dist)+0.5), col='red',cex.lab=clab, main="Distances to next exact match - 3'", cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")
hist(b.dist, breaks=seq(-0.5,max(b.dist)+0.5), col='red',cex.lab=clab, main="Distances to next exact match - 5'",cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")

# distribution of how far you need to travel to find a MATCH WITHIN 1 NT (wobble = 1, offset=10000)
hist(a.dist1, breaks=seq(-0.5,max(a.dist1)+0.5), col='red',cex.lab=clab, main="Distances to next match (1 nt) - 3'", cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")
hist(b.dist1, breaks=seq(-0.5,max(b.dist1)+0.5), col='red',cex.lab=clab, main="Distances to next match (1 nt) - 5'",cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")


# distribution of how far you need to travel to find a MATCH WITHIN 2 NT (wobble = 2, offset=10000)
hist(a.dist2, breaks=seq(-0.5,max(a.dist2)+0.5), col='red',cex.lab=clab, main="Distances to next match (2 nt) - 3'", cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")
hist(b.dist2, breaks=seq(-0.5,max(b.dist2)+0.5), col='red',cex.lab=clab, main="Distances to next match (2 nt) - 5'",cex.axis=caxis, cex.main=cmain, xlab="Distance from Insertion Site (nt)")

for (wobble in 0:1){
  
}

a.dist <- a.dist[!is.na(a.dist)] 

m.before <- matrix(nrow=2, ncol=10)
m.after <- matrix(nrow=2,ncol=10)
for (wobble in 0:1){
  
  for (offset in 0:9){
    flanking <- unname(mapply(insCheck, indel=ins$Seq, pos=ins$Pos, vseq=ins$Vseq, wobble=wobble, offset=offset))
    flanking <- as.data.frame(t(flanking))
    flanking <- cbind( ins[,c(1,7)], len=nchar(ins$Seq), flanking)
    colnames(flanking) <- c("accno","pos", "len", "indel", "vseq","before.bool", "before.offset", "before.diff", "before.seq","after.bool", "after.offset", "after.diff",  "after.seq")    
    flanking$before.bool <- as.logical(flanking$before.bool)
    flanking$after.bool <- as.logical(flanking$after.bool)
    
    
    before.prop <- sum(flanking$before.bool) / nrow(flanking)
    after.prop <- sum(flanking$after.bool) / nrow(flanking)
    
    m.before[wobble+1, offset+1] <- before.prop
    m.after[wobble+1,offset+1] <- after.prop
  }
}

rownames(m.before) <- 0:9
colnames(m.before) <- 0:9
rownames(m.after) <- 0:9
colnames(m.after) <- 0:9

# proportion of duplicates in 5' and 3' positions
sum(!is.na(flanking$before)) / nrow(flanking)
sum(!is.na(flanking$after)) / nrow(flanking)
subset <- flanking[flanking$len >= 6 & flanking$len < 12, ]
sum(!is.na(subset$before)) / nrow(subset)
sum(!is.na(subset$after)) / nrow(subset)



# HISTOGRAMS OF NUCLEOTIDE DIFFERENCES STRATIFIED BY INSERTION LENGTH
# ---------------------------------------------
toTest <- data.frame(counts=c(ins$before,ins$after), len=rep(nchar(ins$Seq),2))
caxis=1.1
clab=1.3
cmain=1.6
main=
  
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
ins[which(nchar(ins$Seq)< 6),"before"] <- NA
ins[which(nchar(ins$Seq)< 6),"after"] <- NA

# proportion of insertions that have a match within x nucl before, and after the insertion
length(ins[!is.na(ins$before), 'before'])/nrow(ins[which(nchar(ins$Seq) >=6),])
length(ins[!is.na(ins$after), 'after'])/nrow(ins[which(nchar(ins$Seq) >=6),])
