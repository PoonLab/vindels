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

ATcontent <- function(indel, pos, vseq, sample_length){
  pos <- as.numeric(pos)
  nucleotides <- c("A", "C","G","T")
  
  counts <- function(indel,pos,vseq){
    len <- nchar(indel)
    nucleotides <- c("A", "C","G","T")
    
    if (pos - sample_length >= 0){
      upstream <- substr(vseq, pos-sample_length, pos)
      #print(upstream)
      ucounts <- c()
      for (n in 1:4){
        ucounts[n] <- str_count(upstream,nucleotides[n])
      }
    }else{
      ucounts <- c(NA,NA,NA,NA)
    }
    
    if ((pos + len + sample_length) <= nchar(vseq)){
      downstream <- substr(vseq, pos+len+1, pos+len+sample_length)
      dcounts <- c()
      for (n in 1:4){
        dcounts[n] <- str_count(downstream,nucleotides[n])
      }
      
    }else{
      dcounts <- c(NA,NA,NA,NA)
    }
    c(ucounts,dcounts)
  }
  nucl <- unname(mapply(counts, indel, pos, vseq))
  up <- as.data.frame(t(nucl))[,1:4]
  print(nrow(up[is.na(up[,1]),]))
  print(nrow(up))
  up <- up[!is.na(up[,1]),]
  colnames(up) <- nucleotides
  
  down <- as.data.frame(t(nucl))[,5:8]
  print(nrow(down[is.na(down[,1]),]))
  print(nrow(down))
  down <- down[!is.na(down[,1]),]
  colnames(down) <- nucleotides
  
  data.frame(upstream=colSums(up)/(nrow(up)*sample_length),downstream=colSums(down)/(nrow(down)*sample_length))
  
}

# RANDOMIZATION TEST 
# -----------------------------------------
sampleString <- function(vloop, pos){
  idx <- sample(1:(nchar(vloop)-14),100, replace=TRUE)
  strings <- sapply(idx, function(x){substr(vloop, x, x+14)})
  props <- list()
  nucl <- c("A","C","G","T")
  for (n in 1:4){
    props[[n]] <- unname(sapply(strings, function(x){str_count(x, nucl[n])/nchar(x)}))
  }
  props
}

iSample <- list(numeric(nrow(flanking)*100),numeric(nrow(flanking)*100),numeric(nrow(flanking)*100),numeric(nrow(flanking)*100))

# generates the randomly sampled substrings for each indel
for (row in 1:nrow(flanking)){
  itemp <- sampleString(flanking[row,"vseq"])
  for (i in 1:4){
    iSample[[i]][((row-1)*100+1):(row*100)] <- itemp[[i]]
  }
}

# compares the observed proportion to the overall distribution of each nucleotide 
usign <- c()
dsign <- c()
u.pct <-c()
d.pct <-c()
for (i in 1:4){
  idist <- iSample[[i]]
  u.prop <- at.props[i,1]
  d.prop <- at.props[i,2]
  
  # find the percentile of the upstream and downstream AT content
  percentile <- ecdf(idist)
  u.pct[i] <- percentile(u.prop)
  d.pct[i] <- percentile(d.prop)
  
  # check for significance
  iQT <- quantile(idist, probs=c(0.025,0.975))
  
  # highlight significant differences 
  if (u.prop < iQT[[1]]){
    usign[i] <- "lower"
  }else if(u.prop > iQT[[2]]){
    usign[i] <- "higher"
  }else{
    usign[i] <- ""
  }
  
  # highlight significant differences 
  if (d.prop < iQT[[1]]){
    dsign[i] <- "lower"
  }else if(d.prop > iQT[[2]]){
    dsign[i] <- "higher"
  }else{
    dsign[i] <- ""
  }
  
}
print(u.pct)
ins.nt$sign <- isign
del.nt$sign <- dsign



# HIDDEN MARKOV MODEL 
# ---------------------------------

HMM <- function(t.mat, emi.mat, p.init, seq)

  
slips <- function(vseq, pos, len){
  pos <- as.numeric(pos)+1
  len <- as.numeric(len)
  
  chars <- strsplit(vseq,"")[[1]]
  
}



path <- '~/PycharmProjects/hiv-withinhost/'
ins <- read.csv(paste0(path,"10_nucleotide/ins-sep.csv"), stringsAsFactors = F, row.names = 1)
all <- read.csv(paste0(path,"10_nucleotide/ins-all.csv"), stringsAsFactors = F, row.names=1)

# FIXED AND NO LONGER NEEDED 
#ins$Accno <- unname(mapply(labels, ins$Accno, ins$Pat, ins$Vloop))
#all$Accno <- unname(mapply(labels, all$Accno, all$Pat, all$Vloop))


# apply flankCheck 
# parameters can be changed here to get different results 
flanking <- unname(mapply(flankCheck, indel=ins$Seq, pos=ins$Pos, vseq=ins$Vseq, wobble=1/6, offset=0))

# modify flanking data.frame 
flanking <- as.data.frame(t(flanking), stringsAsFactors = F)
flanking <- cbind(ins[,c(1,2,7)], len=nchar(ins$Seq), flanking)
colnames(flanking) <- c("header","vloop","pos", "len", "indel", "vseq","before.bool", "before.offset", "before.diff", "before.seq","after.bool", "after.offset", "after.diff",  "after.seq")
flanking[,"before.bool"] <- as.logical(flanking[,"before.bool"] )
flanking[,"after.bool"] <- as.logical(flanking[,"after.bool"] )
flanking[,"before.offset"] <- as.numeric(flanking[,"before.offset"] )
flanking[,"after.offset"] <- as.numeric(flanking[,"after.offset"] )
flanking[,"before.diff"] <- as.numeric(flanking[,"before.diff"] )
flanking[,"after.diff"] <- as.numeric(flanking[,"after.diff"] )

# used to remove 
all <- all[-as.numeric(rownames(all[all$Vlength==0,])),]

# retrieves insertions that have at least one instance of flanking sequence 
tabs <- table(flanking[flanking$before.bool | flanking$after.bool, "header"])
all$count.flanking <- 0
all[all$Header %in% names(tabs),"count.flanking"] <- tabs


all$new.count <- 0
all[all$count.flanking!=0, "new.count"] <- nchar(all[all$count.flanking!=0,"Seq"])

all$subs <- unname(mapply(subs, all$Vseq, all$Anc))


# MAXIMUM LIKELIHOOD ESTIMATION 
# ------------------------------------------------------

# POISSON DISTRIBUTION

pll <- function(rate, count, len){
  lam <- rate * len 
  res <- -lam + count*log(lam) 
  res <- res[!is.na(res)]
  sum(res)
}
obj.f <- function(rate) -pll(rate, all$new.count, all$Vlength)
mle.result <- bbmle::mle2(obj.f, start=list(rate=1), method = "Brent", lower=1e-12, upper = 1)


# BINOMIAL DISTRIBUTION 
binomll <- function(prob, count, len){
  N <- len
  k <- count
  p <- prob
  chs <- factorial(N) / (factorial(k) * factorial(N - k))
  sum(log(chs) +  k*log(p) +  (N - k)*log(1-p))
}
obj.f2 <- function(prob) -binomll(prob, all$new.count, all$Vlength)
mle.result2 <- bbmle::mle2(obj.f2, start=list(prob=1), method = "Brent", lower = 1e-12, upper=1)


# GEOMETRIC LIKELIHOOD FUNCTION 

geomll <- function(forward, count, N){
  N * log(forward) + (sum(count)-N) * log(1-forward)
}
slips.whole <- c(all$count.flanking, rep(0,sum(nchar(all$Vseq))-nrow(all)))
slips.nt <- c(all$new.count, rep(0,sum(nchar(all$Vseq))-nrow(all)))
slips.nt2 <- slips.nt + 1
objf3 <- function(forward) -geomll(forward, slips.nt2, length(slips.nt2))
mle3 <- bbmle::mle2(objf3, start=list(forward=1), method="Brent" , lower=0.9, upper=1)


count.subs <- function()


# CUSTOM GEOMETRIC LIKELIHOOD FUNCTION 
custom1 <- function(slip, mut, len, subs, N){
  N * log(slip) + sum(count) * log(1-slip) + log(choose(count, subs)) + subs * log(mut) + (counts - subs) * log(1-mut)
}

obj4 <- function(slip, mut) -custom1(slip, mut, slips.nt, length(slips.nt))
mle4 <- bbmle::mle2(obj4, start=list(slip=1, mut=1), method="L-BFGS-B",lower=1e-12, upper=1)

x <- runif(5000,min=1e-7, max=1)
y <- unname(sapply(x, objf3))
plot(x=x, y=y)

obs <- all$new.count
lens <- c()
dists <- list()
for (i in 1:25){
  sim <- unname(sapply(all$Vseq, simulateData))
  diff <- nchar(sim) - nchar(all$Vseq)
  simdiff <- diff[diff!=0]
  lens[i] <- length(simdiff)
  dists[[i]] <- simdiff
}




# proportion of insertions ACROSS ALL VARIABLE LOOPS that contain a match with 1/6 wobble directly 
#   adjacent (EITHER 5' or 3', no offset) 
nrow(flanking[flanking$before.bool | flanking$after.bool,]) / nrow(flanking)

# split 'flanking' into 5 different variable loops 
sub.v <- split(flanking, flanking$vloop)

# calculate the probability of slippage occurring in each of the V-loops 
require(R.utils)
p.slip <- c()
for (i in 1:4){
  p.slip[i] <-  nrow(sub.v[[i]][sub.v[[i]]$before.bool | sub.v[[i]]$after.bool,]) / nrow(sub.v[[i]])
}
p.slip <- insert(p.slip,3,0)


# probability of a slippage event occurring per nucleotide per v-loop (flanking contains a single insertion event per row)
nt.slip <- c()
for (i in 1:4){
  nt.slip[i] <-nrow(sub.v[[i]][sub.v[[i]]$before.bool | sub.v[[i]]$after.bool,]) / sum(nchar(all[all$Vloop==i & all$Seq == "","Vseq"]))
}

# creates a transitional probability matrix for EACH V-LOOP (1,2,4,5 ; not V3) 
all$Vseq <- mapply(addX, all$Vseq, all$Pos)
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






# 3 x 2 HISTOGRAM PLOT 
# -------------------------------------
# set flanking to infinite offset, 0 wobble
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
toTest <- data.frame(counts=c(ins$before,ins$after), len=rep(nchar(ins$Seq),2))
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
ins[which(nchar(ins$Seq)< 6),"before"] <- NA
ins[which(nchar(ins$Seq)< 6),"after"] <- NA

# proportion of insertions that have a match within x nucl before, and after the insertion
length(ins[!is.na(ins$before), 'before'])/nrow(ins[which(nchar(ins$Seq) >=6),])
length(ins[!is.na(ins$after), 'after'])/nrow(ins[which(nchar(ins$Seq) >=6),])
