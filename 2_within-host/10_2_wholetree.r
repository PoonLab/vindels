require(bbmle)
require(stringr)
require(ape)
require(phangorn)
require(data.table)
source("~/vindels/2_within-host/utils.r")

path <- "~/PycharmProjects/hiv-withinhost/"

ifolder <- Sys.glob(paste0(path,"9Indels/supp/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"9Indels/supp/del/*.tsv"))

# ifolder <- ifolder[grepl(newreg,ifolder)]
# dfolder <- dfolder[grepl(newreg,dfolder)]


ins.sep <-list()
del.sep <- list()
count <- 0

ins.nosep <- list()
del.nosep <- list()

ins.final <- list()
del.final <- list()

for (file in 1:length(ifolder)){
  print(file)
  filename <- basename(ifolder[file])
  
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F, sep="\t")
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F, sep="\t")

  # extracts info from the indel column and puts it into two separate columns
  insInfo <- t(unname(sapply(iCSV$indel, extractInfo)))
  insInfo <- as.data.frame(insInfo)
  insInfo$indel <- as.character(insInfo$V1)
  insInfo$pos <- as.character(insInfo$V2)
  iCSV <- cbind(iCSV, insInfo[,c('indel','pos')])
  iCSV$indel <- NULL
  
  delInfo <- t(unname(sapply(dCSV$indel, extractInfo)))
  delInfo <- as.data.frame(delInfo)
  delInfo$indel <- as.character(delInfo$V1)
  delInfo$pos <- as.character(delInfo$V2)
  dCSV <- cbind(dCSV, delInfo[,c('indel','pos')])
  dCSV$indel <- NULL
  
  if(all(is.na(iCSV$indel))){
    iCSV$indel <- ""
    iCSV$pos <- ""
  }
  if(all(is.na(dCSV$indel))){
    dCSV$indel <- ""
    dCSV$pos <- ""
  }
  

  
  # Load time-based branch lengths from the time-scaled trees
  tre <- read.tree(paste0(path,"7_5_MCC/main/prelim/", strsplit(filename, "\\.tsv")[[1]], ".tree.sample"))

  res <- unname(sapply(iCSV$header, findAncestor, tree=tre)) 
  
  rttlens <- node.depth.edgelength(tre)
  midpoints <- (rttlens[res] + rttlens[tre$edge[match(res, tre$edge[,2]),1]]) / 2
  
  iCSV$length <- midpoints # tre$edge.length[match(res, tre$edge[,2])]
  dCSV$length <- iCSV$length
  
  iCSV$count <- unname(sapply(iCSV$indel, csvcount))
  dCSV$count <- unname(sapply(dCSV$indel, csvcount))
  
  iCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(iCSV))
  dCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(dCSV))
  
  # iCSV$header <- unname(mapply(labels, iCSV$header, iCSV$pat))
  # dCSV$header <- unname(mapply(labels, dCSV$header, dCSV$pat))
  
  res <- as.data.frame(t(unname(sapply(iCSV$pat, extractPat))))
  colnames(res) <- c("pat", "rep")
  
  iCSV <- cbind(iCSV, res)
  dCSV <- cbind(dCSV, res)
  
  iCSV$pat <- NULL
  dCSV$pat <- NULL
  
  iCSV <- iCSV[,c(10,11,2,3,8,9,6,7,4,5)]
  dCSV <- dCSV[,c(10,11,2,3,8,9,6,7,4,5)]
  
  # Store DF with no separations
  ins.nosep[[file]] <-  iCSV
  del.nosep[[file]] <-  dCSV
  
  # COMMA SEPARATION FIX
  # make a new data.frame for each CSV df
  # transport over all rows which do NOT contain a comma
  iTemp <- iCSV[!grepl(",",iCSV$indel),]
  dTemp <- dCSV[!grepl(",",dCSV$indel),]
  
  # handle comma rows separately with a function
  iCommas <- iCSV[grepl(",",iCSV$indel),]
  dCommas <- dCSV[grepl(",",dCSV$indel),]
  
  # APPLY THE SPLIT ROWS TO GET ONE INDEL PER ROW
  colindel <- which(colnames(iCSV)=='indel')
  colpos <-  which(colnames(iCSV)=='pos')
  if (nrow(iCommas) > 0){
    newrows <- apply(iCommas,1,splitRows, c(colindel,colpos))
    for (i in 1:length(newrows)){
      idx <- as.double(names(newrows)[i])
      len <- nrow(newrows[[i]])
      rownames(newrows[[i]]) <- seq(0,0.1*(len-1),length=len) + idx
      iTemp <- rbind(iTemp, newrows[[i]])
    }
  }
  if (nrow(dCommas) > 0){
    newrows <- apply(dCommas,1,splitRows, c(colindel,colpos))
    for (i in 1:length(newrows)){
      idx <- as.double(names(newrows)[i])
      len <- nrow(newrows[[i]])
      rownames(newrows[[i]]) <- seq(0,0.1*(len-1),length=len) + idx
      dTemp <- rbind(dTemp, newrows[[i]])
    }
  }
  
  # Add the V position column into the two final data frames 
  #iTemp$vpos <- as.numeric(unname(mapply(addPos, pos=iTemp$pos, header=iTemp$header, vloop=iTemp$vloop)))
  #dTemp$vpos <- as.numeric(unname(mapply(addPos, pos=dTemp$pos, header=dTemp$header, vloop=dTemp$vloop)))
  
  # OUTPUT 
  # for other analyses
  # -----------------------------
  ins.sep[[file]] <- iTemp
  del.sep[[file]] <- dTemp
}

# CLEAN UP
rm(iTemp)
rm(dTemp)
rm(iCSV)
rm(dCSV)
rm(iCommas)
rm(dCommas)

# FINALIZE DATAFRAMES 

ins.sep <- as.data.frame(rbindlist(ins.sep))
del.sep <- as.data.frame(rbindlist(del.sep))

ins.nosep <- as.data.frame(rbindlist(ins.nosep))
del.nosep <- as.data.frame(rbindlist(del.nosep))

ins <- ins.sep[ins.sep$indel!="",]
del <- del.sep[del.sep$indel!="",]


# ---- Data Cleaning ----
iprop <- nchar(ins$indel) / nchar(ins$anc)
dprop <- nchar(del$indel) / nchar(del$anc)

#ins <- ins[-which(iprop >= 1),]
del <- del[-which(dprop >= 1),]

# Remove deletions that are greater than 100 nt
ins <- ins[-which(nchar(ins$indel) > 80),]
del <- del[-which(nchar(del$indel) > 80),]

#ins <- ins[-which(as.numeric(ins$pos) > nchar(ins$anc)),]
del <- del[-which(as.numeric(del$pos) > nchar(del$anc)),]



#------ Checkpoint: 10_2_sample.RData

# N - GLYC SITE OUTPUTS 
# ---------------------------------------------
write.table(ins,paste0(path, "13_nglycs/all/ins-sep.csv"), row.names=F, sep="\t", quote=F)
write.table(del,paste0(path, "13_nglycs/all/del-sep.csv"), row.names=F, sep="\t", quote=F)

# INDEL LENGTHS OUTPUT 
# ---------------------------------------------

write.csv(ins, paste0(path,"12_lengths/all/ins-all.csv"))
write.csv(del, paste0(path,"12_lengths/all/del-all.csv"))


# ---- Indel Nucleotide Analysis ----
# focus on only the columns needed 
nt.ins <- ins[,c(2,3,4,7,8,10)]
nt.del <- del[,c(2,3,4,7,8,10)]

nt.ins$pos <- as.numeric(nt.ins$pos)
nt.del$pos <- as.numeric(nt.del$pos)

nt.ins$len <- nchar(nt.ins$indel)
nt.del$len <- nchar(nt.del$indel)

nt.ins$vlen <- as.numeric(nt.ins$vlen)
nt.del$vlen <- as.numeric(nt.del$vlen)

nt.del$anc <- paste0(substr(nt.del$anc, 0, nt.del$pos - nchar(nt.del$indel)), substr(nt.del$anc, nt.del$pos+1, nchar(nt.del$anc)))

nt.ins$anc <- gsub("-","",nt.ins$anc)
nt.del$anc <- gsub("-","",nt.del$anc)


# DINUCLEOTIDE PROPORTIONS OUTPUT 
# ------------------------------------

write.csv(nt.ins[nt.ins$len>1,], paste0(path,"10_nucleotide/all/dinucl-ins.csv"))
write.csv(nt.del[nt.del$len>1,], paste0(path,"10_nucleotide/all/dinucl-del.csv"))

# # FLANKING INSERTIONS (14) OUTPUT 
# # ------------------------------------
# write.csv(ins, paste0(path,"/10_nucleotide/all/flanking/ins-sep.csv"))
# write.csv(del, paste0(path,"/10_nucleotide/all/flanking/del-sep.csv"))
# 
# write.csv(ins.nosep, paste0(path,"/10_nucleotide/all/flanking/ins-nosep-all.csv"))
# write.csv(del.nosep, paste0(path,"/10_nucleotide/all/flanking/del-nosep-all.csv"))
# 
# # --- Modeling 2 ----
# write.csv(ins.sep, paste0(path,"/10_nucleotide/all/flanking/ins-sep-all.csv"))
# write.csv(del.sep, paste0(path,"/10_nucleotide/all/flanking/del-sep-all.csv"))



# NT PROPORTIONS -- ALL
# --------------------------------------------
nucl <- c("A","C","G","T")
ntprop <- function(df){
  
  indels <- df$indel
  anc <- df$anc
  
  props <- c()
  vprops <- c()
  counts <- c()
  totals <- c(sum(unname(sapply(indels, nchar))), sum(unname(sapply(anc, nchar))))
  
  for (n in 1:4){
    counts[n] <- sum(str_count(indels, nucl[n]))
    props[n] <- counts[n] / totals[1]
    vprops[n] <- sum(str_count(anc, nucl[n])) / totals[2]
  }
  return(list(indel=props, vloop=vprops, count=counts))
}


# ---- FULL ----
ires <- ntprop(nt.ins)
dres <- ntprop(nt.del)


ins.nt <- data.frame(nt=nucl,props=ires$indel,vprops=ires$vloop, count=ires$count)
del.nt <- data.frame(nt=nucl,props=dres$indel,vprops=dres$vloop, count=dres$count)
indel.nt <- rbind(ins.nt, del.nt)
indel.nt$indel <- c(rep(0,4),rep(1,4))


# -----------------------------------------
# RANDOMIZATION TEST (gave negative result)
sampleString <- function(len, vloop){
  len <- len-1
  idx <- sample(1:(nchar(vloop)-len),100, replace=TRUE)
  strings <- sapply(idx, function(x){substr(vloop, x, x+len)})
  a <- mean(unname(sapply(strings, function(x){str_count(x, "A")/nchar(x)})))
  c <- mean(unname(sapply(strings, function(x){str_count(x, "C")/nchar(x)})))
  g <- mean(unname(sapply(strings, function(x){str_count(x, "G")/nchar(x)})))
  t <- mean(unname(sapply(strings, function(x){str_count(x, "T")/nchar(x)})))
  list(a,c,g,t)
}


# generates the randomly sampled substrings for each indel
isample <- mapply(sampleString, nt.ins[,"len"], nt.ins[,"anc"])
dsample <- mapply(sampleString, nt.del[,"len"], nt.del[,"anc"])


# compares the observed proportion to the overall distribution of each nucleotide 
isign <- c()
dsign <- c()
for (i in 1:4){
  idist <- unlist(isample[i,])
  ddist <- unlist(dsample[i,])
  print(i)
  
  if(any(is.na(idist))){
    #print(i)
    idist <- idist[-which(is.na(idist))]
  }
  if(any(is.na(ddist))){
    #print(i)
    ddist <- ddist[-which(is.na(ddist))]
  }
  
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
  print(iQT[[1]])
  print(dQT[[1]])
}

ins.nt$sign <- isign
del.nt$sign <- dsign

# BOOTSTRAPS (95% Confidence Intervals)
# ---------------- 



proportions <- function(df1, df2){
  total1 <- c(sum(nchar(df1[,'indel'])), sum(nchar(df1[,'anc'])))
  total2 <- c(sum(nchar(df2[,'indel'])), sum(nchar(df2[,'anc'])))
  
  nucl <- c("A","C","G","T")
  result <- unname(sapply(nucl, function(nuc){
    iProps <-  sum(str_count(df1$indel, nuc)) / total1[1]
    dProps <-  sum(str_count(df2$indel, nuc)) / total2[1]
    
    iVProps <-  sum(str_count(df1$anc, nuc)) / total1[2]
    dVProps <-  sum(str_count(df2$anc, nuc)) / total2[2]
    
    return(c(iProps, dProps, iVProps, dVProps))
  }))
  colnames(result) <- nucl
  return(result)
}

bs.props <- list()
for (n in 1:100){
  
  # create the bootstrap sample 
  sam1 <- sample(nrow(nt.ins), nrow(nt.ins), replace=T)
  sam2 <- sample(nrow(nt.del), nrow(nt.del), replace=T)
  
  df1.bs <- nt.ins[sam1,]
  df2.bs <- nt.del[sam2,]
  
  props <- proportions(df1.bs, df2.bs)
  
  for (n in 1:4){
    bs.props[[paste0("ins-",n)]] <- c(bs.props[[paste0("ins-",n)]], unname(props[1,n]))
    bs.props[[paste0("del-",n)]] <- c(bs.props[[paste0("del-",n)]], unname(props[2,n]))
    bs.props[[paste0("v-ins-",n)]] <- c(bs.props[[paste0("v-ins-",n)]], unname(props[3,n]))
    bs.props[[paste0("v-del-",n)]] <- c(bs.props[[paste0("v-del-",n)]], unname(props[4,n]))
  }
}

m.tab <- t(proportions(nt.ins, nt.del))
colnames(m.tab) <- c("ins","del","v-ins","v-del")

med.x <- as.vector(m.tab[,c(3,4)])
med.y <- as.vector(m.tab[,c(1,2)])

means <- as.vector(t(m.tab))
bounds <- lapply(1:length(bs.props), function(i){
  quantile(bs.props[[i]] - means[i], c(0.025,0.975))
})
con.int <- t(unname(sapply(1:length(bounds), function(j){
  c(means[j] - bounds[[j]][[2]], means[j] - bounds[[j]][[1]])
})))

cols <- rep(1:4,4, each=2)
rows <- rep(1:8,each=4)

lower.x <- con.int[c(3,4,7,8,11,12,15,16),1]
upper.x <- con.int[c(3,4,7,8,11,12,15,16),2]

lower.y <- con.int[c(1,2,5,6,9,10,13,14),1]
upper.y <- con.int[c(1,2,5,6,9,10,13,14),2]
# lower.x <- con.int[which(grepl("v",names(con.int)) & grepl("2.5",names(con.int)))]
# upper.x <- con.int[which(grepl("v",names(con.int)) & grepl("97.5",names(con.int)))]
# lower.y <- con.int[which(!grepl("v",names(con.int)) & grepl("2.5",names(con.int)))]
# upper.y <- con.int[which(!grepl("v",names(con.int)) & grepl("97.5",names(con.int)))]




require(reshape2)
# ---- SPLIT BY REPLICATE ---- 
ins.rep <- split(nt.ins, nt.ins$rep)
del.rep <- split(nt.del, nt.del$rep)

ins.rep <- t(sapply(ins.rep, ntprop))
del.rep <- t(sapply(del.rep, ntprop))

iprop <- do.call(rbind, ins.rep[,1])
ivprop <- do.call(rbind, ins.rep[,2])
icount <- do.call(rbind, ins.rep[,3])

dprop <- do.call(rbind, del.rep[,1])
dvprop <- do.call(rbind, del.rep[,2])
dcount <- do.call(rbind, del.rep[,3])

colnames(iprop) <- nucl
colnames(ivprop) <- nucl
colnames(dprop) <- nucl
colnames(dvprop) <- nucl
colnames(icount) <- nucl
colnames(dcount) <- nucl

iprop <- reshape2::melt(iprop)
ivprop <- reshape2::melt(ivprop)
dprop <- reshape2::melt(dprop)
dvprop <- reshape2::melt(dvprop)
icount <- reshape2::melt(icount)
dcount <- reshape2::melt(dcount)

iprop <- iprop[order(iprop$Var2, iprop$Var1),]
ivprop <- ivprop[order(ivprop$Var2, ivprop$Var1),]
dprop <- dprop[order(dprop$Var2, dprop$Var1),]
dvprop <- dvprop[order(dvprop$Var2, dvprop$Var1),]
icount <- icount[order(icount$Var2, icount$Var1),]
dcount <- dcount[order(dcount$Var2, dcount$Var1),]


ins.nt <- data.frame(nt=iprop$Var2,props=iprop$value,vprops=ivprop$value, count=as.numeric(icount$value))
del.nt <- data.frame(nt=iprop$Var2,props=dprop$value,vprops=dvprop$value, count=as.numeric(dcount$value))
indel.nt2 <- rbind(ins.nt, del.nt)
indel.nt2$indel <- c(rep(0,80),rep(1,80))


# NT ALL INDELS PLOT 
# broken down by variable loop and nucleotide
# -------------------------------------

colors <- c( "limegreen","dodgerblue","red", "purple")
require(scales)

par(pty="s", xpd=F, mar=c(6,8,2,1),las=0)
clab = 1.9
ctext = 1.7


lim = c(0.15,0.42)
plot(x=med.x, y=med.y, pch=indel.nt[,5]+1, col=rep(colors,2),xlim=lim,ylim=lim,cex=0.03*sqrt(indel.nt$count),
     cex.lab=clab, cex.axis=1.6,lwd=6, ylab='', xlab='',las=1)#, main="Nucleotide Proportions")
title(ylab="Proportion In Indels", line=5,cex.lab=clab)
title(xlab="Proportion in Variable Loops", line=4,cex.lab=clab)

points(x=indel.nt2$vprops, y=indel.nt2$props, pch=1, col=alpha(rep(rep(colors,each=20), times=2),0.6), ylab='', xlab='',las=1)#, main="Nucleotide Proportions")

sqn <- c(seq(1,8,2),seq(1,8,2)+1)

abline(0,1)
xpos <- c(0.18,0.162,0.28,0.375)
ypos <- c(0.16, 0.215, 0.245, 0.395)
# y error bars 
arrows(med.x, lower.y[sqn], med.x, upper.y[sqn], length=0, angle=90, code=3,lwd=2)
#arrows(0.175,0.153,0.167,0.162, length=0, lwd=1.2)
# x error bars 
arrows(lower.x[sqn], med.y, upper.x[sqn], med.y, length=0, angle=90, code=3,lwd=2)
legend(0.368,0.24,legend=nucl, pch=21, cex=1.9, pt.lwd=1, pt.bg=colors, col='black',x.intersp = 1.6,y.intersp=1.1, pt.cex=4)
text(xpos, ypos, cex=ctext, labels=c("C", "G", "T", "A"))
#legend(0.14,0.42,legend=c("Insertions", "Deletions"), pch=c(1,2),cex=1.3, lwd=2, col="black",x.intersp = 1.0,y.intersp=1.3, pt.cex=3)
pos <- c(0.36,0.24)
rect(pos[1]+0.008 , pos[2]+0.003,pos[1]+0.058,pos[2]+0.053)
text(pos[1]+0.043, pos[2]+0.04, labels="Ins", cex=ctext)
text(pos[1]+0.043, pos[2]+0.016, labels="Del", cex=ctext)
points(c(pos[1]+0.021,pos[1]+0.021), c(pos[2]+0.04, pos[2]+0.016), pch=c(1,2), cex=4, lwd=5, col='black', bg='black')






# NT PROPORTIONS BY VARIABLE LOOP 
# --------------------------------------------
vloops <- c(1, 2,3,4,5)
vloops2 <- c("V1","V2","V3","V4","V5")
ins.props <- data.frame()
del.props <- data.frame()

for (i in c(1,2,3,4,5)){
  iTemp <- nt.ins[nt.ins$vloop==i,]
  dTemp <- nt.del[nt.del$vloop==i,]

  iProps <- c()
  dProps <- c()
  
  iVProps <- c()
  dVProps <- c()
  
  # a vector of two totals
  # iTotals[1] = total number of nucleotides in insertion sequences
  # iTotals[2] = total nuimebr of nucleotides in the vloops associated with insertion sequences
  iTotals <- c(sum(unname(sapply(iTemp$indel, nchar))), sum(unname(sapply(iTemp$anc, nchar))))
  dTotals <- c(sum(unname(sapply(dTemp$indel, nchar))),sum(unname(sapply(dTemp$anc, nchar))))
  
  for (nuc in nucl){
    iProps <- c(iProps, sum(str_count(iTemp$indel, nuc)) / iTotals[1])
    dProps <- c(dProps, sum(str_count(dTemp$indel, nuc)) / dTotals[1])
    
    iVProps <- c(iVProps, sum(str_count(iTemp$anc, nuc)) / iTotals[2])
    dVProps <- c(dVProps, sum(str_count(dTemp$anc, nuc)) / dTotals[2])
  }
  ins.props <- rbind(ins.props, data.frame(nt=nucl, iprops=iProps, vprops=iVProps, vloop=rep(vloops[i],4)))
  del.props <- rbind(del.props, data.frame(nt=nucl, dprops=dProps, vprops=dVProps, vloop=rep(vloops[i],4)))
}




# NT PROP INSERTIONS PLOT 
# broken down by variable loop and nucleotide
# -------------------------------------
require(RColorBrewer)
colors <- brewer.pal(4, "Set1")

cex=2
par(pty="s", xpd=F, mar=c(6,6,2,1),las=1)

lim = c(0,0.55)
plot(NA,xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='')
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
abline(0,1)
points(ins.props[,c(3,2)], pch=ins.props[,4]+20, bg=ins.props[,1],cex=3.5)
title(ylab="Proportion In Insertions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.47,0.20,legend=nucl, pch=21,cex=1.5, pt.bg=ins.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.36,0.20,legend=vloops2, pch=c(21,22,23,24,25),cex=1.5, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=NA)
text(-0.135, 0.55, labels="a)",cex=2)



require(RColorBrewer)
colors <- brewer.pal(4, "Set1")
colors[1] <- 'white'

cex=2
par(pty="s", xpd=F, mar=c(6,6,2,1),las=1)

lim = c(0.1,0.5)
plot(NA,xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='')
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
abline(0,1)
points(del.props[,c(3,2)], pch=del.props[,4]+20, bg=del.props[,1],cex=3.5)
title(ylab="Proportion In Deletions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.44,0.25,legend=nucl, pch=21,cex=1.5, pt.bg=del.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.36,0.25,legend=vloops2, pch=c(21,22,23,24,25),cex=1.5, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=NA)
text(0, 0.5, labels="b)",cex=2)





#TRASH
# -----------------
cex=2
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0.10,0.45)
plot(del.props[,c(3,2)], pch=21, bg=colors,xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Deletions")
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Inside Deletions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.38,0.22,legend=nucl, pch=21,cex=1.9, pt.bg=colors,x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=F)
abline(0,1)

require(ggplot2)

plot <- ggplot(ntcount, aes(x=nt, 
                               y=props,
                               width=1)) + geom_bar(colour="black",
                                                    stat="identity", 
                                                    fill="dodgerblue",
                                                    position="dodge", 
                                                    show.legend=F) 

plot <- plot + labs(x="Nucleotide", 
                    y="Proportion in Insertions")+scale_y_continuous(expand = c(0, 0),
                                                                                                                                                     limits = c(0, 1))+theme(panel.grid.major.y = element_line(color="black",size=0.3),
                                                                                                                                                                               panel.grid.major.x = element_blank(),
                                                                                                                                                                               panel.grid.minor.y = element_blank(),
                                                                                                                                                                               panel.grid.minor.x = element_blank(),
                                                                                                                                                                               panel.spacing=unit(1, "mm"),
                                                                                                                                                                               #panel.background=element_rect(fill="gray88",colour="white",size=0),
                                                                                                                                                                               plot.margin =margin(t = 20, r = 20, b = 20, l = 8, unit = "pt"),
                                                                                                                                                                               axis.line = element_line(colour = "black"), 
                                                                                                                                                                               axis.title.y=element_text(size=20,margin=margin(t = 0, r = 15, b = 0, l = 0)),
                                                                                                                                                                               axis.title.x=element_text(size=20,margin=margin(t = 15, r = 0, b = 0, l = 0)),
                                                                                                                                                                               strip.text.x = element_text(size=16),
                                                                                                                                                                               axis.text=element_text(size=14),
                                                                                                                                                                               legend.position="none")




# RANDOMIZATION TEST 


