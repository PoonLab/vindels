require(bbmle)
require(stringr)
require(ape)
require(phangorn)
require(data.table)
source("~/vindels/2_within-host/utils.r")
# Lio
path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
 
ifolder <- Sys.glob(paste0(path,"9Indels/mcc/wholetree/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"9Indels/mcc/wholetree/del/*.tsv"))

all.ins <-list()
all.del <- list()
count <- 0

ins.nosep <- list()
del.nosep <- list()

for (file in 1:length(ifolder)){
  print(file)
  filename <- basename(ifolder[file])
  
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F, sep="\t")
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F, sep="\t")

  # extracts info from the indel column and puts it into two separate columns
  insInfo <- sapply(iCSV$indel, extractInfo)
  insInfo <- unname(insInfo)
  insInfo <- t(insInfo)
  insInfo <- as.data.frame(insInfo)
  insInfo$V1 <- as.character(insInfo$V1)
  insInfo$V2 <- as.character(insInfo$V2)
  iCSV <- cbind(iCSV, insInfo)
  iCSV$indel <- NULL
  
  delInfo <- sapply(dCSV$indel, extractInfo)
  delInfo <- unname(delInfo)
  delInfo <- t(delInfo)
  delInfo <- as.data.frame(delInfo)
  delInfo$V1 <- as.character(delInfo$V1)
  delInfo$V2 <- as.character(delInfo$V2)
  dCSV <- cbind(dCSV, delInfo)
  dCSV$indel <- NULL
  
  if(all(is.na(iCSV$V1))){
    iCSV$V1 <- ""
    iCSV$V2 <- ""
  }
  
  
  iCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(iCSV))
  dCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(dCSV))
  
  # Load time-based branch lengths from the time-scaled trees
  tre <- read.tree(paste0(path,"7_5_MCC/prelim/", strsplit(filename, "\\.tsv")[[1]], ".tree"))

  # remove the root from the rtt length vector because it is NOT found in the reconstruction or the indel extraction (deprecated)
  #lens <- lens[-(Ntip(tre)+1)]
  res <- unname(sapply(iCSV$header, function(x){
    # this expression will return results for NODES ONLY
    # second column provides the CAPTURED TIP LABELS from within the node label
    x <- substr(x, 1, nchar(x)-2)
    tips <- str_match_all(x,"([^\\)\\(,\n:]+):")[[1]][,2]
    if (length(tips) == 0){
      # no colons; this means its a TIP 
      # the index in the tre$tip.label vector is the final result
      index <- match(x, tre$tip.label)
    }else{
      # retreive all descendants of every node and tip in the tree
      desc <- Descendants(tre)
      
      # find the numeric labels of all extracted tips 
      matches <- match(tips, tre$tip.label)
      
      # find the SINGLE node in the descendants list that contains the exact same subset of tips
      index <- which(sapply(desc, function(x){ifelse(length(x) == length(matches) && all(x==matches),T,F)}))
    }
    if (length(index)!=1){
      return(paste0("PROBLEM:",as.character(index)))
    }
    return(c(index))
  }))
  
  iCSV$length <- tre$edge.length[match(res, tre$edge[,2])]
  dCSV$length <- iCSV$length
  
  # lens <- node.depth.edgelength(tre)
  # iCSV$rtt.mid <- (lens[res] + lens[tre$edge[match(res, tre$edge[,2]),1]]) / 2
  # dCSV$rtt.mid <- iCSV$rtt.mid 
  
  iCSV$count <- unname(sapply(iCSV$V1, csvcount))
  dCSV$count <- unname(sapply(dCSV$V1, csvcount))
  
  iCSV$header <- unname(mapply(labels, iCSV$header, iCSV$pat))
  dCSV$header <- unname(mapply(labels, dCSV$header, dCSV$pat))
  
  iCSV <- iCSV[,c(1,2,3,9,10,6,7,4,5,8)]
  dCSV <- dCSV[,c(1,2,3,9,10,6,7,4,5,8)]
  
  colnames(iCSV) <- c("header","vloop", "vlen", "length","count", "indel", "pos", "tip","anc","pat")
  colnames(dCSV) <- c("header", "vloop", "vlen", "length","count",  "indel", "pos", "tip","anc","pat")
  
  ins.nosep[[file]] <-  iCSV
  del.nosep[[file]] <-  dCSV
  # COMMA SEPARATION FIX
  # make a new data.frame for each CSV df
  # transport over all rows which do NOT contain a comma
  ins.sep <- iCSV[!grepl(",",iCSV$indel),]
  del.sep <- dCSV[!grepl(",",dCSV$indel),]
  
  # handle comma rows separately with a function
  iCommas <- iCSV[grepl(",",iCSV$indel),]
  dCommas <- dCSV[grepl(",",dCSV$indel),]
  # APPLY THE SPLIT ROWS TO GET ONE INDEL PER ROW
  if (nrow(iCommas) > 0){
    newrows <- apply(iCommas,1,splitRows, c(6,7))
    for (i in 1:length(newrows)){
      idx <- as.double(names(newrows)[i])
      len <- nrow(newrows[[i]])
      rownames(newrows[[i]]) <- seq(0,0.1*(len-1),length=len) + idx
      ins.sep <- rbind(ins.sep, newrows[[i]])
    }
  }
  if (nrow(dCommas) > 0){
    newrows <- apply(dCommas,1,splitRows, c(6,7))
    for (i in 1:length(newrows)){
      idx <- as.double(names(newrows)[i])
      len <- nrow(newrows[[i]])
      rownames(newrows[[i]]) <- seq(0,0.1*(len-1),length=len) + idx
      del.sep <- rbind(del.sep, newrows[[i]])
    }
  }
  
  # Add the V position column into the two final data frames 
  #ins.sep$vpos <- as.numeric(unname(mapply(addPos, pos=ins.sep$pos, header=ins.sep$header, vloop=ins.sep$vloop)))
  #del.sep$vpos <- as.numeric(unname(mapply(addPos, pos=del.sep$pos, header=del.sep$header, vloop=del.sep$vloop)))
  
  

  # OUTPUT 
  # for other analyses
  # -----------------------------
  all.ins[[file]] <- ins.sep
  all.del[[file]] <- del.sep

}
all.ins <- as.data.frame(rbindlist(all.ins))
all.del <- as.data.frame(rbindlist(all.del))

ins.nosep <- as.data.frame(rbindlist(ins.nosep))
del.nosep <- as.data.frame(rbindlist(del.nosep))

ins <- all.ins[all.ins$indel!="",]
del <- all.del[all.del$indel!="",]

# N - GLYC SITE OUTPUTS 
# ---------------------------------------------
write.table(ins,paste0(path, "13_nglycs/all/ins-sep.csv"), row.names=F, sep="\t", quote=F)
write.table(del,paste0(path, "13_nglycs/all/del-sep.csv"), row.names=F, sep="\t", quote=F)

# INDEL LENGTHS OUTPUT 
# ---------------------------------------------

write.csv(all.ins[,c(1,2,5,6)], paste0(path,"12_lengths/all/ins-all.csv"))
write.csv(all.del[,c(1,2,5,6)], paste0(path,"12_lengths/all/del-all.csv"))


# ---- Indel Nucleotide Analysis ----
# focus on only the columns needed 
total.ins <- ins[,c(1,2,6,8)]
total.del <- del[,c(1,2,6,8)]


total.ins$len <- sapply(total.ins$indel, nchar)
total.del$len <- sapply(total.del$indel, nchar)


write.csv(total.ins, paste0(path,"10_nucleotide/all/total-ins.csv"))
write.csv(total.del, paste0(path,"10_nucleotide/all/total-del.csv"))

# DINUCLEOTIDE PROPORTIONS OUTPUT 
# ------------------------------------
total.ins2 <- total.ins[total.ins$len>1, ]
total.del2 <- total.del[total.del$len>1, ]
write.csv(total.ins2, paste0(path,"10_nucleotide/all/total-ins2.csv"))
write.csv(total.del2, paste0(path,"10_nucleotide/all/total-del2.csv"))

# FLANKING INSERTIONS (14) OUTPUT 
# ------------------------------------
write.csv(ins, paste0(path,"/10_nucleotide/all/ins-sep.csv"))
write.csv(del, paste0(path,"/10_nucleotide/all/del-sep.csv"))

write.csv(ins.nosep, paste0(path,"/10_nucleotide/all/ins-nosep-all.csv"))
write.csv(del.nosep, paste0(path,"/10_nucleotide/all/del-nosep-all.csv"))

# --- Modeling 2 ----
write.csv(all.ins, paste0(path,"/10_nucleotide/all/ins-sep-all.csv"))
write.csv(all.ins, paste0(path,"/10_nucleotide/all/del-sep-all.csv"))



# NT PROPORTIONS -- ALL
# --------------------------------------------

nucleotides <- c("A","C","G","T")
iProps <- c()
dProps <- c()
iVProps <- c()
dVProps <- c()
counts <- data.frame()
iTotals <- c(sum(unname(sapply(total.ins$indel, nchar))), sum(unname(sapply(total.ins$tip, nchar))))
dTotals <- c(sum(unname(sapply(total.del$indel, nchar))),sum(unname(sapply(total.del$tip, nchar))))

for (nuc in nucleotides){
  icount <- sum(str_count(total.ins$indel, nuc))
  dcount <- sum(str_count(total.del$indel, nuc))
  counts <- rbind(counts, data.frame(nucl=nuc, ins=icount, del=dcount))
  
  iProps <- c(iProps, icount / iTotals[1])
  dProps <- c(dProps, dcount / dTotals[1])
  
  iVProps <- c(iVProps, sum(str_count(total.ins$tip, nuc)) / iTotals[2])
  dVProps <- c(dVProps, sum(str_count(total.del$tip, nuc)) / dTotals[2])
}
require(reshape)
counts2 <- reshape2::melt(counts)

ins.nt <- data.frame(nt=nucleotides,props=iProps,vprops=iVProps)
del.nt <- data.frame(nt=nucleotides,props=dProps,vprops=dVProps)
indel.nt <- rbind(ins.nt, del.nt)
indel.nt$indel <- c(rep(0,4),rep(1,4))
indel.nt$counts <- counts2$value

# -----------------------------------------
# RANDOMIZATION TEST (gave negative result)
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



iSample <- list(c(),c(),c(),c())
dSample <- list(c(),c(),c(),c())
# generates the randomly sampled substrings for each indel
for (row in 1:nrow(total.ins)){
  itemp <- sampleString(total.ins[row,"len"], total.ins[row,"tip"])
  for (i in 1:4){
    iSample[[i]] <- c(iSample[[i]], itemp[[i]])
  }
}
for (row in 1:nrow(total.del)){
  dtemp <- sampleString(total.del[row,"len"], total.del[row,"tip"])
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

ins.nt$sign <- isign
del.nt$sign <- dsign

# BOOTSTRAPS
# ----------------

df1 <- total.ins
df2 <- total.del
bs.props <- list()

for (n in 1:1000){
  
  # create the bootstrap sample 
  sam1 <- sample(nrow(df1), nrow(df1), replace=T)
  sam2 <- sample(nrow(df2), nrow(df2), replace=T)
  
  df1.bs <- df1[sam1,c('indel','tip')]
  df2.bs <- df2[sam2,c('indel','tip')]
  
  total1 <- c(sum(unname(sapply(df1.bs[,"indel"], nchar))), sum(unname(sapply(df1.bs[,"tip"], nchar))))
  total2 <- c(sum(unname(sapply(df2.bs[,"indel"], nchar))), sum(unname(sapply(df2.bs[,"tip"], nchar))))
  
  props <- sapply(nucleotides, function(nuc){
    iProps <-  sum(str_count(df1.bs$indel, nuc)) / total1[1]
    dProps <-  sum(str_count(df2.bs$indel, nuc)) / total2[1]
    
    iVProps <-  sum(str_count(df1$tip, nuc)) / total1[2]
    dVProps <-  sum(str_count(df2$tip, nuc)) / total2[2]
    
    return(c(iProps, dProps, iVProps, dVProps))
  })
 
  for (n in nucleotides){
    bs.props[[paste0("ins-",n)]] <- c(bs.props[[paste0("ins-",n)]], unname(props[1,n]))
    bs.props[[paste0("del-",n)]] <- c(bs.props[[paste0("del-",n)]], unname(props[2,n]))
    bs.props[[paste0("v-ins-",n)]] <- c(bs.props[[paste0("v-ins-",n)]], unname(props[3,n]))
    bs.props[[paste0("v-del-",n)]] <- c(bs.props[[paste0("v-del-",n)]], unname(props[4,n]))
  }
}
medians <- unlist(lapply(bs.props,median))
m.tab <- matrix(nrow=4,ncol=4, dimnames=list(nucleotides,c("ins","del","v-ins","v-del")))
cols <- rep(1:4,4)
rows <- rep(1:4,each=4)
for (i in 1:16){
  m.tab[rows[i], cols[i]] <- medians[[i]]
}

med.x <- as.vector(m.tab[,c(3,4)])
med.y <- as.vector(m.tab[,c(1,2)])

con.int <- unlist(lapply(bs.props, function(x){quantile(x, c(0.025,0.975))}))
cols <- rep(1:4,4, each=2)
rows <- rep(1:8,each=4)
lower.x <- con.int[which(grepl("v",names(con.int)) & grepl("2.5",names(con.int)))]
upper.x <- con.int[which(grepl("v",names(con.int)) & grepl("97.5",names(con.int)))]
lower.y <- con.int[which(!grepl("v",names(con.int)) & grepl("2.5",names(con.int)))]
upper.y <- con.int[which(!grepl("v",names(con.int)) & grepl("97.5",names(con.int)))]



# NT ALL INDELS PLOT 
# broken down by variable loop and nucleotide
# -------------------------------------
require(RColorBrewer)
#colors <- brewer.pal(4, 'Set1')
colors <- c( "limegreen","dodgerblue","red", "magenta")

cex=1
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

xpos <- c(0.17,0.155,0.24,0.37)
ypos <- c(0.15, 0.21, 0.26, 0.39)

lim = c(0.14,0.42)
plot(indel.nt[,c(3,2)], pch=indel.nt[,4]+1, col=rep(colors,2),xlim=lim,ylim=lim,cex=0.08*sqrt(indel.nt$counts),
     cex.lab=1.3, cex.axis=1.2,cex.main=1.8,lwd=5, ylab='', xlab='', main="Nucleotide Proportions")
title(ylab="Proportion In Indels", line=3,cex.lab=1.3)
title(xlab="Proportion in Variable Loops", line=3,cex.lab=1.3)
arrows(indel.nt[1:8,3], lower.y[c(seq(1,8,2),seq(1,8,2)+1)], indel.nt[1:8,3], upper.y[c(seq(1,8,2),seq(1,8,2)+1)], length=0, angle=90, code=3,lwd=1.5)
arrows(lower.x[c(seq(1,8,2),seq(1,8,2)+1)], indel.nt[1:8,2], upper.x[c(seq(1,8,2),seq(1,8,2)+1)], indel.nt[1:8,2], length=0, angle=90, code=3,lwd=1.5)
#legend(0.38,0.24,legend=nucleotides, pch=22,cex=1.3, pt.bg=indel.nt[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
text(xpos, ypos, cex=1.3, labels=c("C", "G", "T", "A"))
#legend(0.14,0.42,legend=c("Insertions", "Deletions"), pch=c(1,2),cex=1.3, lwd=2, col="black",x.intersp = 1.0,y.intersp=1.3, pt.cex=3)
par(xpd=F)
abline(0,1)
rect(0.14,0.35,0.21,0.42)
text(0.19, 0.40, labels="Ins", cex=1.5)
text(0.19, 0.37, labels="Del", cex=1.5)
points(c(0.16,0.16), c(0.40, 0.37), pch=c(1,2), cex=3, lwd=5, col='black', bg='black')






# NT PROPORTIONS BY VARIABLE LOOP 
# --------------------------------------------
vloops <- c(1, 2,3,4,5)
vloops2 <- c("V1","V2","V4","V5")
ins.props <- data.frame()
del.props <- data.frame()

for (i in c(1,2,4,5)){
  iTemp <- total.ins[total.ins$vloop==i,]
  dTemp <- total.del[total.del$vloop==i,]

  iProps <- c()
  dProps <- c()
  
  iVProps <- c()
  dVProps <- c()
  
  # a vector of two totals
  # iTotals[1] = total number of nucleotides in insertion sequences
  # iTotals[2] = total nuimebr of nucleotides in the vloops associated with insertion sequences
  iTotals <- c(sum(unname(sapply(iTemp$indel, nchar))), sum(unname(sapply(iTemp$tip, nchar))))
  dTotals <- c(sum(unname(sapply(dTemp$indel, nchar))),sum(unname(sapply(dTemp$tip, nchar))))
  
  for (nuc in nucleotides){
    iProps <- c(iProps, sum(str_count(iTemp$indel, nuc)) / iTotals[1])
    dProps <- c(dProps, sum(str_count(dTemp$indel, nuc)) / dTotals[1])
    
    iVProps <- c(iVProps, sum(str_count(iTemp$tip, nuc)) / iTotals[2])
    dVProps <- c(dVProps, sum(str_count(dTemp$tip, nuc)) / dTotals[2])
  }
  ins.props <- rbind(ins.props, data.frame(nt=nucleotides, iprops=iProps, vprops=iVProps, vloop=rep(vloops[i],4)))
  del.props <- rbind(del.props, data.frame(nt=nucleotides, dprops=dProps, vprops=dVProps, vloop=rep(vloops[i],4)))
}




# NT PROP INSERTIONS PLOT 
# broken down by variable loop and nucleotide
# -------------------------------------
require(RColorBrewer)
colors <- brewer.pal(4, "Set1")

cex=2
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0,0.5)
plot(ins.props[,c(3,2)], pch=ins.props[,4]+20, bg=ins.props[,1],xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Insertions")
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Inside Insertions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.43,0.18,legend=nucleotides, pch=21,cex=1.5, pt.bg=ins.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.33,0.18,legend=vloops2, pch=c(21,22,24,25),cex=1.5, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
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

lim = c(0.1,0.50)
plot(del.props[,c(3,2)], pch=del.props[,4]+20, bg=del.props[,1],xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Deletions")
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Inside Deletions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.44,0.24,legend=nucleotides, pch=21,cex=1.5, pt.bg=del.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.35,0.24,legend=vloops2, pch=c(21,22,24,25),cex=1.5, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=F)
abline(0,1)




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
legend(0.38,0.22,legend=nucleotides, pch=21,cex=1.9, pt.bg=colors,x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
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


