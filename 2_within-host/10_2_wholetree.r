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
  insInfo$V1 <- as.character(insInfo$V1)
  insInfo$V2 <- as.character(insInfo$V2)
  iCSV <- cbind(iCSV, insInfo)
  iCSV$indel <- NULL
  
  delInfo <- t(unname(sapply(dCSV$indel, extractInfo)))
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
  tre <- read.tree(paste0(path,"7_5_MCC/new_skygrid/prelim/", strsplit(filename, "\\.tsv")[[1]], ".tree.sample"))

  # remove the root from the rtt length vector because it is NOT found in the reconstruction or the indel extraction (deprecated)
  #lens <- lens[-(Ntip(tre)+1)]
  res <- unname(sapply(iCSV$header, findAncestor)) 
  
  #iCSV$length <- tre$edge.length[match(res, tre$edge[,2])]
  #dCSV$length <- iCSV$length
  
  # lens <- node.depth.edgelength(tre)
  # iCSV$rtt.mid <- (lens[res] + lens[tre$edge[match(res, tre$edge[,2]),1]]) / 2
  # dCSV$rtt.mid <- iCSV$rtt.mid 
  
  iCSV$count <- unname(sapply(iCSV$V1, csvcount))
  dCSV$count <- unname(sapply(dCSV$V1, csvcount))
  
  iCSV$header <- unname(mapply(labels, iCSV$header, iCSV$pat))
  dCSV$header <- unname(mapply(labels, dCSV$header, dCSV$pat))
  
  iCSV <- iCSV[,c(1,2,3,9,10,6,7,4,5,8)]
  dCSV <- dCSV[,c(1,2,3,9,10,6,7,4,5,8)]
  
  colnames(iCSV) <- c("header","vloop", "vlen","count", "indel", "pos", "tip","anc","pat")
  colnames(dCSV) <- c("header", "vloop", "vlen", "count",  "indel", "pos", "tip","anc","pat")
  
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
  if (nrow(iCommas) > 0){
    newrows <- apply(iCommas,1,splitRows, c(6,7))
    for (i in 1:length(newrows)){
      idx <- as.double(names(newrows)[i])
      len <- nrow(newrows[[i]])
      rownames(newrows[[i]]) <- seq(0,0.1*(len-1),length=len) + idx
      iTemp <- rbind(iTemp, newrows[[i]])
    }
  }
  if (nrow(dCommas) > 0){
    newrows <- apply(dCommas,1,splitRows, c(6,7))
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

rm(iTemp)
rm(dTemp)
rm(iCSV)
rm(dCSV)

ins.sep <- as.data.frame(rbindlist(ins.sep))
del.sep <- as.data.frame(rbindlist(del.sep))

ins.nosep <- as.data.frame(rbindlist(ins.nosep))
del.nosep <- as.data.frame(rbindlist(del.nosep))

ins <- ins.sep[ins.sep$indel!="",]
del <- del.sep[del.sep$indel!="",]


#------ Checkpoint: 10_2_finished.RData



# N - GLYC SITE OUTPUTS 
# ---------------------------------------------
write.table(ins,paste0(path, "13_nglycs/all/ins-sep.csv"), row.names=F, sep="\t", quote=F)
write.table(del,paste0(path, "13_nglycs/all/del-sep.csv"), row.names=F, sep="\t", quote=F)

# INDEL LENGTHS OUTPUT 
# ---------------------------------------------

write.csv(ins.sep[,c(1,2,5,6)], paste0(path,"12_lengths/all/ins-all.csv"))
write.csv(del.sep[,c(1,2,5,6)], paste0(path,"12_lengths/all/del-all.csv"))


# ---- Indel Nucleotide Analysis ----
# focus on only the columns needed 
total.ins <- ins[,c(2,3,6,9)]
total.del <- del[,c(2,3,6,9)]


total.ins$len <- sapply(total.ins$indel, nchar)
total.del$len <- sapply(total.del$indel, nchar)

total.ins$vlen <- as.numeric(total.ins$vlen)
total.del$vlen <- as.numeric(total.del$vlen)

total.ins$anc <- gsub("-","",total.ins$anc)
total.del$anc <- gsub("-","",total.del$anc)


write.csv(total.ins, paste0(path,"10_nucleotide/all/total-ins.csv"))
write.csv(total.del, paste0(path,"10_nucleotide/all/total-del.csv"))

# DINUCLEOTIDE PROPORTIONS OUTPUT 
# ------------------------------------

write.csv(total.ins[total.ins$len>1,], paste0(path,"10_nucleotide/all/dinucl-ins.csv"))
write.csv(total.del[total.del$len>1,], paste0(path,"10_nucleotide/all/dinucl-del.csv"))

# FLANKING INSERTIONS (14) OUTPUT 
# ------------------------------------
write.csv(ins, paste0(path,"/10_nucleotide/all/flanking/ins-sep.csv"))
write.csv(del, paste0(path,"/10_nucleotide/all/flanking/del-sep.csv"))

write.csv(ins.nosep, paste0(path,"/10_nucleotide/all/flanking/ins-nosep-all.csv"))
write.csv(del.nosep, paste0(path,"/10_nucleotide/all/flanking/del-nosep-all.csv"))

# --- Modeling 2 ----
write.csv(ins.sep, paste0(path,"/10_nucleotide/all/flanking/ins-sep-all.csv"))
write.csv(del.sep, paste0(path,"/10_nucleotide/all/flanking/del-sep-all.csv"))



# NT PROPORTIONS -- ALL
# --------------------------------------------

nucl <- c("A","C","G","T")
iProps <- c()
dProps <- c()
iVProps <- c()
dVProps <- c()
counts <- data.frame()
iTotals <- c(sum(unname(sapply(total.ins$indel, nchar))), sum(unname(sapply(total.ins$anc, nchar))))
dTotals <- c(sum(unname(sapply(total.del$indel, nchar))),sum(unname(sapply(total.del$anc, nchar))))

for (n in 1:4){
  icount <- sum(str_count(total.ins$indel, nucl[n]))
  dcount <- sum(str_count(total.del$indel, nucl[n]))
  counts <- rbind(counts, data.frame(nucl=nucl[n], ins=icount, del=dcount))
  
  iProps[n] <- icount / iTotals[1]
  dProps[n] <- dcount / dTotals[1]
  
  iVProps[n] <- sum(str_count(total.ins$anc, nucl[n])) / iTotals[2]
  dVProps[n] <- sum(str_count(total.del$anc, nucl[n])) / dTotals[2]
}
require(reshape)
counts2 <- reshape2::melt(counts)

ins.nt <- data.frame(nt=nucl,props=iProps,vprops=iVProps)
del.nt <- data.frame(nt=nucl,props=dProps,vprops=dVProps)
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


# generates the randomly sampled substrings for each indel
ires <- mapply(sampleString, total.ins[,"len"], total.ins[,"anc"])
dres <- mapply(sampleString, total.del[,"len"], total.del[,"anc"])


# compares the observed proportion to the overall distribution of each nucleotide 
isign <- c()
dsign <- c()
for (i in 1:4){
  idist <- unlist(ires[i,])
  ddist <- unlist(dres[i,])
  print(i)
  
  if(any(is.na(idist))){
    #print(i)
    idist <- idist[-which(is.na(idist))]
  }else if (any(is.na(ddist))){
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

# BOOTSTRAPS
# ----------------

bs.props <- list()

for (n in 1:1000){
  
  # create the bootstrap sample 
  sam1 <- sample(nrow(total.ins), nrow(total.ins), replace=T)
  sam2 <- sample(nrow(total.del), nrow(total.del), replace=T)
  
  df1.bs <- total.ins[sam1,]
  df2.bs <- total.del[sam2,]
  
  total1 <- c(sum(nchar(df1.bs[,'indel'])), sum(nchar(df1.bs[,'anc'])))
  total2 <- c(sum(nchar(df2.bs[,'indel'])), sum(nchar(df2.bs[,'anc'])))
  
  props <- sapply(nucl, function(nuc){
    iProps <-  sum(str_count(df1.bs$indel, nuc)) / total1[1]
    dProps <-  sum(str_count(df2.bs$indel, nuc)) / total2[1]
    
    iVProps <-  sum(str_count(df1.bs$anc, nuc)) / total1[2]
    dVProps <-  sum(str_count(df2.bs$anc, nuc)) / total2[2]
    
    return(c(iProps, dProps, iVProps, dVProps))
  })
 
  for (nuc in nucl){
    bs.props[[paste0("ins-",nuc)]] <- c(bs.props[[paste0("ins-",nuc)]], unname(props[1,nuc]))
    bs.props[[paste0("del-",nuc)]] <- c(bs.props[[paste0("del-",nuc)]], unname(props[2,nuc]))
    bs.props[[paste0("v-ins-",nuc)]] <- c(bs.props[[paste0("v-ins-",nuc)]], unname(props[3,nuc]))
    bs.props[[paste0("v-del-",nuc)]] <- c(bs.props[[paste0("v-del-",nuc)]], unname(props[4,nuc]))
  }
}
medians <- unlist(lapply(bs.props,median))
m.tab <- matrix(nrow=4,ncol=4, dimnames=list(nucl,c("ins","del","v-ins","v-del")))
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
colors <- c( "limegreen","dodgerblue","red", "purple")


par(pty="s", xpd=NA, mar=c(6,8,2,1),las=0)
clab = 1.9
ctext = 1.9
xpos <- c(0.18,0.162,0.245,0.38)
ypos <- c(0.147, 0.22, 0.26, 0.395)

lim = c(0.14,0.42)
plot(x=med.x, y=med.y, pch=indel.nt[,4]+1, col=rep(colors,2),xlim=lim,ylim=lim,cex=0.10*sqrt(indel.nt$counts),
     cex.lab=clab, cex.axis=1.6,lwd=6, ylab='', xlab='',las=1)#, main="Nucleotide Proportions")
title(ylab="Proportion In Indels", line=5,cex.lab=clab)
title(xlab="Proportion in Variable Loops", line=4,cex.lab=clab)

sqn <- c(seq(1,8,2),seq(1,8,2)+1)
# y error bars 
arrows(med.x, lower.y[sqn], med.x, upper.y[sqn], length=0, angle=90, code=3,lwd=1.5)
arrows(0.175,0.153,0.167,0.162, length=0, lwd=1.2)
# x error bars 
arrows(lower.x[sqn], med.y, upper.x[sqn], med.y, length=0, angle=90, code=3,lwd=1.5)
#legend(0.38,0.24,legend=nucleotides, pch=22,cex=1.3, pt.bg=indel.nt[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
text(xpos, ypos, cex=ctext, labels=c("C", "G", "T", "A"))
#legend(0.14,0.42,legend=c("Insertions", "Deletions"), pch=c(1,2),cex=1.3, lwd=2, col="black",x.intersp = 1.0,y.intersp=1.3, pt.cex=3)
par(xpd=F)
abline(0,1)
rect(0.145,0.36,0.206,0.415)
text(0.19, 0.40, labels="Ins", cex=ctext)
text(0.19, 0.376, labels="Del", cex=ctext)
points(c(0.16,0.16), c(0.40, 0.376), pch=c(1,2), cex=4, lwd=5, col='black', bg='black')






# NT PROPORTIONS BY VARIABLE LOOP 
# --------------------------------------------
vloops <- c(1, 2,3,4,5)
vloops2 <- c("V1","V2","V3","V4","V5")
ins.props <- data.frame()
del.props <- data.frame()

for (i in c(1,2,3,4,5)){
  iTemp <- total.ins[total.ins$vloop==i,]
  dTemp <- total.del[total.del$vloop==i,]

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


