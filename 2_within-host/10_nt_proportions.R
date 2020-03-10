require(bbmle)
require(stringr)
require(ape)

source("~/vindels/2_within-host/utils.r")
# Lio
path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
 
ifolder <- Sys.glob(paste0(path,"9Indels/mcc/ins/*.csv"))
dfolder <- Sys.glob(paste0(path,"9Indels/mcc/del/*.csv"))

# INSERTION PARSING ----------
#ifolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/ins_mcc/*.csv")
#dfolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/del_mcc/*.csv")
all.ins <- data.frame()
all.del <- data.frame()
iTotal <- list()
dTotal <- list()
count <- 0
ifull <- list()
dfull <- list()

ins.glycs <- data.frame(stringsAsFactors = F)
del.glycs <- data.frame(stringsAsFactors = F)

for (file in 1:length(ifolder)){
  print(file)
  filename <- basename(ifolder[file])
  
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F)
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F)
  
  # used for handling cases where there are no indels
  if (all(is.na(iCSV$ins))){
    iCSV$ins <- ""
  }
  if (all(is.na(dCSV$del))){
    dCSV$del <- ""
  }
  
  print("10% complete")
  # retrieving subtype field from the header
  iCSV$Subtype <- unname(sapply(iCSV$header, getSubtype))
  dCSV$Subtype <- unname(sapply(dCSV$header, getSubtype))
  
  
  # store the sequences from these two data frames for nucleotide analysis
  # remove them as they arent needed for this analysis
  ifull$var <- as.character(iCSV$Vseq)
  dfull$var <- as.character(dCSV$Vseq)
  ifull$anc <- as.character(dCSV$anc)
  dfull$anc <- as.character(dCSV$anc)
  
  dCSV$Vseq <- NULL
  iCSV$Vseq <- NULL
  dCSV$anc <- NULL
  iCSV$anc <- NULL
  
  # creates the counts column
  iCSV$Count <- sapply(iCSV$ins, csvcount) 
  dCSV$Count <- sapply(dCSV$del, csvcount)
  
  # extracts info from the indel column and puts it into two separate columns
  insInfo <- sapply(iCSV$ins, extractInfo)
  insInfo <- unname(insInfo)
  insInfo <- t(insInfo)
  insInfo <- as.data.frame(insInfo)
  insInfo$V1 <- as.character(insInfo$V1)
  insInfo$V2 <- as.character(insInfo$V2)
  iCSV <- cbind(iCSV, insInfo)
  iCSV$ins <- NULL
  
  delInfo <- sapply(dCSV$del, extractInfo)
  delInfo <- unname(delInfo)
  delInfo <- t(delInfo)
  delInfo <- as.data.frame(delInfo)
  delInfo$V1 <- as.character(delInfo$V1)
  delInfo$V2 <- as.character(delInfo$V2)
  dCSV <- cbind(dCSV, delInfo)
  dCSV$del <- NULL
  
  iCSV$Vseq <- ifull$var
  dCSV$Vseq <- dfull$var
  iCSV$Anc <- ifull$anc
  dCSV$Anc <- dfull$anc
  print("50% complete")
  iCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(iCSV))
  dCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(dCSV))
  
  colnames(iCSV) <- c("Header","Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq","Anc","Pat")
  colnames(dCSV) <- c("Header", "Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq","Anc","Pat")

  # ************************
  # Retrieve variable loop positions from file 
  treename <- strsplit(filename, "\\.csv")[[1]]
  
  # Load time-based branch lengths from the time-scaled trees
  tre <- read.nexus(paste0(path,"7_5_MCC/prelim/", treename, ".tree"))
  branches <- tre$edge.length[tre$edge[,2] <=Ntip(tre)]  
  
  iCSV$Date <- branches[match(gsub("_\\d$","",iCSV$Header), tre$tip.label)]
  dCSV$Date <- branches[match(gsub("_\\d$","",dCSV$Header), tre$tip.label)]
  
  # Load subs/site branch lengths from the rescaled trees
  tre <- read.tree(paste0(path,"7_5_MCC/final/",treename , ".tree"))
  branches <- tre$edge.length[tre$edge[,2] <=Ntip(tre)]  
  
  iCSV$length <- branches[match(gsub("_\\d$","",iCSV$Header), tre$tip.label)]
  dCSV$length <- branches[match(gsub("_\\d$","",dCSV$Header), tre$tip.label)]
  
  iCSV$Header <- unname(mapply(labels, iCSV$Header, iCSV$Pat, iCSV$Vloop))
  dCSV$Header <- unname(mapply(labels, dCSV$Header, dCSV$Pat, dCSV$Vloop))
  
  ins.glycs <- rbind(ins.glycs, iCSV)
  del.glycs <- rbind(del.glycs, dCSV)
  # COMMA SEPARATION FIX
  # make a new data.frame for each CSV df
  # transport over all rows which do NOT contain a comma
  new.ins <- iCSV[!grepl(",",iCSV$Seq),]
  new.del <- dCSV[!grepl(",",dCSV$Seq),]
  
  # handle comma rows separately with a function
  iCommas <- iCSV[grepl(",",iCSV$Seq),]
  dCommas <- dCSV[grepl(",",dCSV$Seq),]
  # APPLY THE SPLIT ROWS TO GET ONE INDEL PER ROW
  if (nrow(iCommas) > 0){
    newrows <- apply(iCommas,1,splitRows, colnum=12)
    for (i in 1:length(newrows)){
      idx <- as.double(names(newrows)[i])
      len <- nrow(newrows[[i]])
      rownames(newrows[[i]]) <- seq(0,0.1*(len-1),length=len) + idx
      colnames(newrows[[i]]) <- c("Header", "Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq", "Anc", "Pat", "Date","length")
      new.ins <- rbind(new.ins, newrows[[i]])
    }
  }
  if (nrow(dCommas) > 0){
    newrows <- apply(dCommas,1,splitRows, colnum=12)
    for (i in 1:length(newrows)){
      idx <- as.double(names(newrows)[i])
      len <- nrow(newrows[[i]])
      rownames(newrows[[i]]) <- seq(0,0.1*len-0.1,length=len) + idx
      colnames(newrows[[i]]) <- c("Header", "Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq", "Anc","Pat", "Date","length")
      newnew.del <- rbind(new.del, newrows[[i]])

    }
  }
  print("80% complete")
  
  new.ins[is.na(new.ins$Pos),"Pos"] <- ""
  new.del[is.na(new.del$Pos),"Pos"] <- ""
  
  # Add the V position column into the two final data frames 
  new.ins$Vpos <- as.numeric(unname(mapply(addPos, pos=new.ins$Pos, header=new.ins$Header, vloop=new.ins$Vloop)))
  new.del$Vpos <- as.numeric(unname(mapply(addPos, pos=new.del$Pos, header=new.del$Header, vloop=new.del$Vloop)))
  
  # # ADJUST POSITIONS TO MATCH THE PLACE WHERE THE INSERTION WAS
  # new.ins$Pos <- as.numeric(new.ins$Pos) - nchar(new.ins$Seq)
  # new.ins$Vpos <- new.ins$Vpos - nchar(new.ins$Seq)
  
  # no adjustment needed for deletions
  # new.del$Pos <- as.numeric(new.del$Pos)
  
  
  # OUTPUT 
  # for other analyses
  # -----------------------------
  all.ins <- rbind(all.ins, new.ins)
  all.del <- rbind(all.del, new.del)

}
ins <- all.ins[all.ins$Seq!="",]
del <- all.del[all.del$Seq!="",]

# N - GLYC SITE OUTPUTS 
# ---------------------------------------------
ins.glycs2 <- ins.glycs[ins.glycs$Seq!="",-c(3,4,5)]
del.glycs2 <- del.glycs[del.glycs$Seq!="",-c(3,4,5)]
write.table(ins,paste0(path, "13_nglycs/ins-sep.csv"), row.names=F, sep="\t", quote=F)
write.table(del,paste0(path, "13_nglycs/del-sep.csv"), row.names=F, sep="\t", quote=F)

# INDEL LENGTHS OUTPUT 
# ---------------------------------------------

write.csv(all.ins[,c(1,2,4,5,6)], paste0(path,"12_lengths/ins-full.csv"))
write.csv(all.del[,c(1,2,4,5,6)], paste0(path,"12_lengths/del-full.csv"))



nucleotides <- c("A","C","G","T")
ntcount <- c()
total.ins <- data.frame()
total.del <- data.frame()
# Read through every CSV file in the ins and del folders 
for (n in c(1,2,4,5)){
  ins.df <- all.ins[all.ins$Vloop==n,c(1,2,6,8)]
  del.df <- all.del[all.del$Vloop==n,c(1,2,6,8)]
  # if the csv is not entirely blank with NAs 
  if (all(!is.na(ins.df$Seq))){
    ins.df$Seq <- sapply(ins.df$Seq, removeNA)
    del.df$Seq <- sapply(del.df$Seq, removeNA)
    #colnames(ins.df) <- c("Accno", "Vloop", "Seq", "VSeq", "Run")
    #colnames(del.df) <- c("Accno", "Vloop", "Seq", "VSeq", "Run")
    
    total.ins <- rbind(total.ins, ins.df[ins.df$Seq!="",])
    total.del <- rbind(total.del, del.df[del.df$Seq!="",])
  }
}
total.ins$len <- sapply(total.ins$Seq, nchar)
total.del$len <- sapply(total.del$Seq, nchar)
write.csv(total.ins, paste0(path,"10_nucleotide/total-ins.csv"))
write.csv(total.del, paste0(path,"10_nucleotide/total-del.csv"))




# DINUCLEOTIDE PROPORTIONS OUTPUT 
# ------------------------------------
total.ins2 <- total.ins[total.ins$len>1, ]
total.del2 <- total.del[total.del$len>1, ]
write.csv(total.ins2, paste0(path,"10_nucleotide/total-ins2.csv"))
write.csv(total.del2, paste0(path,"10_nucleotide/total-del2.csv"))

# FLANKING INSERTIONS PROPORTIONS OUTPUT 
# ------------------------------------
write.csv(ins, paste0(path,"/10_nucleotide/ins-sep-only.csv"))
write.csv(del, paste0(path,"/10_nucleotide/del-sep-only.csv"))

write.csv(all.ins, paste0(path,"/10_nucleotide/ins-sep-all.csv"))
write.csv(all.ins, paste0(path,"/10_nucleotide/del-sep-all.csv"))

write.csv(ins.glycs, paste0(path,"/10_nucleotide/ins-nosep-all.csv"))
write.csv(del.glycs, paste0(path,"/10_nucleotide/del-nosep-all.csv"))



# NT PROPORTIONS -- ALL
# ---------------------------------------------
iProps <- c()
dProps <- c()
iVProps <- c()
dVProps <- c()
counts <- data.frame()
iTotals <- c(sum(unname(sapply(total.ins$Seq, nchar))), sum(unname(sapply(total.ins$Vseq, nchar))))
dTotals <- c(sum(unname(sapply(total.del$Seq, nchar))),sum(unname(sapply(total.del$Vseq, nchar))))

for (nuc in nucleotides){
  icount <- sum(str_count(total.ins$Seq, nuc))
  dcount <- sum(str_count(total.del$Seq, nuc))
  counts <- rbind(counts, data.frame(nucl=nuc, ins=icount, del=dcount))
  
  iProps <- c(iProps, icount / iTotals[1])
  dProps <- c(dProps, dcount / dTotals[1])
  
  iVProps <- c(iVProps, sum(str_count(total.ins$Vseq, nuc)) / iTotals[2])
  dVProps <- c(dVProps, sum(str_count(total.del$Vseq, nuc)) / dTotals[2])
}
require(reshape)
counts <- melt(counts)

ins.nt <- data.frame(nt=nucleotides,props=iProps,vprops=iVProps)
del.nt <- data.frame(nt=nucleotides,props=dProps,vprops=dVProps)
indel.nt <- rbind(ins.nt, del.nt)
indel.nt$indel <- c(rep(0,4),rep(3,4))
indel.nt$counts <- counts$value

# BOOTSTRAPS
# ----------------

df1 <- total.ins
df2 <- total.del
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
    iProps <-  sum(str_count(df1.bs$Seq, nuc)) / total1[1]
    dProps <-  sum(str_count(df2.bs$Seq, nuc)) / total2[1]
    
    iVProps <-  sum(str_count(df1$Vseq, nuc)) / total1[2]
    dVProps <-  sum(str_count(df2$Vseq, nuc)) / total2[2]
    
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
colors <- brewer.pal(4, "Set1")

cex=1
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0.1,0.45)
plot(indel.nt[,c(3,2)], pch=indel.nt[,4]+21, bg=indel.nt[,1],xlim=lim,ylim=lim,cex=0.12*sqrt(indel.nt$counts),
     cex.lab=1.3, cex.axis=1.2,cex.main=1.8, ylab='', xlab='', main="Nucleotide Proportions")
title(ylab="Proportion Inside Indels", line=3,cex.lab=1.3)
title(xlab="Proportion in Variable Loops", line=3,cex.lab=1.3)
arrows(indel.nt[1:8,3], lower.y[c(seq(1,8,2),seq(1,8,2)+1)], indel.nt[1:8,3], upper.y[c(seq(1,8,2),seq(1,8,2)+1)], length=0.05, angle=90, code=3)
arrows(lower.x[c(seq(1,8,2),seq(1,8,2)+1)], indel.nt[1:8,2], upper.x[c(seq(1,8,2),seq(1,8,2)+1)], indel.nt[1:8,2], length=0.05, angle=90, code=3)
legend(0.38,0.24,legend=nucleotides, pch=22,cex=1.3, pt.bg=indel.nt[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.10,0.45,legend=c("Insertions", "Deletions"), pch=c(21,24),cex=1.3, pt.bg="black",x.intersp = 1.0,y.intersp=1.3, pt.cex=3)
par(xpd=F)
abline(0,1)






# NT PROPORTIONS BY VARIABLE LOOP 
# --------------------------------------------
vloops <- c(1, 2,3,4,5)
vloops2 <- c("V1","V2","V4","V5")
ins.props <- data.frame()
del.props <- data.frame()

for (i in c(1,2,4,5)){
  iTemp <- total.ins[total.ins$Vloop==i,]
  dTemp <- total.del[total.del$Vloop==i,]

  iProps <- c()
  dProps <- c()
  
  iVProps <- c()
  dVProps <- c()
  
  # a vector of two totals
  # iTotals[1] = total number of nucleotides in insertion sequences
  # iTotals[2] = total nuimebr of nucleotides in the vloops associated with insertion sequences
  iTotals <- c(sum(unname(sapply(iTemp$Seq, nchar))), sum(unname(sapply(iTemp$Vseq, nchar))))
  dTotals <- c(sum(unname(sapply(dTemp$Seq, nchar))),sum(unname(sapply(dTemp$Vseq, nchar))))
  
  for (nuc in nucleotides){
    iProps <- c(iProps, sum(str_count(iTemp$Seq, nuc)) / iTotals[1])
    dProps <- c(dProps, sum(str_count(dTemp$Seq, nuc)) / dTotals[1])
    
    iVProps <- c(iVProps, sum(str_count(iTemp$Vseq, nuc)) / iTotals[2])
    dVProps <- c(dVProps, sum(str_count(dTemp$Vseq, nuc)) / dTotals[2])
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



iSample <- list(c(),c(),c(),c())
dSample <- list(c(),c(),c(),c())
# generates the randomly sampled substrings for each indel
for (row in 1:nrow(total.ins)){
  itemp <- sampleString(total.ins[row,"len"], total.ins[row,"Vseq"])
  for (i in 1:4){
    iSample[[i]] <- c(iSample[[i]], itemp[[i]])
  }
}
for (row in 1:nrow(total.del)){
  dtemp <- sampleString(total.del[row,"len"], total.del[row,"Vseq"])
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

