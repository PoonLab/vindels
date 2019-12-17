require(bbmle)
require(stringr)
require(ape)

source("~/vindels/2_within-host/10_nt_utils.r")
# Lio
path <- "~/PycharmProjects/hiv-withinhost/"
 
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
  
  # retrieving the accno from the header
  iHeader <- unname(sapply(iCSV$header, getAccno))
  dHeader <- unname(sapply(dCSV$header, getAccno))
  
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

  ins.glycs <- rbind(ins.glycs, iCSV)
  del.glycs <- rbind(del.glycs, dCSV)
  
  # COMMA SEPARATION FIX
  # make a new data.frame for each CSV df
  # transport over all rows which do NOT contain a comma
  new.ins <- iCSV[!grepl(",",iCSV$Seq),]
  new.del <- dCSV[!grepl(",",dCSV$Seq),]
  
  # # handle comma rows separately with a function 
  # iCommas <- iCSV[grepl(",",iCSV$Seq),]
  # dCommas <- dCSV[grepl(",",dCSV$Seq),]
  # # APPLY THE SPLIT ROWS TO GET ONE INDEL PER ROW
  # if (nrow(iCommas) > 0){
  #   newrows <- apply(iCommas,1,splitRows)
  #   for (i in 1:length(newrows)){
  #     idx <- as.double(names(newrows)[i])
  #     len <- nrow(newrows[[i]])
  #     rownames(newrows[[i]]) <- seq(0,0.1*(len-1),length=len) + idx
  #     colnames(newrows[[i]]) <- c("Header", "Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq", "Anc", "Pat")
  #     new.ins <- rbind(new.ins, newrows[[i]])
  #   }
  # }
  # if (nrow(dCommas) > 0){
  #   newrows <- apply(dCommas,1,splitRows)
  #   for (i in 1:length(newrows)){
  #     idx <- as.double(names(newrows)[i])
  #     len <- nrow(newrows[[i]])
  #     rownames(newrows[[i]]) <- seq(0,0.1*len-0.1,length=len) + idx
  #     colnames(newrows[[i]]) <- c("Header", "Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq", "Anc","Pat")
  #     newnew.del <- rbind(new.del, newrows[[i]])
  #     
  #   }
  # }
  print("80% complete")
  # Retrieve variable loop positions from file 
  var.pos <- read.csv(paste0(path,"3RegionSequences/variable/", strsplit(filename, "-")[[1]][1], ".csv"), stringsAsFactors = F)
  var.pos <- var.pos[,-c(2,5,8,11,14)]
  
  new.ins[is.na(new.ins$Pos),"Pos"] <- ""
  new.del[is.na(new.del$Pos),"Pos"] <- ""
  
  # Add the V position column into the two final data frames 
  new.ins$Vpos <- as.numeric(unname(mapply(addPos, pos=new.ins$Pos, accno=new.ins$Header, vloop=new.ins$Vloop)))
  new.del$Vpos <- as.numeric(unname(mapply(addPos, pos=new.del$Pos, accno=new.del$Header, vloop=new.del$Vloop)))
  
  # ADJUST POSITIONS TO MATCH THE PLACE WHERE THE INSERTION WAS
  new.ins$Pos <- as.numeric(new.ins$Pos) - nchar(new.ins$Seq)
  new.ins$Vpos <- new.ins$Vpos - nchar(new.ins$Seq)
  
  # no adjustment needed for deletions
  new.del$Pos <- as.numeric(new.del$Pos)
  
  
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
write.table(ins.glycs2[,c(1,2,6,3,4,5,7)],paste0(path, "13_nglycs/ins.csv"), row.names=F, sep="\t", quote=F)
write.table(del.glycs2[,c(1,2,6,3,4,5,7)],paste0(path, "13_nglycs/del.csv"), row.names=F, sep="\t", quote=F)

# INDEL LENGTHS OUTPUT 
# ---------------------------------------------

write.csv(all.ins[,c(1,2,4,5,6)], "~/Lio/12_lengths/ins-full.csv")
write.csv(all.del[,c(1,2,4,5,6)], "~/Lio/12_lengths/del-full.csv")


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

total.ins <- total.ins[total.ins$len>1, ]
total.del <- total.del[total.del$len>1, ]


# DINUCLEOTIDE PROPORTIONS OUTPUT 
# ------------------------------------
write.csv(ins.glycs2[,c(1,2,6,3,4,5,7)], "~/Lio/10_nucleotide/total-ins.csv")
write.csv(del.glycs2[,c(1,2,6,3,4,5,7)], "~/Lio/10_nucleotide/total-del.csv")

# FLANKING INSERTIONS PROPORTIONS OUTPUT 
# ------------------------------------
write.csv(ins, "~/PycharmProjects/hiv-withinhost/10_nucleotide/ins-sep-only.csv")
write.csv(del, "~/PycharmProjects/hiv-withinhost/10_nucleotide/del-sep-only.csv")

write.csv(all.ins, "~/PycharmProjects/hiv-withinhost/10_nucleotide/ins-sep-all.csv")
write.csv(all.ins, "~/PycharmProjects/hiv-withinhost/10_nucleotide/del-sep-all.csv")

write.csv(ins.glycs, "~/PycharmProjects/hiv-withinhost/10_nucleotide/ins-nosep-all.csv")
write.csv(del.glycs, "~/PycharmProjects/hiv-withinhost/10_nucleotide/del-nosep-all.csv")



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
indel.nt$indel <- c(rep(1,4),rep(2,4))
indel.nt$counts <- counts$value

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
  dtemp <- sampleString(total.del[row,"len"], total.del[row,"Vseq"])
  for (i in 1:4){
    iSample[[i]] <- c(iSample[[i]], itemp[[i]])
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




# NT ALL INDELS PLOT 
# broken down by variable loop and nucleotide
# -------------------------------------
require(RColorBrewer)
colors <- brewer.pal(4, "Set1")

cex=1
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0.1,0.45)
plot(indel.nt[,c(3,2)], pch=indel.nt[,4]+21, bg=indel.nt[,1],xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.2,cex.main=1.8, ylab='', xlab='',cex=2.5, main="Nucleotide Proportions")
title(ylab="Proportion Inside Indels", line=3,cex.lab=1.3)
title(xlab="Proportion in Variable Loops", line=3,cex.lab=1.3)
legend(0.38,0.24,legend=nucleotides, pch=21,cex=1.3, pt.bg=indel.nt[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.10,0.45,legend=c("Insertions", "Deletions"), pch=c(22,23),cex=1.3, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
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



#A
cex=2
par(pty="s", mfrow=c(2,2), xpd=NA, mar=c(3,8,4,1),las=0)

lim = c(0.24,0.46)
plot(list.df[[1]][,1:2], cex=sizes.v, pch=(list.df[[1]]$vregion+20), bg=rep(colors, 7),xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=0.9,cex.main=1.5, ylab='', xlab='')
text(0.187,0.475,labels="a)", cex=1.5)
text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Within Indels", line=2.5,cex.lab=1.15)
title(xlab="Proportion Outside Indels", line=2.1,cex.lab=1.15)
par(xpd=F)
abline(0,1)
