require(bbmle)
require(stringr)
require(ape)

csvcount <- function(input){
  commas <- str_count(input, ",")
  if (commas > 0){
    result <- commas + 1  
  }else if(input == ""){
    result <- 0
  }else{
    result <- 1
  }
  result
}

# used for extracting condensed CSV information 
extractInfo <- function(input){
  if (length(input)==1 && input == ""){
    return(c("",""))
  }else{
    insertions <- strsplit(input, ":")
  }
  seq <- c()
  pos <- c()
  
  for (ins in insertions[[1]]){
    fields <- strsplit(ins, "-")
    seq <- c(seq, fields[[1]][1])
    pos <- c(pos, as.numeric(fields[[1]][2]))
  }
  return(c(paste(seq,collapse=","), paste(pos,collapse=",")))
}

# for changing headers from ACCNO_DATE format to ACCNO
getSubtype <- function(header){
  newheader <- str_split(as.character(header),"\\.")[[1]][1]
  newheader
}

# used for handling entire columns of NA values
removeNA <- function(input){
  if (all(is.na(input))){
    input <- ""
  }
  input
}

# retrieves the accno field from the full scale header 
getAccno <- function(input){
  accno <- strsplit(input, "\\.")[[1]][5]
  accno
}


# specifically handles fields containing a comma
splitRows <- function(row){
  row <- data.frame(t(row),stringsAsFactors = F)
  seqs <- str_split(row[1,6], ",")[[1]]
  pos <- str_split(row[1,7],",")[[1]]
  len <- length(seqs)
  #print(seqs)
  data.frame(row[rep(1,len),1:5], Seq=seqs, Pos=pos, row[rep(1,len),8])
  
}



# INSERTION PARSING ----------
ifolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/ins_mcc/*.csv")
dfolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/del_mcc/*.csv")
all.ins <- data.frame()
all.del <- data.frame()
iTotal <- list()
dTotal <- list()
count <- 0
sequences <- list()


for (file in 1:length(ifolder)){
  print(file)
  filename <- basename(ifolder[file])
  
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F)
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F)
  
  # used for handling cases where there are no indels
  if (all(is.na(iCSV$Ins))){
    iCSV$Ins <- ""
  }
  if (all(is.na(dCSV$Del))){
    dCSV$Del <- ""
  }
  
  # retrieving subtype field from the header
  iCSV$Subtype <- unname(sapply(iCSV$Accno, getSubtype))
  dCSV$Subtype <- unname(sapply(dCSV$Accno, getSubtype))
  
  # retrieving the accno from the header
  iAccno <- unname(sapply(iCSV$Accno, getAccno))
  dAccno <- unname(sapply(dCSV$Accno, getAccno))
  
  # 
  iCSV$Accno <- iAccno
  dCSV$Accno <- dAccno
  
  # store the sequences from these two data frames for nucleotide analysis
  sequences$ins <- as.character(iCSV$Seq)
  sequences$del <- as.character(dCSV$Seq)
  
  # remove them as they arent needed for this analysis
  dCSV$Seq <- NULL
  iCSV$Seq <- NULL
  
  # creates the counts column
  iCSV$Count <- sapply(iCSV$Ins, csvcount) 
  dCSV$Count <- sapply(dCSV$Del, csvcount)
  
  # extracts info from the indel column and puts it into two separate columns
  insInfo <- sapply(iCSV$Ins, extractInfo)
  insInfo <- unname(insInfo)
  insInfo <- t(insInfo)
  insInfo <- as.data.frame(insInfo)
  insInfo$V1 <- as.character(insInfo$V1)
  insInfo$V2 <- as.character(insInfo$V2)
  iCSV <- cbind(iCSV, insInfo)
  iCSV$Ins <- NULL
  
  delInfo <- sapply(dCSV$Del, extractInfo)
  delInfo <- unname(delInfo)
  delInfo <- t(delInfo)
  delInfo <- as.data.frame(delInfo)
  delInfo$V1 <- as.character(delInfo$V1)
  delInfo$V2 <- as.character(delInfo$V2)
  dCSV <- cbind(dCSV, delInfo)
  dCSV$Del <- NULL
  
  iCSV$Vseq <- sequences$ins
  dCSV$Vseq <- sequences$del
  
  colnames(iCSV) <- c("Accno","Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq")
  colnames(dCSV) <- c("Accno", "Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq")
  
  # COMMA SEPARATION FIX
  
  new.ins <- data.frame()
  new.del <- data.frame()
  # make a new data.frame for each CSV df
  # transport over all rows which do NOT contain a comma
  new.ins <- iCSV[!grepl(",",iCSV$Seq),]
  new.del <- dCSV[!grepl(",",dCSV$Seq),]
  
  # handle comma rows separately with a function 
  iCommas <- iCSV[grepl(",",iCSV$Seq),]
  dCommas <- dCSV[grepl(",",dCSV$Seq),]
  #c()
  if (nrow(iCommas) > 0){
    newrows <- apply(iCommas,1,splitRows)
    for (i in 1:length(newrows)){
      idx <- as.double(names(newrows)[i])
      len <- nrow(newrows[[i]])
      rownames(newrows[[i]]) <- seq(0,0.1*len-0.1,length=len) + idx
      colnames(newrows[[i]]) <- c("Accno", "Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq")
      new.ins <- rbind(new.ins, newrows[[i]])
    }
  }
  if (nrow(dCommas) > 0){
    newrows <- apply(dCommas,1,splitRows)
    for (i in 1:length(newrows)){
      idx <- as.double(names(newrows)[i])
      len <- nrow(newrows[[i]])
      rownames(newrows[[i]]) <- seq(0,0.1*len-0.1,length=len) + idx
      colnames(newrows[[i]]) <- c("Accno", "Vloop", "Vlength","Subtype", "Count", "Seq", "Pos", "Vseq")
      new.del <- rbind(new.del, newrows[[i]])
    }
  }
  
  # OUTPUT 
  # for other analyses
  # -----------------------------
  
  all.ins <- rbind(all.ins, new.ins)
  all.del <- rbind(all.del, new.del)

}


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
nucleotides <- c("A","C","G","T")


# NT PROPORTIONS -- ALL
# ---------------------------------------------
iProps <- c()
dProps <- c()
iVProps <- c()
dVProps <- c()
iTotals <- c(sum(unname(sapply(total.ins$Seq, nchar))), sum(unname(sapply(total.ins$Vseq, nchar))))
dTotals <- c(sum(unname(sapply(total.del$Seq, nchar))),sum(unname(sapply(total.del$Vseq, nchar))))

for (nuc in nucleotides){
  iProps <- c(iProps, sum(str_count(total.ins$Seq, nuc)) / iTotals[1])
  dProps <- c(dProps, sum(str_count(total.del$Seq, nuc)) / dTotals[1])
  
  iVProps <- c(iVProps, sum(str_count(total.ins$Vseq, nuc)) / iTotals[2])
  dVProps <- c(dVProps, sum(str_count(total.del$Vseq, nuc)) / dTotals[2])
}

indel.nt <- data.frame(nt=rep(nucleotides,2),indel=c(rep(1,4),rep(2,4)),props=c(iProps, dProps),vprops=c(iVProps,dVProps))



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


total.ins$len <- sapply(total.ins$Seq, nchar)
total.del$len <- sapply(total.del$Seq, nchar)

#total.ins$indel <- rep(TRUE, nrow(total.ins))
#total.del$indel <- rep(FALSE, nrow(total.del))

#all.ins$len <- sapply(all.ins$Seq, getLength)
#all.del$len <- sapply(all.del$Seq, getLength)
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
  
  ins.p <- indel.nt[i,3]
  del.p <- indel.nt[i+4,3]
  
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


ins.final <- cbind(total.ins[,c(1,2,3,5)], isign)
del.final <- cbind(total.del[,c(1,2,3,5)], dsign)




dinucleotide <- function(seq){
  ditotal <- nchar(seq)-1
  dinucl <- matrix(nrow=4,ncol=4)
  colnames(dinucl) <- c("A","C","G","T")
  rownames(dinucl) <- c("A","C","G","T")
  for (n in 1:(nchar(seq)-1)){
    di <- substr(seq,n,n+1)
    pos1 <- substr(di,1,1)
    pos2 <- substr(di,2,2)
    
    if (is.na(dinucl[pos2,pos1])){
      dinucl[pos2,pos1] <- 1
    }else{
      dinucl[pos2,pos1] <- dinucl[pos2,pos1] + 1
    }
  }
  dinucl
}

# DINUCLEOTIDE PROPORTIONS 

for (row in 1:nrow(total.ins)){
  
  
}




# NT ALL INDELS PLOT 
# broken down by variable loop and nucleotide
# -------------------------------------
require(RColorBrewer)
colors <- brewer.pal(4, "Set1")

cex=2
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0.1,0.45)
plot(indel.nt[,c(4,3)], pch=indel.nt[,2]+21, bg=indel.nt[,1],xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Nucleotide Proportions")
title(ylab="Proportion Inside Indels", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.38,0.24,legend=nucleotides, pch=21,cex=1.9, pt.bg=ins.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.10,0.45,legend=c("Insertions", "Deletions"), pch=c(22,23),cex=1.9, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
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
legend(0.43,0.18,legend=nucleotides, pch=21,cex=1.9, pt.bg=ins.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.33,0.18,legend=vloops2, pch=c(21,22,24,25),cex=1.9, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
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

lim = c(0.1,0.45)
plot(del.props[,c(3,2)], pch=del.props[,4]+20, bg=del.props[,1],xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Deletions")
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Inside Deletions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.43,0.18,legend=nucleotides, pch=21,cex=1.9, pt.bg=del.props[,1],x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
legend(0.33,0.18,legend=vloops2, pch=c(21,22,24,25),cex=1.9, pt.bg="black",x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
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
