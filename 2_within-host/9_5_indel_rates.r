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
  if (is.na(input)){
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
  seqs <- str_split(row[1,7], ",")[[1]]
  pos <- str_split(row[1,8],",")[[1]]
  len <- length(seqs)
  #print(seqs)
  data.frame(row[rep(1,len),1:6], Seq=seqs, Pos=pos, row[rep(1,len),9:10])

}


# INSERTION PARSING ----------
ifolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/rep/ins/*.csv")
dfolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/rep/ins/*.csv")
all.ins <- data.frame()
all.del <- data.frame()
csv.ins <- data.frame()
csv.del <- data.frame()

iTotal <- list()
dTotal <- list()
count <- 0
sequences <- list()


for (file in 1:length(ifolder)){
  print(file)
  filename <- strsplit(basename(ifolder[file]),"\\.")[[1]][1]
  runno <- strsplit(filename, "_")[[1]][2]
  count <- count + 1
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F)
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F)
  
  # used for handling cases where there are no indels
  if (all(is.na(iCSV$ins))){
    iCSV$ins <- ""
  }
  if (all(is.na(dCSV$del))){
    dCSV$del <- ""
  }
  
  # retrieving subtype field from the header
  iCSV$Subtype <- unname(sapply(iCSV$header, getSubtype))
  dCSV$Subtype <- unname(sapply(dCSV$header, getSubtype))
  
  # retrieving the accno from the header
  #iAccno <- unname(sapply(iCSV$Accno, getAccno))
  #dAccno <- unname(sapply(dCSV$Accno, getAccno))
  
  # 
  #iCSV$Accno <- iAccno
  #dCSV$Accno <- dAccno

  
  # creates the counts column
  iCSV$Count <- sapply(iCSV$ins, csvcount) 
  dCSV$Count <- sapply(dCSV$del, csvcount)
  
  # reads in the tree
  tre <- read.tree(paste0("~/PycharmProjects/hiv-withinhost/7SampleTrees/prelim/",filename , ".tree.sample"))
  
  # adjusts the tre tip labels to match the accession numbers
  #tre$tip.label <- unname(sapply(tre$tip.label, function(x){strsplit(x,"_")[[1]][1]}))
  
  # retrieves branch lengths from the tree
  branches <- tre$edge.length[tre$edge[,2] <=Ntip(tre)]   
  
  # matches the branch length to each of the sequences 
  iCSV$Date <- branches[match(sub("_\\d*$","",iCSV$header), tre$tip.label)]
  dCSV$Date <- branches[match(sub("_\\d*$","",dCSV$header), tre$tip.label)]
  
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
  
  iCSV$Run <- rep(runno, nrow(iCSV))
  dCSV$Run <- rep(runno, nrow(dCSV))
  
  iCSV <- iCSV[c(1,2,3,6,7,8,9,10,4,5,11)]
  dCSV <- dCSV[c(1,2,3,6,7,8,9,10,4,5,11)]
  
  colnames(iCSV) <- c("header","Vloop", "Vlength","Subtype", "Count","Date", "Seq", "Pos", "Vseq", "Anc", "Run")
  colnames(dCSV) <- c("header", "Vloop", "Vlength","Subtype", "Count","Date", "Seq", "Pos", "Vseq", "Anc", "Run")
  
  # COMMA SEPARATION FIX
  
  # new.ins <- data.frame()
  # new.del <- data.frame()
  # # make a new data.frame for each CSV df
  # # transport over all rows which do NOT contain a comma
  # new.ins <- iCSV[!grepl(",",iCSV$Seq),]
  # new.del <- dCSV[!grepl(",",dCSV$Seq),]
  # 
  # # handle comma rows separately with a function 
  # iCommas <- iCSV[grepl(",",iCSV$Seq),]
  # dCommas <- dCSV[grepl(",",dCSV$Seq),]
  # #c()
  # if (nrow(iCommas) > 0){
  #   newrows <- apply(iCommas,1,splitRows)
  #   for (i in 1:length(newrows)){
  #     idx <- as.double(names(newrows)[i])
  #     len <- nrow(newrows[[i]])
  #     rownames(newrows[[i]]) <- seq(0,0.1*len-0.1,length=len) + idx
  #     new.ins <- rbind(new.ins, newrows[[i]])
  #   }
  #   #new.ins <- new.ins[order(as.double(rownames(new.ins)))]
  # }
  # if (nrow(dCommas) > 0){
  #   newrows <- apply(dCommas,1,splitRows)
  #   for (i in 1:length(newrows)){
  #     idx <- as.double(names(newrows)[i])
  #     len <- nrow(newrows[[i]])
  #     rownames(newrows[[i]]) <- seq(0,0.1*len-0.1,length=len) + idx
  #     new.del <- rbind(new.del, newrows[[i]])
  #   }
  #   #new.del <- new.del[order(as.double(rownames(new.del)))]
  # }
  # 
  # # OUTPUT 
  # # for other analyses
  # # -----------------------------
  # 
  # all.ins <- rbind(all.ins, new.ins)
  # all.del <- rbind(all.del, new.del)
  # 
  
  # OUTPUT 2 
  # used for indel rates 
  
  csv.ins <- rbind(csv.ins, iCSV)
  csv.del <- rbind(csv.del, dCSV)

  # if (!is.null(iTotal[[runno]])){
  #   iTotal[[runno]] <- rbind(iTotal[[runno]], iCSV)
  # }else{
  #   iTotal[[runno]] <- iCSV
  # }
  # 
  # if (!is.null(dTotal[[runno]])){
  #   dTotal[[runno]] <- rbind(dTotal[[runno]], dCSV)
  # }else{
  #   dTotal[[runno]] <- dCSV
  # }
}
iTotal <- split(csv.ins, csv.ins$Run)
dTotal <- split(csv.del, csv.del$Run)



require(BSDA)
# RATE ANALYSIS -------------
ins.df <- data.frame()
del.df <- data.frame()

#median(as.numeric(all.ins[all.ins$Vloop==2,3])) # used to determine the median lengths of the variable loops
vlengths <- c(72,126,105,87,33)
all.df <- data.frame()
for (run in 1:20){
  iData <- iTotal[[run]]
  dData <- dTotal[[run]]
  
  irates <- c()
  drates <- c()
  
  for (vloop in 1:5){
    itemp <- iData[iData$Vloop==vloop,]
    dtemp <- dData[dData$Vloop==vloop,]
    
    iFinal <- itemp[itemp$Date < 325 & itemp$Count < 2,]
    dFinal <- dtemp[dtemp$Date < 325 & dtemp$Count < 2,]
    #print(nrow(current) - nrow(iFinal))
    
    ifit <- glm(iFinal$Count ~ 1, offset=log(iFinal$Date), family="poisson")
    irate <- exp(coef(ifit)[[1]])*365/vlengths[vloop]
    irates <- c(irates, irate)
    print(summary(ifit))
    
    dfit <- glm(dFinal$Count ~ 1, offset=log(dFinal$Date), family="poisson")
    drate <- exp(coef(dfit)[[1]])*365/vlengths[vloop]
    drates <- c(drates, drate)
    print(summary(dfit))
    
    #print(1 - (fit$deviance/fit$null.deviance))
    #all.df <- rbind(all.df, iFinal)
    #EDA(residuals(fit))
    
    #to.remove <- c(order(residuals(fit), decreasing = T)[1:52])
    #iFinal2 <- iFinal[-to.remove,]
    #fit2 <- glm(iFinal2$Count ~ iFinal2$Date, family="poisson")
    #EDA(residuals(fit2))
  }
  
  irates <- irates*10^3
  drates <- drates*10^3
  
  # contain the 20 
  ins.df <- rbind(ins.df, data.frame(V1=irates[1],V2=irates[2],V3=irates[3],V4=irates[4],V5=irates[5]))
  del.df <- rbind(del.df, data.frame(V1=drates[1],V2=drates[2],V3=drates[3],V4=drates[4],V5=drates[5]))
}

vloops <- c("V1","V2","V3","V4","V5")


# BOXPLOTS 
# --------------------------
par(mar=c(6,6,3,2))
boxplot(ins.df, main="Insertions", xlab="Variable Loop", ylab="Events/Nt/Year x 10^-3",cex.main=1.8,cex.lab=1.6,cex.axis=1.4)
boxplot(del.df, main="Deletions", xlab="Variable Loop", ylab="Events/Nt/Year x 10^-3",cex.main=1.8,cex.lab=1.6,cex.axis=1.4)
ins.df <- t(ins.df)
del.df <- t(del.df)



#insrates <- data.frame(VLoop=vloops, iRate=irates, AdjRate=irates*10^3)
#delrates <- data.frame(VLoop=vloops, dRate=drates, AdjRate=drates*10^3)

insrates <- data.frame(vloop=vloops,
                       rate=apply(ins.df, 1, median), 
                       lower=apply(ins.df,1,function(x){quantile(x, c(0.025,0.975))[1]}), 
                       upper=apply(ins.df,1,function(x){quantile(x, c(0.025,0.975))[2]}))
delrates <- data.frame(vloop=vloops,
                       rate=apply(del.df, 1, median), 
                       lower=apply(del.df,1,function(x){quantile(x, c(0.025,0.975))[1]}), 
                       upper=apply(del.df,1,function(x){quantile(x, c(0.025,0.975))[2]}))

#indels <- cbind(insrates, delrates[,c(2,3)])

# INDEL ABOVE/BELOW MULTIPLOTS 
# -------------------------------------

require(Rmisc)
require(ggplot2)

g1 <- ggplot(insrates, aes(x=vloop, y=rate,width=0.8)) + 
  geom_bar(colour="black", stat="identity",fill="dodgerblue",position="dodge",show.legend=F) +
  geom_errorbar(aes(ymax = insrates$upper, ymin = insrates$lower), 
                width = 0.25, size=1.1) +
  labs(x="Variable Loop", 
       y=expression(paste("       Insertion Rate \n(Events/Nt/Year x  ", 10^-3, ")", sep = "")))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 5)) +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 42, r = 10, b = 8, l = 18, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=18,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        axis.title.x=element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        legend.position="none")+ geom_text(aes(y=0.3,x=3 ),
                                           label="N/A", 
                                           size=6)
g1


g2 <- ggplot(delrates, aes(x=vloop, y=rate,width=0.8)) + 
  geom_bar(colour="black", stat="identity",fill="firebrick1",position="dodge",show.legend=F) + 
  geom_errorbar(aes(ymax = delrates$upper, ymin = delrates$lower), 
                width = 0.25, size=1.1) +
  labs(x="Variable Loop", 
       y="Deletion Rate")+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_y_reverse(lim=c(5,0))+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = -10, r = 10, b = 8, l = 24, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=18,margin=margin(t = 0, r = 11, b = 0, l = 6)),
        axis.title.x=element_text(size=18,margin=margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size=16),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=14),
        legend.position="none") + geom_text(aes(y=0.3,x=3 ),
                                            label="N/A", 
                                            size=6)
  
multiplot(g1,g2)



# SUBTYPE STRATIFICATION 
# ----------------------------------------------
ins <- list()
ins$B <- csv.ins[csv.ins$Subtype=="B",]
ins$C <- csv.ins[csv.ins$Subtype=="C",]
del <- list()
del$B <- csv.del[csv.del$Subtype=="B",]
del$C <- csv.del[csv.del$Subtype=="C",]

sub_irates <- data.frame()
sub_drates <- data.frame()

for (i in 1:2){
  
  iData <- ins[[i]]
  dData <- del[[i]]
  irates <- c()
  drates <- c()
  
  st <- c("B", "C")
  
  for (vloop in 1:5){
    iFinal <- iData[iData$Vloop==vloop & iData$Date < 325 & iData$Count < 2,]
    dFinal <- dData[dData$Vloop==vloop & dData$Date < 325 & dData$Count < 2,]
    #print(nrow(current) - nrow(iFinal))
    
    ifit <- glm(iFinal$Count ~ 1, offset=log(iFinal$Date), family="poisson")
    irate <- exp(coef(ifit)[[1]])*365/vlengths[vloop]
    irates <- c(irates, irate)
    print(summary(ifit))
    
    dfit <- glm(dFinal$Count ~ 1, offset=log(dFinal$Date), family="poisson")
    drate <- exp(coef(dfit)[[1]])*365/vlengths[vloop]
    drates <- c(drates, drate)
    print(summary(dfit))
    
    #print(1 - (fit$deviance/fit$null.deviance))
    #all.df <- rbind(all.df, iFinal)
    #EDA(residuals(fit))
    
    #to.remove <- c(order(residuals(fit), decreasing = T)[1:52])
    #iFinal2 <- iFinal[-to.remove,]
    #fit2 <- glm(iFinal2$Count ~ iFinal2$Date, family="poisson")
    #EDA(residuals(fit2))
  }
  sub_irates <- rbind(sub_irates, data.frame(rate=irates, adjrate=(irates*10^3),subtype=rep(st[[i]],5), VLoop=vloops))
  sub_drates <- rbind(sub_drates, data.frame(rate=drates, adjrate=(drates*10^3),subtype=rep(st[[i]],5), VLoop=vloops))
}


# SUBTYPE PLOT INSERTIONS
g1 <- ggplot(sub_irates, aes(x=VLoop, y=adjrate,fill=subtype)) + 
  geom_bar(stat="identity",position="dodge") + 
  labs(title="Insertions",x="Variable Loop", 
       y="Insertion Rate")+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 3))+
  scale_fill_brewer(palette="Set1")+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 12, r = 12, b = 15, l = 24, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=20,margin=margin(t = 0, r = 11, b = 0, l = 6)),
        axis.title.x=element_text(size=20,margin=margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size=16),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=25,hjust=0.5),
        legend.position="right", legend.text=element_text(size=18), legend.title=element_text(size=20)) + geom_text(aes(y=0.1,x=3 ),
                                            label="N/A", 
                                            size=6)
g1


# SUBTYPE PLOT DELETIONS 
g2 <- ggplot(sub_drates, aes(x=VLoop, y=adjrate,fill=subtype)) + 
  geom_bar(stat="identity",position="dodge") + 
  labs(title="Deletions",x="Variable Loop", 
       y="Deletion Rate")+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_fill_brewer(palette="Set1")+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 12, r = 12, b = 15, l = 24, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=20,margin=margin(t = 0, r = 11, b = 0, l = 6)),
        axis.title.x=element_text(size=20,margin=margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size=16),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=25,hjust=0.5),
        legend.position="right", legend.text=element_text(size=18), legend.title=element_text(size=20)) + geom_text(aes(y=0.1,x=3 ),
                                                                                                                    label="N/A", 
                                                                                                                    size=6)
g2

