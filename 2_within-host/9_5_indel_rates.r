require(bbmle)
require(stringr)
require(ape)
source("~/vindels/2_within-host/utils.r")

vloops <- c("V1","V2","V3","V4","V5")
# csvcount <- function(input){
#   commas <- str_count(input, ",")
#   if (commas > 0){
#     result <- commas + 1  
#   }else if(input == ""){
#     result <- 0
#   }else{
#     result <- 1
#   }
#   result
# }



# used for extracting condensed CSV information 
# specifically handles fields containing a comma



# INSERTION PARSING ----------
path <- "~/Lio/"
path <- "~/PycharmProjects/hiv-withinhost/"
ifolder <- Sys.glob(paste0(path,"9Indels/rep/ins/*.csv"))
dfolder <- Sys.glob(paste0(path,"9Indels/rep/del/*.csv"))
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
  pat <- str_split(filename, "_")[[1]][1]
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
  
  
  # reads in the tree
  tre <- read.tree(paste0(paste0(path,"7SampleTrees/prelim/",filename , ".tree.sample")))
  
  tip.rtt <- node.depth.edgelength(tre)[1:Ntip(tre)]
  names(tip.rtt) <- tre$tip.label
  
  tip.rows <- which(tre$edge[,2] <= Ntip(tre))
  tip.nodes <- tre$edge[tip.rows,1]
  anc.rtt <- node.depth.edgelength(tre)[tip.nodes]
  names(anc.rtt) <- tre$tip.label
  
  iheaders <- unname(sapply(iCSV$header, function(x){gsub("_\\d+$","",x)[[1]]}))
  dheaders <- unname(sapply(dCSV$header, function(x){gsub("_\\d+$","",x)[[1]]}))
  
  iCSV$mid.rtt <- (tip.rtt[iheaders] + anc.rtt[iheaders]) / 2
  dCSV$mid.rtt <- (tip.rtt[dheaders] + anc.rtt[dheaders]) / 2

  
  # adjusts the tre tip labels to match the accession numbers
  #tre$tip.label <- unname(sapply(tre$tip.label, function(x){strsplit(x,"_")[[1]][1]}))
  #tre <- read.tree(paste0(paste0(path,"7SampleTrees/prelim/",filename , ".tree.sample")))
  # retrieves branch lengths from the tree
  branches <- tre$edge.length[tre$edge[,2] <=Ntip(tre)]   
  
  # matches the branch length to each of the sequences 
  iCSV$Date <- branches[match(sub("_\\d*$","",iCSV$header), tre$tip.label)]
  dCSV$Date <- branches[match(sub("_\\d*$","",dCSV$header), tre$tip.label)]
  
  # extracts info from the indel column and puts it into two separate columns
  insInfo <- unname(sapply(iCSV$ins, extractInfo))
  insInfo <- t(insInfo)
  insInfo <- as.data.frame(insInfo)
  insInfo$V1 <- as.character(insInfo$V1)
  insInfo$V2 <- as.character(insInfo$V2)
  iCSV <- cbind(iCSV, insInfo)
  iCSV$ins <- NULL
  
  delInfo <- unname(sapply(dCSV$del, extractInfo))
  delInfo <- t(delInfo)
  delInfo <- as.data.frame(delInfo)
  delInfo$V1 <- as.character(delInfo$V1)
  delInfo$V2 <- as.character(delInfo$V2)
  dCSV <- cbind(dCSV, delInfo)
  dCSV$del <- NULL
  
  iCSV$Run <- rep(runno, nrow(iCSV))
  dCSV$Run <- rep(runno, nrow(dCSV))
  
  # creates the counts column
  iCSV$Count <- sapply(iCSV$V1, csvcount) 
  dCSV$Count <- sapply(dCSV$V1, csvcount)
  
  iCSV$Pat <- paste0(str_split(rep(pat, nrow(iCSV)),"-")[[1]][1],"-", iCSV$Run)
  dCSV$Pat <- paste0(str_split(rep(pat, nrow(dCSV)),"-")[[1]][1] ,"-", dCSV$Run)
  
  #rearrange
  iCSV <- iCSV[c(1,2,3,6,12,8,9,10,7,4,5,13,11)]
  dCSV <- dCSV[c(1,2,3,6,12,8,9,10,7,4,5,13,11)]
  
  c.headers <- c("header","Vloop", "Vlength","Subtype", "Count","Date", "Seq", "Pos","mid.rtt","Vseq", "Anc", "Pat", "Run")
  colnames(iCSV) <- c.headers
  colnames(dCSV) <- c.headers
  

  
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

csv.ins$header <- getPat(csv.ins$header, csv.ins$Pat)
csv.del$header <- getPat(csv.del$header, csv.del$Pat)

iTotal <- split(csv.ins, csv.ins$Pat)
dTotal <- split(csv.del, csv.del$Pat)


require(BSDA)
# RATE ANALYSIS -------------
ins.df <- data.frame()
del.df <- data.frame()

#median(as.numeric(all.ins[all.ins$Vloop==2,3])) # used to det ermine the median lengths of the variable loops
vlengths <- c(72,126,105,87,33)
all.df <- data.frame()

irtt <- list()
drtt <- list()

all.ins <- list()
all.del <- list()
for (name in unique(csv.ins$Pat)){
  iData <- iTotal[[name]]
  dData <- dTotal[[name]]
  
  irates <- c()
  drates <- c()
  
  for (vloop in 1:5){
    itemp <- iData[iData$Vloop==vloop,]
    dtemp <- dData[dData$Vloop==vloop,]
    
    iFinal <- itemp[itemp$Date < 325 & itemp$Count < 3,]
    dFinal <- dtemp[dtemp$Date < 325 & dtemp$Count < 3,]
    
    # ADDED THIS FOR RTT ANALYSIS
    # used to retrieve the midpoint of each branch length along which an indel occurred
    irtt[[as.character(vloop)]] <- c(irtt[[as.character(vloop)]],iFinal[iFinal$Count>0,'mid.rtt'])
    drtt[[as.character(vloop)]] <- c(drtt[[as.character(vloop)]], dFinal[dFinal$Count>0,'mid.rtt'])
    #print(nrow(current) - nrow(iFinal))
    
    # INSERTION RATES 
    ifit <- glm(iFinal$Count ~ 1, offset=log(iFinal$Date), family="poisson")
    irate <- exp(coef(ifit)[[1]])*365/vlengths[vloop]
    irates <- c(irates, irate)
    print(summary(ifit))
    
    # DELETION RATES
    dfit <- glm(dFinal$Count ~ 1, offset=log(dFinal$Date), family="poisson")
    drate <- exp(coef(dfit)[[1]])*365/vlengths[vloop]
    drates <- c(drates, drate)
    print(summary(dfit))
    
    # COUNTING NUCLEOTIDES: change the formula to this : nchar(gsub(",","",iFinal$Seq))*iFinal$Count ~ 1
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
  
  pat <- str_split(name, "-")[[1]][1]
  run <- str_split(name, "-")[[1]][2]
  
  # contain the 20 
  ins.df <- rbind(ins.df, data.frame(pat=pat, run=run, V1=irates[1],V2=irates[2],V3=irates[3],V4=irates[4],V5=irates[5]))
  del.df <- rbind(del.df, data.frame(pat=pat, run=run,V1=drates[1],V2=drates[2],V3=drates[3],V4=drates[4],V5=drates[5]))
}

ins.list <- split(ins.df, ins.df$pat)
del.list <- split(del.df, del.df$pat)

lapply(ins.list, function(x) if (mean(x$V1)>10^-2){mean(x$V1)}else{NA})

require(data.table)
iTotal2 <- rbindlist(iTotal)
iTotal3 <- split(iTotal2, iTotal2$Vloop)


# looking at rtt branch lengths and when indels tend to occur
ires <- c()
dres <- c()
for (i in 1:length(irtt)){
  ires <- c(ires, median(irtt[[i]]))
  dres <- c(dres, median(drtt[[i]]))
}


# --------------------------------
# INSERTIONS VS DELETIONS 
# Comparing insertion and deletion rates to each other 
# Shown that Deletions are singificantly higher than insertions
ir <- c(ins.df[,1],ins.df[,2],ins.df[,3],ins.df[,4],ins.df[,5])
dr <- c(del.df[,1],del.df[,2],del.df[,3],del.df[,4],del.df[,5])
wilcox.test(ir,dr,paired=T)




# ----------------------------
#  BETWEEN HOST COMPARISON 
# Comparing combined within host indel rates to between host indel rates 
subtypes <- c("A1", "B", "C")
ins.sub <- list()
del.sub <- list()

all.df <- data.frame()
for (run in 1:20){
  iData <- iTotal[[run]]
  dData <- dTotal[[run]]
  
  for (sub in subtypes){
    ins.sub[[sub]] <- rbind(ins.sub[[sub]], iData[iData$Subtype==sub,])
    del.sub[[sub]] <- rbind(del.sub[[sub]], dData[dData$Subtype==sub,])
  }
}

isub <- data.frame()
dsub <- data.frame()

for (sub in subtypes){
  iData <- ins.sub[[sub]]
  dData <- del.sub[[sub]]
  
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
  }
  
  irates <- irates*10^3
  drates <- drates*10^3
  isub <- rbind(isub, data.frame(t(irates)))
  dsub <- rbind(dsub, data.frame(t(drates)))
}

rownames(isub) <- subtypes
rownames(dsub) <- subtypes
colnames(isub) <- vloops
colnames(dsub) <- vloops

indel.df <- isub + dsub
within <- c()

for (row in 1:3){
  within <- c(within, as.double(indel.df[row,]))
}

# max.llh is retrieved from indel_analysis.r in the ~/indelrates/ repo
btwrates <- max.llh[-c(6:15,26:35),]
wilcox.test(within, btwrates$adj.rate, paired=T)



# wilcoxon statistical test comparing median indel rates 
indel.df <- ins.df + del.df
btw.rates <- c()
for (i in 1:5){
  btw.rates <- c(btw.rates, median(max.llh[max.llh$vloop==i,"adj.rate"]))
}
wth.rates <- as.double(apply(indel.df, 2, median))
wilcox.test(btw.rates,wth.rates, paired=T)



# --------------------------
# BOOTSTRAPS - over 20 replicates (simple)
ins.rep <- matrix(nrow=2, ncol=5)
del.rep <- matrix(nrow=2, ncol=5)
for (v in 1:5){
  itemp <- ins.df[,v]
  dtemp <- del.df[,v]
  
  imeds <- sapply(1:1000, function(x){
    rand <- sample(1:20, 20, replace=T)
    median(itemp[rand])
  })
  
  dmeds <- sapply(1:1000, function(x){
    rand <- sample(1:20, 20, replace=T)
    median(dtemp[rand])
  })
  
  ins.rep[,v] <- quantile(imeds, c(0.025,0.975))
  del.rep[,v] <- quantile(dmeds, c(0.025,0.975))
}



# ---------------
# STANDARD INDEL RATES 

insrates <- data.frame(VLoop=vloops, iRate=irates, AdjRate=irates*10^3)
delrates <- data.frame(VLoop=vloops, dRate=drates, AdjRate=drates*10^3)

insrates <- data.frame(vloop=vloops,
                       rate=apply(ins.df, 2, median),
                       lower=ins.rep[1,],
                       upper=ins.rep[2,])
                       #lower=apply(ins.df,2,function(x){quantile(x, c(0.025,0.975))[1]}), 
                       #upper=apply(ins.df,2,function(x){quantile(x, c(0.025,0.975))[2]}))
delrates <- data.frame(vloop=vloops,
                       rate=apply(del.df, 2, median), 
                       lower=del.rep[1,],
                       upper=del.rep[2,])
                       #lower=apply(del.df,2,function(x){quantile(x, c(0.025,0.975))[1]}), 
                       #upper=apply(del.df,2,function(x){quantile(x, c(0.025,0.975))[2]}))








# MIDPOINTS
# Used for looking at RTT midpoints (time when indels occur)
insrates <- data.frame(vloop=vloops,
                       rate=unname(unlist(lapply(irtt, median))), 
                       lower=unname(unlist(lapply(irtt,function(x){quantile(x, c(0.25,0.75))[1]}))), 
                       upper=unname(unlist(lapply(irtt,function(x){quantile(x, c(0.25,0.75))[2]}))))
delrates <- data.frame(vloop=vloops,
                       rate=unname(unlist(lapply(drtt, median))), 
                       lower=unname(unlist(lapply(drtt,function(x){quantile(x, c(0.25,0.75))[1]}))), 
                       upper=unname(unlist(lapply(drtt,function(x){quantile(x, c(0.25,0.75))[2]}))))


print(insrates)
print(delrates)
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
       y=expression(paste("       Insertion Rate \n (Events/Nt/Year x  ",10^-3 ,")", sep = "")))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) + # breaks=c(300,600,900,1200,1500)) +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 42, r = 10, b = 4, l = 20, unit = "pt"),
        axis.line = element_line(colour = "black"),
        axis.title.y=element_text(size=18,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        axis.title.x=element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        legend.position="none")#+ geom_text(aes(y=0.4,x=3 ),
                                           #label="N/A",
                                           #size=6)
#g1


g2 <- ggplot(delrates, aes(x=vloop, y=rate,width=0.8)) +
  geom_bar(colour="black", stat="identity",fill="firebrick1",position="dodge",show.legend=F) +
  geom_errorbar(aes(ymax = delrates$upper, ymin = delrates$lower),
                width = 0.25, size=1.1) +
  geom_errorbar(aes(ymax = delrates$upper, ymin = delrates$lower),
                width = 0.25, size=1.1) +
  labs(x="Variable Loop",
       y="Deletion Rate")+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_y_reverse(lim=c(6,0)) + #breaks=c(300,600,900,1200,1500))+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = -2, r = 10, b = 8, l = 26, unit = "pt"),
        axis.line = element_line(colour = "black"),
        axis.title.y=element_text(size=18,margin=margin(t = 0, r = 11, b = 0, l = 6)),
        axis.title.x=element_text(size=18,margin=margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size=16),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=14),
        legend.position="none")
#+ geom_text(aes(y=0.5,x=3 ),label="N/A", size=6)

multiplot(g1,g2)
# 
# 
# 
# # SUBTYPE STRATIFICATION 
# # ----------------------------------------------
# ins <- list()
# ins$B <- csv.ins[csv.ins$Subtype=="B",]
# ins$C <- csv.ins[csv.ins$Subtype=="C",]
# del <- list()
# del$B <- csv.del[csv.del$Subtype=="B",]
# del$C <- csv.del[csv.del$Subtype=="C",]
# 
# sub_irates <- data.frame()
# sub_drates <- data.frame()
# 
# for (i in 1:2){
#   
#   iData <- ins[[i]]
#   dData <- del[[i]]
#   irates <- c()
#   drates <- c()
#   
#   st <- c("B", "C")
#   
#   for (vloop in 1:5){
#     iFinal <- iData[iData$Vloop==vloop & iData$Date < 325 & iData$Count < 2,]
#     dFinal <- dData[dData$Vloop==vloop & dData$Date < 325 & dData$Count < 2,]
#     #print(nrow(current) - nrow(iFinal))
#     
#     ifit <- glm(iFinal$Count ~ 1, offset=log(iFinal$Date), family="poisson")
#     irate <- exp(coef(ifit)[[1]])*365/vlengths[vloop]
#     irates <- c(irates, irate)
#     print(summary(ifit))
#     
#     dfit <- glm(dFinal$Count ~ 1, offset=log(dFinal$Date), family="poisson")
#     drate <- exp(coef(dfit)[[1]])*365/vlengths[vloop]
#     drates <- c(drates, drate)
#     print(summary(dfit))
#     
#     #print(1 - (fit$deviance/fit$null.deviance))
#     #all.df <- rbind(all.df, iFinal)
#     #EDA(residuals(fit))
#     
#     #to.remove <- c(order(residuals(fit), decreasing = T)[1:52])
#     #iFinal2 <- iFinal[-to.remove,]
#     #fit2 <- glm(iFinal2$Count ~ iFinal2$Date, family="poisson")
#     #EDA(residuals(fit2))
#   }
#   sub_irates <- rbind(sub_irates, data.frame(rate=irates, adjrate=(irates*10^3),subtype=rep(st[[i]],5), VLoop=vloops))
#   sub_drates <- rbind(sub_drates, data.frame(rate=drates, adjrate=(drates*10^3),subtype=rep(st[[i]],5), VLoop=vloops))
# }
# 
# 
# # SUBTYPE PLOT INSERTIONS
# g1 <- ggplot(sub_irates, aes(x=VLoop, y=adjrate,fill=subtype)) + 
#   geom_bar(stat="identity",position="dodge") + 
#   labs(title="Insertions",x="Variable Loop", 
#        y="Insertion Rate")+
#   #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
#   scale_y_continuous(expand = c(0, 0),limits = c(0, 3))+
#   scale_fill_brewer(palette="Set1")+
#   theme(panel.grid.major.y = element_line(color="black",size=0.3),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.spacing=unit(1, "mm"),
#         #panel.background=element_rect(fill="gray88",colour="white",size=0),
#         plot.margin =margin(t = 12, r = 12, b = 15, l = 24, unit = "pt"),
#         axis.line = element_line(colour = "black"), 
#         axis.title.y=element_text(size=20,margin=margin(t = 0, r = 11, b = 0, l = 6)),
#         axis.title.x=element_text(size=20,margin=margin(t = 5, r = 0, b = 0, l = 0)),
#         strip.text.x = element_text(size=16),
#         axis.text.x=element_text(size=18),
#         axis.text.y=element_text(size=18),
#         plot.title = element_text(size=25,hjust=0.5),
#         legend.position="right", legend.text=element_text(size=18), legend.title=element_text(size=20)) + geom_text(aes(y=0.1,x=3 ),
#                                             label="N/A", 
#                                             size=6)
# g1
# 
# 
# # SUBTYPE PLOT DELETIONS 
# g2 <- ggplot(sub_drates, aes(x=VLoop, y=adjrate,fill=subtype)) + 
#   geom_bar(stat="identity",position="dodge") + 
#   labs(title="Deletions",x="Variable Loop", 
#        y="Deletion Rate")+
#   #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
#   scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
#   scale_fill_brewer(palette="Set1")+
#   theme(panel.grid.major.y = element_line(color="black",size=0.3),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.spacing=unit(1, "mm"),
#         #panel.background=element_rect(fill="gray88",colour="white",size=0),
#         plot.margin =margin(t = 12, r = 12, b = 15, l = 24, unit = "pt"),
#         axis.line = element_line(colour = "black"), 
#         axis.title.y=element_text(size=20,margin=margin(t = 0, r = 11, b = 0, l = 6)),
#         axis.title.x=element_text(size=20,margin=margin(t = 5, r = 0, b = 0, l = 0)),
#         strip.text.x = element_text(size=16),
#         axis.text.x=element_text(size=18),
#         axis.text.y=element_text(size=18),
#         plot.title = element_text(size=25,hjust=0.5),
#         legend.position="right", legend.text=element_text(size=18), legend.title=element_text(size=20)) + geom_text(aes(y=0.1,x=3 ),
#                                                                                                                     label="N/A", 
#                                                                                                                     size=6)
