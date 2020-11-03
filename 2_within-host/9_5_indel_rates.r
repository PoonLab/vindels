require(ape)
require(stringr)
require(phangorn)
require(data.table)
require(bbmle)
source("~/vindels/2_within-host/utils.r")

# INSERTION PARSING ----------
path <- "~/PycharmProjects/hiv-withinhost/"
ifolder <- Sys.glob(paste0(path,"9Indels/rep2/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"9Indels/rep2/del/*.tsv"))
sep <- "\t"
trefolder <- paste0(path,"7SampleTrees/prelim200/")

# CASE: Removed patient 56552 because could not complete Historian runs 
# CASE: Removed patients 49641 and 56549 because they are SUBTYPE B
# CASE: Removed patient 28376 and B because of very bad Rsquared value 
reg <- "56552|49641|56549|28376|B"

ifolder <- ifolder[!grepl(reg,ifolder)]
dfolder <- dfolder[!grepl(reg,dfolder)]

tally <- function(infolder){
  name <- basename(infolder)
  name <- gsub("-.+","",name)
  return (table(name))
}

findAncestor <- function(header){
  # this expression will return results for NODES ONLY
  # second column provides the CAPTURED TIP LABELS from within the node label
  header <- substr(header, 1, nchar(header)-4)
  tips <- str_match_all(header,"([^\\)\\(,\n:]+):")[[1]][,2]
  if (length(tips) == 0){
    # no colons; this means its a TIP 
    # the index in the tre$tip.label vector is the final result
    index <- match(header, tre$tip.label)
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
  return(index)
}

all.ins <- c()
all.del <- c()

maxes <- c()

iint <- list()
itip <- list()
dint <- list()
dtip <- list()


for (file in 1:length(ifolder)){
  print(file)
  filename <- basename(ifolder[file])
  full.id <- gsub("_\\d+\\.tsv$","",filename)
  
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F, sep=sep)
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F, sep=sep)
  
  iCSV$count <- unname(sapply(iCSV$indel, csvcount, delim=":"))
  dCSV$count <- unname(sapply(dCSV$indel, csvcount, delim=":"))
  
  iCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(iCSV))
  dCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(dCSV))
  
  iCSV$header <- unname(mapply(labels, iCSV$header, iCSV$pat))
  dCSV$header <- unname(mapply(labels, dCSV$header, dCSV$pat))
  
  # reads in the tree
  treename <- strsplit(filename, "\\.")[[1]][1]
  tre <- read.tree(paste0(paste0(trefolder,treename,".tree.sample")))
  
  res <- unname(sapply(iCSV$header, findAncestor)) 
  
  iCSV$length <- tre$edge.length[match(res, tre$edge[,2])]
  dCSV$length <- iCSV$length
  
  lens <- node.depth.edgelength(tre)
  # midpoint = (rtt length of tip) + (rtt length of ancestor) / 2
  iCSV$rtt.mid <- (lens[res] + lens[tre$edge[match(res, tre$edge[,2]),1]]) / 2
  dCSV$rtt.mid <- iCSV$rtt.mid
  
  # use the full tree lengths as the maximum cutoff
  maxes[file] <- max(lens,na.rm=T)
  
  # ----- TIP + INTERIOR INDEL COUNTS ---
  # Cumulative data frame split by interior vs tip
  id <- gsub("[ab]_","",full.id)   # 16362-100
  
  # this regexp matches TIP SEQUENCES (NOT CONTAINING LEFT AND RIGHT BRACKETS)
  tips <-  which(grepl("^[^\\(\\):\n]+$", iCSV$header))
  nodes <- which(!grepl("^[^\\(\\):\n]+$", iCSV$header))
  
  iCSV <- iCSV[,-c(1,2,5,6,8)]
  dCSV <- dCSV[,-c(1,2,5,6,8)]
  
  if (is.null(iint[[id]])){
    iint[[id]] <- iCSV[nodes,]
    itip[[id]]  <- iCSV[tips,]
    dint[[id]]  <- dCSV[nodes,]
    dtip[[id]]  <- dCSV[tips,]
  }else{
    iint[[id]] <- rbind(iint[[id]], iCSV[nodes,])
    itip[[id]]  <- rbind(itip[[id]], iCSV[tips,])
    dint[[id]]  <- rbind(dint[[id]] , dCSV[nodes,])
    dtip[[id]]  <- rbind(dtip[[id]], dCSV[tips,])
  }
}

# CHECKPOINT : 9_6_finished.RData 

patnames <- unname(sapply(names(iint), function(x){strsplit(x, "-")[[1]][1]}))
pat.idx <- table(sapply(ifolder, function(x){
  strsplit(basename(x), "-")[[1]][1]
}))
toRemove <- which(grepl(reg, patnames))




# checks for which patient did not fully complete all 200 replicates
table(unname(sapply(names(csv.ins), function(x)strsplit(x, "-")[[1]][1])))
iTotal <- csv.ins[-which(unname(sapply(names(csv.ins), function(x)strsplit(x, "-")[[1]][1]))== '56552')]
dTotal <- csv.del[-which(unname(sapply(names(csv.del), function(x)strsplit(x, "-")[[1]][1]))== '56552')]

rm(csv.ins)
rm(csv.del)

require(data.table)
all.ins <- as.data.frame(rbindlist(csv.ins))
all.del <- as.data.frame(rbindlist(csv.del))

#csv.ins$header <- getPat(csv.ins$header, csv.ins$Pat)
#csv.del$header <- getPat(csv.del$header, csv.del$Pat)


irtt <- list()
drtt <- list()

require(BSDA)
# RATE ANALYSIS -------------
ins.list <- list()
del.list <- list()

for (i in 1:length(iTotal)){
  iData <- iTotal[[i]]
  dData <- dTotal[[i]]
  print(i)
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
    irate <- exp(coef(ifit)[[1]])*365/median(iFinal$Vlength)
    irates <- c(irates, irate)
    #print(summary(ifit))
    
    # DELETION RATES
    dfit <- glm(dFinal$Count ~ 1, offset=log(dFinal$Date), family="poisson")
    drate <- exp(coef(dfit)[[1]])*365/median(dFinal$Vlength)
    drates <- c(drates, drate)
    #print(summary(dfit))
    
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
  
  pat <- str_split(names(iTotal)[i], "-")[[1]][1]
  run <- str_split(names(dTotal)[i], "-")[[1]][2]
  
  # contain the 20 
  ins.list[[i]] <- data.frame(pat=pat, run=run, V1=irates[1],V2=irates[2],V3=irates[3],V4=irates[4],V5=irates[5])
  del.list[[i]] <- data.frame(pat=pat, run=run,V1=drates[1],V2=drates[2],V3=drates[3],V4=drates[4],V5=drates[5])
}
ins.list <- as.data.frame(rbindlist(ins.list))
ins.list <- as.data.frame(rbindlist(del.list))

ins.final <- split(ins.list, ins.list$pat)
del.final <- split(del.list, del.list$pat)


pat <- names(ins.final)


for (i in 1:26){
  p.data <- ins.final[[i]]
  png(file=paste0("~/vindels/Figures/within-host/rates/pat/",names(ins.final)[i],".png"), width=1200,height=1000)
  par(mfrow=c(3,2))
  for (v in 3:7){
    vloop <- p.data[p.data[,v] > 10^-3,v]
    if (length(vloop) > 0 ){
      hist(vloop, main=paste0(pat[i],"-V",v-2), breaks=30, col="skyblue")
    }else{
      hist(c(0))
    }
  }
  dev.off()
}


V1 <- lapply(ins.final, function(x){if (median(x[,3] > 1e-2))x[,3]})
V2 <- lapply(ins.final, function(x){if (median(x[,4] > 1e-2))x[,4]})
V3 <- lapply(ins.final, function(x){if (median(x[,5] > 1e-2))x[,5]})
V4 <- lapply(ins.final, function(x){if (median(x[,6] > 1e-2))x[,6]})
V5 <- lapply(ins.final, function(x){if (median(x[,7] > 1e-2))x[,7]})

require(plyr)
V1 <- compact(V1)
V2 <- compact(V2)
V3 <- compact(V3)
V4 <- compact(V4)
V5 <- compact(V5)


dev.off()

ins.list <- split(ins.df, ins.df$pat)
del.list <- split(del.df, del.df$pat)

lapply(ins.list, function(x) if (mean(x$V1)>10^-2){mean(x$V1)}else{NA})


# used to randomly sample a single rate from each patient, to get a sense of the variation among patients 

res <- lapply(1:5, function(v){
  c(sapply(1:100, function(n){
    unname(unlist(lapply(ins.final, function(x){
      x[sample(1:nrow(x),1), v+2]
    })))
  }))
})
# rates <- list(V1,V2,V3,V4,V5)
# boxplot(rates)
# View(rates)


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
