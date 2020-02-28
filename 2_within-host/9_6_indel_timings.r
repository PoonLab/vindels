source("~/vindels/2_within-host/utils.r")

# INSERTION PARSING ----------
path <- "~/Lio/"
path <- "~/PycharmProjects/hiv-withinhost/"
ifolder <- Sys.glob(paste0(path,"9Indels/rep/wholetree/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"9Indels/rep/wholetree/del/*.tsv"))

require(stringr)
require(phangorn)
require(data.table)
require(bbmle)

vlist <- list(V1=numeric(),V2=numeric(),V3=numeric(),V4=numeric(),V5=numeric())
all.ins <- c()
all.del <- c()
count <- 0
iseqcount <- 1
dseqcount <- 1
imaxes <- c()
dmaxes <- c()
itotal <- 0
dtotal <- 0

#ipatlist <- list()
#dpatlist <- list()
iint <- data.frame()
itip <- data.frame()
dint <- data.frame()
dtip <- data.frame()
for (file in 1:length(ifolder)){
  print(file)
  filename <- strsplit(basename(ifolder[file]),"\\.")[[1]][1]
  pat <- gsub("_\\d+$","",filename)
  runno <- strsplit(filename, "_")[[1]][2]
  count <- count + 1
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F, sep='\t')
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F, sep='\t')
  
  if (grepl("B.+$", filename) || grepl("OS.+$", filename) || grepl("G.+$", filename) || grepl("56549.+$", filename)){
    next
  }
  
  for (i in 2:6){
    res <- unname(sapply(iCSV[,i], function(x){csvcount(x,":")}))
    iCSV[,i] <- res
    res <- unname(sapply(dCSV[,i], function(x){csvcount(x,":")}))
    dCSV[,i] <- res
  }
  
  # reads in the tree
  tre <- read.tree(paste0(paste0(path,"7SampleTrees/prelim/",filename,".tree.sample")))
  lens <- node.depth.edgelength(tre)    # [(length(tre$tip.label)+1):(length(tre$edge.length)+1)]  #used if you want to only access internal nodes and not tips

  # remove the root from the rtt length vector because it is NOT found in the reconstruction or the indel extraction
  #lens <- lens[-(Ntip(tre)+1)]
  res <- unname(sapply(iCSV$header, function(x){
    tips <- str_match_all(x,"([^\\)\\(,\n:]+):")[[1]][,2]
    if (length(tips) == 0){
      # the index in the tre$tip.label vector is the final result
      index <- match(x, tre$tip.label)
    }else{
      desc <- Descendants(tre)

      # find the numeric labels of all extracted tips 
      matches <- match(tips, tre$tip.label)
      
      # find the SINGLE node in the descendants list that contains the exact same subset of tips
      index <- which(sapply(desc, function(x){ifelse(length(x) == length(matches) && all(x==matches),T,F)}))
    }
    if (length(index)!=1){
      return(paste0("PROBLEM:",as.character(index)))
    }
    return(lens[index])
  }))
  
  iCSV$length <- res
  dCSV$length <- res
  
  icounts <- rowSums(iCSV[,2:6])
  dcounts <- rowSums(dCSV[,2:6])
  
  if (sum(icounts)==0){
    next
  }
  if (sum(dcounts)==0){
    next
  }
  
  idates <- rep(iCSV$length, icounts)
  ddates <- rep(dCSV$length, dcounts)
  
  if (is.na(max(idates))){
    print(idates)
  }
  if (is.na(max(ddates))){
    print(ddates)
  }
  # Regular patient-based data frame
  #ipatlist[[pat]] <- rbind(ipatlist[[pat]], iCSV[,2:7])
  #dpatlist[[pat]] <- rbind(dpatlist[[pat]], dCSV[,2:7])
  
  # Patient-based data frame split by interior vs tip
  # iint[[pat]] <- rbind(iint[[pat]], iCSV[which(!grepl("^[^\\(\\):\n]+$", iCSV$header)),2:7])
  # itip[[pat]] <- rbind(itip[[pat]], iCSV[which(grepl("^[^\\(\\):\n]+$", iCSV$header)),2:7])
  # dint[[pat]] <- rbind(dint[[pat]], dCSV[which(!grepl("^[^\\(\\):\n]+$", dCSV$header)),2:7])
  # dtip[[pat]] <- rbind(dtip[[pat]], dCSV[which(grepl("^[^\\(\\):\n]+$", dCSV$header)),2:7])
  
  # Cumulative data frame split by interior vs tip
  iint <- rbind(iint, iCSV[which(!grepl("^[^\\(\\):\n]+$", iCSV$header)),2:7])
  itip <- rbind(itip, iCSV[which(grepl("^[^\\(\\):\n]+$", iCSV$header)),2:7])
  dint <- rbind(dint, dCSV[which(!grepl("^[^\\(\\):\n]+$", dCSV$header)),2:7])
  dtip <- rbind(dtip, dCSV[which(grepl("^[^\\(\\):\n]+$", dCSV$header)),2:7])
  #iprop <- idates / max(iCSV$length)
  #dprop <- ddates / max(dCSV$length)
  
  imaxes[count] <- max(iCSV$length)
  dmaxes[count] <- max(dCSV$length)
  
  # load the all.ins and all.del vectors (more efficient algorithm)
  all.ins[iseqcount:(iseqcount+sum(icounts)-1)] <- idates
  all.del[dseqcount:(dseqcount+sum(dcounts)-1)] <- ddates
  # used to maintain the vector loading algorithm above 
  iseqcount <- iseqcount + sum(icounts)
  dseqcount <- dseqcount + sum(dcounts)
  
  itotal <- itotal + nrow(iCSV)
  dtotal <- dtotal + nrow(dCSV)
  
}

# NEED TO FIX THIS 
vlengths <- c(72,126,105,87,33)

# COMPARISON OF INSERTION RATES INTERIOR VS TIP
int.df <- dint
tip.df <- dtip
# int.out <- data.frame()
# tip.out <- data.frame()
# for (run in 1:20){
#   int.df <- iint[[run]]
#   tip.df <- itip[[run]]
#   
  irates2 <- c()
  trates2 <- c()
  for (i in 1:5){
    tip.fit <- glm(tip.df[,i] ~ 1, offset=log(tip.df[,"length"]), family="poisson")
    int.fit <- glm(int.df[,i] ~ 1, offset=log(int.df[,"length"]), family="poisson")
    
    print(summary(tip.fit))
    print(summary(int.fit))
    
    tip.rate <- exp(coef(tip.fit)[[1]])*365/vlengths[i]
    int.rate <- exp(coef(int.fit)[[1]])*365/vlengths[i]
    
    irates2[i] <- int.rate
    trates2[i] <- tip.rate
    
  }
#   irates <- irates*10^3
#   trates <- trates*10^3
#   
#   int.out <- rbind(int.out, data.frame(V1=irates[1],V2=irates[2],V3=irates[3],V4=irates[4],V5=irates[5]))
#   tip.out <- rbind(tip.out, data.frame(V1=trates[1],V2=trates[2],V3=trates[3],V4=trates[4],V5=trates[5]))
# }

  comb1 <- data.frame(rate=irates, vloop=c("V1","V2","V3","V4","V5"), id=rep("Interior",5))
  
require(RColorBrewer)
#pal <- c("gray28", "blue4",  'tomato', 'dodgerblue',  'red',  "skyblue", 'darkred' )

require(ggplot2)
dplot <- ggplot() + 
  geom_bar(aes(x=vloop, y=rate, fill=id), data=comb2, stat='identity', position="dodge") + 
  scale_fill_manual(values=c("red","blue"))+
  labs(x="Variable Loop", 
       y="Deletion Rate (Events/Nt/Year)", title="Deletion Rates") +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 1.3, r = 1, b = 0.7, l = 1.5, unit = "cm"),
        axis.line = element_line(colour = "black"), 
        axis.title=element_text(size=18,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        axis.text = element_text(size=16, colour="black"),
        plot.title = element_text(size=22, hjust = 0.5),
        legend.text=element_text(size=16), 
        legend.background=element_rect(colour="black"),
        legend.title=element_text(size=18))
dplot



# COMPARISON OF DELETION RATES INTERIOR VS TIP 


#imaster <- as.data.frame(rbindlist(ipatlist))
#dmaster <- as.data.frame(rbindlist(dpatlist))

# INEFFICIENT PLEASE REWRITE
# -------------------------------
idates <- lapply(ipatlist, function(list){
  counts <- rowSums(list[,1:5])
  dates <- rep(list[,6], counts)
  dates
})

ibins <- lapply(idates, function(x){
  res <- c()
  for (i in 1:15){
    res[i] <- sum(x > (i-1)*500 & x < i*500)
  }
  as.data.frame(t(res))
})

ddates <- lapply(dpatlist, function(list){
  counts <- rowSums(list[,1:5])
  dates <- rep(list[,6], counts)
  dates
})

dbins <- lapply(ddates, function(x){
  res <- c()
  for (i in 1:15){
    res[i] <- sum(x > (i-1)*500 & x < i*500)
  }
  as.data.frame(t(res))
})


# INSERTIONS 
ibins <- rbindlist(ibins)
colnames(ibins) <- as.character(seq(0,7500,500)[-1])
ifreq <- apply(ibins, 2, mean)

imaxes <- imaxes[!is.na(imaxes)]
newimaxes <- imaxes/500
par(xpd=NA, mar=c(7,6,4,1))
barplot(ifreq, col="dodgerblue", space=0, xaxt = "n",
        #xlab="Days Since Start of Infection",
        ylab="Average Number of Insertions / Patient",
        main="Insertion Timings",
        cex.lab=1.3,cex.main=1.7)
arrows(newimaxes, 0, newimaxes, -0.2, length=0)
axis(1, seq(0,15), labels=seq(0,7500,500), tick=T, line=0.5)
title(xlab="Days Since Start of Infection", line=4, cex.lab=1.3)

# DELETIONS
dbins <- rbindlist(dbins)
colnames(dbins) <- as.character(seq(0,7500,500)[-1])
dfreq <- apply(dbins, 2, mean)

dmaxes <- dmaxes[!is.na(dmaxes)]
newdmaxes <- dmaxes/500
par(xpd=NA, mar=c(7,6,4,1))
barplot(dfreq, col="dodgerblue", space=0, xaxt = "n",
        #xlab="Days Since Start of Infection",
        ylab="Average Number of Deletions / Patient",
        main="Deletion Timings",
        cex.lab=1.3,cex.main=1.7)
arrows(newdmaxes, 0, newdmaxes, -0.5, length=0)
axis(1, seq(0,15), labels=seq(0,7500,500), tick=T, line=0.5)
title(xlab="Days Since Start of Infection", line=4, cex.lab=1.3)







# HISTOGRAMS (used for counts)
# ----------------------
imaxes <- imaxes[!is.na(imaxes)]
par(mar=c(5,5,5,2),xpd=F)
caxis=1.3
clab=1.4
cmain=1.5

hist(all.ins[all.ins<3000], 
     col='red',cex.lab=clab,
     main="Timing of Insertions",
     cex.axis=caxis, 
     cex.main=cmain, 
     xlab="Normalized Time of Infection")
#hist(imaxes[imaxes<3000], 
     #col='blue',add=T,breaks=seq(0,3000,200))

arrows(imaxes, 0, imaxes, -90, length=0)

hist(all.del[all.del<3000], 
     col='red',cex.lab=clab,
     main="Timing of Deletions",
     cex.axis=caxis, 
     cex.main=cmain, 
     xlab="Normalized Time of Infection")
#hist(imaxes[imaxes<3000], 
#col='blue',add=T,breaks=seq(0,3000,200))

arrows(dmaxes, 0, dmaxes, -300, length=0)

# LINE PLOTS 
# ----------------------------
par(mar=c(5,5,5,2))
caxis=1.3
clab=1.4
cmain=1.5

hist(all.ins,
     breaks=40,
     col='red',cex.lab=clab,
     main="Timing of Insertions",
     cex.axis=caxis, 
     cex.main=cmain, 
     xlab="Normalized Time")

hist(all.del, 
     breaks=20,
     col='red',cex.lab=clab,
     main="Timing of Deletions",
     cex.axis=caxis, 
     cex.main=cmain, 
     xlab="Days Since Start of Infection")





