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

ipatlist <- list()
dpatlist <- list()
iint <- list()
itip <- list()
dint <- list()
dtip <- list()

imid <- list()
dmid <- list()

for (i in 1:5){
  imid[[i]] <- numeric()
  dmid[[i]] <- numeric()
}

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

  # remove the root from the rtt length vector because it is NOT found in the reconstruction or the indel extraction (deprecated)
  #lens <- lens[-(Ntip(tre)+1)]
  res <- unname(sapply(iCSV$header, function(x){
    # this expression will return results for NODES ONLY
    # second column provides the CAPTURED TIP LABELS from within the node label
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
    return(c(lens[index], index))
  }))
  
  res <- as.data.frame(t(res))
  
  iCSV$length <- res[,1]
  dCSV$length <- res[,1]
  
  iCSV$rtt.mid <- (res[,1] + lens[tre$edge[match(res[,2], tre$edge[,2]),1]]) / 2
  dCSV$rtt.mid <- iCSV$rtt.mid
  
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
  ipatlist[[pat]] <- rbind(ipatlist[[pat]], iCSV[,2:8])
  dpatlist[[pat]] <- rbind(dpatlist[[pat]], dCSV[,2:8])
  
  # Patient-based data frame split by interior vs tip
  # iint[[pat]] <- rbind(iint[[pat]], iCSV[which(!grepl("^[^\\(\\):\n]+$", iCSV$header)),2:7])
  # itip[[pat]] <- rbind(itip[[pat]], iCSV[which(grepl("^[^\\(\\):\n]+$", iCSV$header)),2:7])
  # dint[[pat]] <- rbind(dint[[pat]], dCSV[which(!grepl("^[^\\(\\):\n]+$", dCSV$header)),2:7])
  # dtip[[pat]] <- rbind(dtip[[pat]], dCSV[which(grepl("^[^\\(\\):\n]+$", dCSV$header)),2:7])
  
  # Cumulative data frame split by interior vs tip
  iint[[file]] <- iCSV[which(!grepl("^[^\\(\\):\n]+$", iCSV$header)),2:8]
  itip[[file]]  <- iCSV[which(grepl("^[^\\(\\):\n]+$", iCSV$header)),2:8]
  dint[[file]]  <- dCSV[which(!grepl("^[^\\(\\):\n]+$", dCSV$header)),2:8]
  dtip[[file]]  <- dCSV[which(grepl("^[^\\(\\):\n]+$", dCSV$header)),2:8]
  #iprop <- idates / max(iCSV$length)
  #dprop <- ddates / max(dCSV$length)
  ires <- sapply(2:6, function(x){
    lens <- rep(iCSV$rtt.mid, iCSV[,x])
    if (length(lens) > 0 ){
      return(lens)
    }else{
      return(c())
    }
  })
  dres <- sapply(2:6, function(x){
    lens <- rep(dCSV$rtt.mid, dCSV[,x])
    if (length(lens) > 0 ){
      return(lens)
    }else{
      return(c())
    }
  })
  for (n in 1:5){
    imid[[n]] <- c(imid[[n]], ires[[n]])
    dmid[[n]] <- c(dmid[[n]], dres[[n]])
  }
  
  imaxes[count] <- max(iCSV$length,na.rm=T)
  dmaxes[count] <- max(dCSV$length,na.rm=T)
  
  # load the all.ins and all.del vectors (more efficient algorithm)
  all.ins[iseqcount:(iseqcount+sum(icounts)-1)] <- idates
  all.del[dseqcount:(dseqcount+sum(dcounts)-1)] <- ddates
  # used to maintain the vector loading algorithm above 
  iseqcount <- iseqcount + sum(icounts)
  dseqcount <- dseqcount + sum(dcounts)
  
  itotal <- itotal + nrow(iCSV)
  dtotal <- dtotal + nrow(dCSV)
  
}

require(data.table)
iint <- as.data.frame(rbindlist(iint))
itip <- as.data.frame(rbindlist(itip))
dint <- as.data.frame(rbindlist(dint))
dtip <- as.data.frame(rbindlist(dtip))


# ---- RIDGES PLOT FOR INDEL TIMINGS --------
library(ggplot2)
library(ggridges)

data <- dmid
par(mfrow=c(5,1), mar=c(3,5,1,1))
for (i in 1:4){
  hist(data[[i]], 
       breaks=30, 
       col="red", 
       xaxt='n',
       xlab="",
       main=paste0("Variable Loop ",i), 
       cex.lab=1.4,
       cex.axis=1.2,
       cex.main=1.4)
}
par(mar=c(4,5,1,1))
hist(data[[5]], 
     breaks=30, 
     col="red",
     main=paste0("Variable Loop ",5), 
     xlab="Days",
     cex.lab=1.4,
     cex.axis=1.2,
     cex.main=1.4)
axis(1, labels=T,at=seq(0, 7000, 1000), line=3)






# NEED TO FIX THIS 
vlengths <- c(72,126,105,87,33)

# COMPARISON OF INSERTION RATES INTERIOR VS TIP
int.df <- iint
tip.df <- itip
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


# ---- INSERTIONS ----
ibins <- lapply(ipatlist, function(df){
  counts <- rowSums(df[,1:5])
  dates <- rep(df[,7], counts)
  
  res <- c()
  for (i in 1:15){
    res[i] <- sum(dates > (i-1)*500 & dates < i*500)
  }
  as.data.frame(t(res))
})

ibin.df <- as.data.frame(rbindlist(ibins))
colnames(ibin.df) <- as.character(seq(0,7500,500)[-1])
ifreq <- apply(ibin.df, 2, mean)

imaxes <- imaxes[!is.na(imaxes)]
# adjust the means for the number of patients
adj.means <- mapply(function(bin, mean){
  # this is a calculation of how many data sets are still active, decreasing as less data is available
  adj.factor <- (length(imaxes) - sum(imaxes <= (bin - 500))) / length(imaxes)
  print(adj.factor)
  mean / adj.factor
}, as.numeric(colnames(ibin.df)), ifreq)


par(xpd=NA, mar=c(7,6,4,1))
barplot(adj.means, col="dodgerblue", space=0, xaxt = "n",
        #xlab="Days Since Start of Infection",
        ylab="Average Number of Insertions / Patient",
        main="Insertion Timings",
        cex.lab=1.3,cex.main=1.7)
axis(1, seq(0,15), labels=seq(0,7500,500), tick=T, line=0.5)
title(xlab="Days Since Start of Infection \n(Branch Midpoints)", line=5, cex.lab=1.3)

# ---- DELETIONS ----
dbins <- lapply(dpatlist, function(df){
  counts <- rowSums(df[,1:5])
  dates <- rep(df[,7], counts)
  
  res <- c()
  for (i in 1:15){
    res[i] <- sum(dates > (i-1)*500 & dates < i*500)
  }
  as.data.frame(t(res))
})

dbin.df <- as.data.frame(rbindlist(dbins))
colnames(dbin.df) <- as.character(seq(0,7500,500)[-1])
dfreq <- apply(dbin.df, 2, mean)

# adjust the means for the number of patients
adj.means <- mapply(function(bin, mean){
  adj.factor <- (length(dmaxes) - sum(dmaxes <= (bin - 500))) / length(dmaxes)
  print(adj.factor)
  mean / adj.factor
}, as.numeric(colnames(dbin.df)), dfreq)


dmaxes <- dmaxes[!is.na(dmaxes)]
par(xpd=NA, mar=c(7,6,4,1))
barplot(adj.means, col="red", space=0, xaxt = "n",
        #xlab="Days Since Start of Infection",
        ylab="Average Number of Deletions / Patient",
        main="Deletion Timings",
        cex.lab=1.3,cex.main=1.7)
axis(1, seq(0,15), labels=seq(0,7500,500), tick=T, line=0.5)
title(xlab="Days Since Start of Infection", line=4, cex.lab=1.3)




# --------- HISTOGRAMS (used for counts) -------------
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


# -------------
# SURVIVAL PLOT FOR INDEL TIMINGS

#imaxes <- imaxes[!is.na(imaxes)]
#dmaxes <- dmaxes[!is.na(dmaxes)]

# amalgamate the data sets 
indel.max <- data.frame(max=c(imaxes,dmaxes), status=rep(1,length(imaxes)+length(dmaxes)), type=c(rep("Insertion",length(imaxes)), rep("Deletion", length(dmaxes))))

imax <- data.frame(max=imaxes, status=rep(1,length(imaxes)))
dmax <- data.frame(max=dmaxes, status=rep(1,length(dmaxes)))

data <- imax

fit <- survfit(Surv(max,status) ~ 1, data=data)
require(survminer)
require(ggfortify)
plot <- autoplot(fit, facets=T, conf.int = F, surv.colour = "red")  + 
  labs(x="Time (Days)",
       y="Survival (%)",title = "Patient Max Dates")+
  theme(panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 42, r = 10, b = 30, l = 20, unit = "pt"),
        axis.line = element_line(colour = "black"),
        axis.title.y=element_text(size=16,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        axis.title.x=element_text(size=16,margin=margin(t = 8, r = 3, b = 0, l = 0)),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.title=element_text(size=18,hjust=0.5)),
axis.title = ,
legend.position="none")#+ geom_text(aes(y=0.4,x=3 ),
#label="N/A",
#size=6)


