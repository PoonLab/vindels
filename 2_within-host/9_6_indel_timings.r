source("~/vindels/2_within-host/utils.r")

# INSERTION PARSING ----------
path <- "~/Lio/"
path <- "~/PycharmProjects/hiv-withinhost/"
ifolder <- Sys.glob(paste0(path,"9Indels/rep/wholetree/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"9Indels/rep/wholetree/del/*.tsv"))

require(stringr)
require(phangorn)

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
  
  ipatlist[[pat]] <- rbind(ipatlist[[pat]], iCSV[,2:7])
  dpatlist[[pat]] <- rbind(dpatlist[[pat]], dCSV[,2:7])
  
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

imaster <- rbindlist(ipatlist)
dmaster <- rbindlist(dpatlist)

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





