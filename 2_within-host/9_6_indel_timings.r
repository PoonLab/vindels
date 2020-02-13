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
for (file in 1:length(ifolder)){
  print(file)
  filename <- strsplit(basename(ifolder[file]),"\\.")[[1]][1]
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
  
  iprop <- idates / max(iCSV$length)
  dprop <- ddates / max(dCSV$length)
  
  imaxes[count] <- max(iCSV$length)
  dmaxes[count] <- max(dCSV$length)
  
  # load the all.ins and all.del vectors (more efficient algorithm)
  all.ins[iseqcount:(iseqcount+sum(icounts)-1)] <- iprop
  all.del[dseqcount:(dseqcount+sum(dcounts)-1)] <- dprop
  # used to maintain the vector loading algorithm above 
  iseqcount <- iseqcount + sum(icounts)
  dseqcount <- dseqcount + sum(dcounts)
  
}

# HISTOGRAMS (used for counts)
# ----------------------

par(mar=c(5,5,5,2))
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





