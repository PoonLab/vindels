source("~/vindels/2_within-host/utils.r")

# INSERTION PARSING ----------
path <- "~/Lio/"
path <- "~/PycharmProjects/hiv-withinhost/"
ifolder <- Sys.glob(paste0(path,"9Indels/rep2/wholetree/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"9Indels/rep2/wholetree/del/*.tsv"))

require(stringr)
require(phangorn)
require(data.table)
require(bbmle)

vlist <- list(V1=numeric(),V2=numeric(),V3=numeric(),V4=numeric(),V5=numeric())
all.ins <- c()
all.del <- c()
count <- 0

maxes <- c()
itotal <- 0
dtotal <- 0

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
  full.id <- gsub("_\\d+$","",filename)

  count <- count + 1
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F, sep='\t')
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F, sep='\t')
  
  # if (grepl("B.+$", filename) || grepl("OS.+$", filename) || grepl("G.+$", filename) || grepl("56549.+$", filename)){
  #   next
  # }
  
  for (i in 2:6){
    res <- unname(sapply(iCSV[,i], function(x){csvcount(x,":")}))
    iCSV[,i] <- res
    res <- unname(sapply(dCSV[,i], function(x){csvcount(x,":")}))
    dCSV[,i] <- res
  }
  
  # reads in the tree
  tre <- read.tree(paste0(paste0(path,"7SampleTrees/prelim200/",filename,".tree.sample")))
      # [(length(tre$tip.label)+1):(length(tre$edge.length)+1)]  #used if you want to only access internal nodes and not tips
  
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
    return(c(index))
  }))
  
  
  
  iCSV$length <- tre$edge.length[match(res, tre$edge[,2])]
  dCSV$length <- iCSV$length
  
  lens <- node.depth.edgelength(tre)
  iCSV$rtt.mid <- (lens[res] + lens[tre$edge[match(res, tre$edge[,2]),1]]) / 2
  dCSV$rtt.mid <- iCSV$rtt.mid
  
  icounts <- rowSums(iCSV[,2:6])
  dcounts <- rowSums(dCSV[,2:6])
  
  if (sum(icounts)==0){
    next
  }
  if (sum(dcounts)==0){
    next
  }
  
  # ----- INDEL TIMINGS ---
  
  idates <- rep(iCSV$rtt.mid, icounts)
  ddates <- rep(dCSV$rtt.mid, dcounts)
  
  if (is.na(max(idates))){
    print(idates)
  }
  if (is.na(max(ddates))){
    print(ddates)
  }

  imid[[full.id]] <- idates
  dmid[[full.id]] <- ddates
  
  # these remain 'lengths' and not 'rtt.mid's because I need to use the full tree length as the maximum cutoff
  maxes[count] <- max(lens,na.rm=T)
  
  # load the all.ins and all.del vectors (more efficient algorithm)
  #all.ins[iseqcount:(iseqcount+sum(icounts)-1)] <- idates
  #all.del[dseqcount:(dseqcount+sum(dcounts)-1)] <- ddates
  # used to maintain the vector loading algorithm above 
  #iseqcount <- iseqcount + sum(icounts)
  #dseqcount <- dseqcount + sum(dcounts)
  
  # ----- TIP + INTERIOR INDEL COUNTS ---
  # Cumulative data frame split by interior vs tip
  id <- gsub("[ab]_","",full.id)   # 16362-100
  
  tips <-  which(grepl("^[^\\(\\):\n]+$", iCSV$header))
  nodes <- which(!grepl("^[^\\(\\):\n]+$", iCSV$header))
  if (is.null(iint[[id]])){
    iint[[id]] <- iCSV[nodes,2:8]
    itip[[id]]  <- iCSV[tips,2:8]
    dint[[id]]  <- dCSV[nodes,2:8]
    dtip[[id]]  <- dCSV[tips,2:8]
  }else{
    iint[[id]] <- rbind(iint[[id]], iCSV[nodes,2:8])
    itip[[id]]  <- rbind(itip[[id]], iCSV[tips,2:8])
    dint[[id]]  <- rbind(dint[[id]] , dCSV[nodes,2:8])
    dtip[[id]]  <- rbind(dtip[[id]], dCSV[tips,2:8])
  }
  
  itotal <- itotal + nrow(iCSV)
  dtotal <- dtotal + nrow(dCSV)
  
}
# determine which patients did not complete fully 
pat.idx <- unname(sapply(names(iint), function(x){strsplit(x, "-")[[1]][1]}))
table(pat.idx)
toRemove <- which(pat.idx == "F" | pat.idx == "H" | pat.idx == "QJ" | pat.idx == "56552")

iint <- iint[-toRemove]
itip <- itip[-toRemove]
dint <- dint[-toRemove]
dtip <- dtip[-toRemove]


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

all.rates <- list('1'=list(),'2'=list(),'3'=list(),'4'=list())


for (i in 1:length(iint)){
  print(i)
  rate.list <- list('1'=c(),'2'=c(),'3'=c(),'4'=c())
  
  data.list <- list(ins.tip = itip[[i]], 
                    ins.node = iint[[i]], 
                    del.tip = dtip[[i]],
                    del.node = dint[[i]])
  for (j in 1:5){

    fit.list <- lapply(data.list, function(list){
      tryCatch(
        {
          glm(list[,j] ~ 1, offset=log(list[,"length"]), family="poisson")
        },
        error=function(cond) {
          return(NA)
        })
    })
    # in.fit <- glm(ins.node[,i] ~ 1, offset=log(ins.node[,"length"]), family="poisson")
    # dt.fit <- glm(del.tip[,i] ~ 1, offset=log(del.tip[,"length"]), family="poisson")
    # dn.fit <- glm(del.node[,i] ~ 1, offset=log(del.node[,"length"]), family="poisson")
    
    rates <- unname(sapply(fit.list, function(x){
      if (all(is.na(x))){
        return(NA)
      }
      rate <- exp(coef(x)[[1]])*365/vlengths[j]
      if (rate > 10^-8){
        rate
      }else{
        NA
      }
    }))
    
    for (x in 1:4){
      rate.list[[x]][j] <- rates[x]
    }

    
    # in.rates[j] <- exp(coef(in.fit)[[1]])*365/vlengths[i]
    # dt.rates[j] <- exp(coef(dt.fit)[[1]])*365/vlengths[i]
    # dn.rates[j] <- exp(coef(dn.fit)[[1]])*365/vlengths[i]
    
  }
  rate.list <- lapply(rate.list, function(x){
    x * 10^3
  })
  patid <- strsplit(names(iint)[i], "-")[[1]][1]
  runno <- strsplit(names(iint)[i], "-")[[1]][2]
  
  for (y in 1:4){
    all.rates[[y]][[i]] <- data.frame(pat=patid, 
                                 runno=runno, 
                                 V1=rate.list[[y]][1],
                                 V2=rate.list[[y]][2],
                                 V3=rate.list[[y]][3],
                                 V4=rate.list[[y]][4],
                                 V5=rate.list[[y]][5])
  }
}

all.rates <- lapply(1:4, function(x){
  df <- as.data.frame(rbindlist(all.rates[[x]]))
  split(df, df[,"pat"])
})


vloop <- lapply(1:5, function(w){
  lapply(1:4, function(x){
    compact(lapply(all.rates[[x]], function(y){
      nas <- sum(is.na(y[,w+2]))
      if (nas > 100){
        NULL
      }else{
        y[,w+2]
      }
    }))
  })
})

vloop.mat <- lapply(1:5, function(w){
  lapply(1:4, function(x){
    data <- vloop[[w]][[x]]
    if (length(data) == 0){
      return(NA)
    }else{
      mat <- matrix(nrow=200, ncol=length(data))
      
      # load the matrix 
      for (i in 1:ncol(mat)){
        mat[,i] <- data[[i]]
      }
      toRemove <- c()
      for (j in 1:ncol(mat)){
        if (any(is.na(mat[,j]))){
          toRemove <- c(toRemove, j)
        }
      }
      if (length(toRemove) > 0){
        mat[,-toRemove]
      }else{
        mat
      }
    }
  })
})

# remove NAs from matrix


require(plyr)
V1 <- compact(V1)
V2 <- compact(V2)
V3 <- compact(V3)
V4 <- compact(V4)
V5 <- compact(V5)

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



# --------INDEL TIMINGS -------
maxes <- maxes[!is.na(maxes)]
interval <- 100
mx <- 4000
# ---- INSERTIONS ----
ibins <- lapply(1:length(imid), function(x){
  res <- c()
  for (i in 1:(mx/interval)){
    res[i] <- sum(imid[[x]] > (i-1)*interval & imid[[x]] < i*interval)
  }
  cbind(pat=strsplit(names(imid)[x],"-")[[1]][1], as.data.frame(t(res)))
})

ibin.df <- as.data.frame(rbindlist(ibins))
colnames(ibin.df) <- c("pat",as.character(seq(0,mx,interval)[-1]))
ifreq <- apply(ibin.df[,2:ncol(ibin.df)], 2, mean)


# adjust the means for the number of patients
iadj.means <- mapply(function(bin, mean){
  # this is a calculation of how many data sets are still active, decreasing as less data is available
  adj.factor <- (length(maxes) - sum(maxes <= (bin - interval))) / length(maxes)
  print(adj.factor)
  mean / adj.factor
}, as.numeric(colnames(ibin.df[-1])), ifreq)

dbins <- lapply(1:length(dmid), function(x){
  res <- c()
  for (i in 1:(mx/interval)){
    res[i] <- sum(dmid[[x]] > (i-1)*interval & dmid[[x]] < i*interval)
  }
  cbind(pat=names(dmid)[x], as.data.frame(t(res)))
})

dbin.df <- as.data.frame(rbindlist(dbins))
colnames(dbin.df) <- c("pat",as.character(seq(0,mx,interval)[-1]))
dfreq <- apply(dbin.df[,2:ncol(dbin.df)], 2, mean)

dmaxes <- dmaxes[!is.na(dmaxes)]
# adjust the means for the number of patients
dadj.means <- mapply(function(bin, mean){
  adj.factor <- (length(maxes) - sum(maxes <= (bin - interval))) / length(maxes)
  print(adj.factor)
  mean / adj.factor
}, as.numeric(colnames(dbin.df)[-1]), dfreq)



#cairo_pdf("~/vindels/Figures/within-host/finalized/ins-timings.pdf",height=8, width=12)
par(xpd=NA, mar=c(0,6,6.5,1), mfrow=c(2,1))
barplot(iadj.means, col="dodgerblue", space=0, xaxt = "n",
        ylab="      Average Insertion Count \n Per Patient",
        cex.lab=1.3,
        cex.axis=1.2,
        cex.main=1.7,
        las=1,
        ylim=c(0,6))
        #ylim=c(0,20))
#axis(1, seq(0,15), labels=F, tick=T, line=0.5)
#text(0:15,rep(-0.7,16), labels=seq(0,7500,500), srt=25, cex=1.1)
#title(xlab="Days After Estimated Start of Infection \n(Branch Midpoints)", line=5, cex.lab=1.4)
# ---- DELETIONS ----

#cairo_pdf("~/vindels/Figures/within-host/finalized/del-timings.pdf",height=8, width=12)
par(xpd=NA, mar=c(6.5,6,0,1))
barplot(dadj.means, col="red", space=0, xaxt = "n",
        ylab="Average Deletion Count  \n Per Patient",
        #main="Deletion Timings",
        cex.lab=1.3,
        cex.axis=1.1,
        cex.main=1.7,
        las=1, ylim=c(6,0))
axis(1, seq(0,40,2), labels=F, tick=T, line=01.0)
text(seq(0,40,2),rep(7.5,21), labels=seq(0,4000,200), srt=25, cex=1.2)
title(xlab="Days After Estimated Start of Infection \n(Branch Midpoints)", line=4.75, cex.lab=1.4)
#dev.off()

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


