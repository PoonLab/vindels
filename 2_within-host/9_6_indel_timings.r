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
  
  res <- unname(sapply(iCSV$header, findAncestor, tree=tre)) 

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

iint <- iint[-toRemove]
itip <- itip[-toRemove]
dint <- dint[-toRemove]
dtip <- dtip[-toRemove]

all.data <- list(itip,iint,dtip,dint)

rm(iint)
rm(itip)
rm(dtip)
rm(dint)
rm(iCSV)
rm(dCSV)

# ---- DATA FLATTENING ----

# COMPARISON OF INSERTION RATES INTERIOR VS TIP
patnames <- unname(sapply(names(all.data[[1]]), function(x){strsplit(x,"-")[[1]][1]}))
all.counts <- list()
all.times <- list()
all.lens <- list()
all.mid <- list(list(),list())
sizes <- list()
for (a in 1:4){
  all.counts[[a]] <- list()
  all.times[[a]] <- list()
  all.lens[[a]] <- list()
  sizes[[a]] <- double()
  for (b in 1:5){
    all.counts[[a]][[b]] <- list()
    all.times[[a]][[b]] <- list()
    all.lens[[a]][[b]] <- double()
  }
}


names <- c()
for (i in 1:4){
  for (j in 1:length(unique(patnames))){
    idx <- which(patnames == unique(patnames)[j])
    if (i == 1){
      names[j] <- unique(patnames)[j] 
    }
    print(j)
    for (k in 1:5){
      mat <- sapply(idx, function(df){
        x <- all.data[[i]][[df]]
        x[x$vloop == k,"count"]
      })
      all.counts[[i]][[k]][[j]] <- mat
      
      mat2 <- sapply(idx, function(df){
        x <- all.data[[i]][[df]]
        x[x$vloop == k, "length"]
      })
      all.times[[i]][[k]][[j]] <- mat2
      
      if (k == 1){
        sizes[[i]] <- c(sizes[[i]],nrow(mat))
        #print(sizes[[i]])
      }
    }
  }
}

all.counts <- lapply(1:4, function(x){
  lapply(1:5, function(y){
    do.call(rbind, all.counts[[x]][[y]])
  })
})
all.times <- lapply(1:4, function(x){
  lapply(1:5, function(y){
    do.call(rbind, all.times[[x]][[y]])
  })
})

# ---- DATA PADDING APPROACH (NEW) ----

# COMPARISON OF INSERTION RATES INTERIOR VS TIP
patnames <- unname(sapply(names(all.data[[1]]), function(x){strsplit(x,"-")[[1]][1]}))
all.counts <- list()
all.times <- list()
all.lens <- list()
all.mid <- list(list(),list())
sizes <- list()
for (a in 1:4){
  all.counts[[a]] <- list()
  all.times[[a]] <- list()
  all.lens[[a]] <- list()
  sizes[[a]] <- double()
  # for (b in 1:5){
  #   all.counts[[a]][[b]] <- list()
  #   all.times[[a]][[b]] <- list()
  #   all.lens[[a]][[b]] <- double()
  # }
}

# --- Find maximum number of branch lengths among patients ---- 
idx <- seq(1,4201,200)
max <- c()
for (a in 1:4){
  #print(a)
  max[a] <- max(sapply(all.data[[a]][idx], nrow) / 5)
  # v3 <- sapply(all.data[[1]][idx], function(df){sum(df$vloop==3)})  # to test
}
max <- max(max)

for(x in 1:4){
  for (v in 1:5){
    counts <- array(data=Inf, dim=c(200, max, length(unique(patnames))))
    times <- array(data=Inf, dim=c(200, max, length(unique(patnames))))
    
    s <- c()
    for (y in 1:length(unique(patnames))){
      range <- idx[y]:(idx[y] + 199)
      temp <- t(sapply(all.data[[x]][range], function(df) df[df$vloop==v,3]))
      counts[1:200, 1:ncol(temp),y] <- temp 
      s[y] <- ncol(temp)
    }
    if (v == 1){
      sizes[[x]] <- s
    }
    all.counts[[x]][[v]] <- counts
    all.times[[x]][[v]] <- times
  }
}


# --- For handling completely different data sizes and formats ---- 

sizes <- list(double(),double(),double(),double())
counts <- lapply(1:4, function(x){
  lapply(1:5, function(y){
    sub <- all.counts[[x]][[y]]
    all.mats <- list()
    lowest.dim <- 10000
    for (i in 1:length(sub)){
      if (class(sub[[i]]) == "list"){
        
        lowest <- min(unname(sapply(sub[[i]], function(z) sum(!is.na(z)))))
        mat <- sapply(sub[[i]], function(a){
          notna <- which(!is.na(a))
          a[notna[sample(lowest)]]
        })
        #print(nrow(mat))
        if ( y == 1){
          print(nrow(mat))
        }
        
        
        if (dim(mat)[2] < lowest.dim){
          lowest.dim <- dim(mat)[2]
        }
        all.mats[[i]] <- mat
      }else{
        if (dim(sub[[i]])[2] < lowest.dim){
          lowest.dim <- dim(sub[[i]])[2]
        }
        if ( y == 1){
          print(nrow(sub[[i]]))
        }
        all.mats[[i]] <- sub[[i]]
      }
    }
    all.mats <- lapply(all.mats, function(b){
      b[,sample(lowest.dim)]
    })
    do.call(rbind, all.mats)
  })
})
times <- lapply(1:4, function(x){
  lapply(1:5, function(y){
    sub <- all.times[[x]][[y]]
    all.mats <- list()
    lowest.dim <- 10000
    for (i in 1:length(sub)){
      if (class(sub[[i]]) == "list"){
        lowest <- min(unname(sapply(sub[[i]], function(z) sum(!is.na(z)))))
        mat <- sapply(sub[[i]], function(a){
          notna <- which(!is.na(a))
          a[notna[sample(lowest)]]
        })
        if (dim(mat)[2] < lowest.dim){
          lowest.dim <- dim(mat)[2]
        }
        all.mats[[i]] <- mat
      }else{
        if (dim(sub[[i]])[2] < lowest.dim){
          lowest.dim <- dim(sub[[i]])[2]
        }
        all.mats[[i]] <- sub[[i]]
      }
    }
    all.mats <- lapply(all.mats, function(b){
      b[,sample(lowest.dim)]
    })
    do.call(rbind, all.mats)
  })
})

all.mid <- lapply(1:2, function(x){
  unlist(all.mid[[x]])
})



# --------INDEL TIMINGS -------
maxes <- maxes[!is.na(maxes)]
interval <- 200
mx <- 4800
# ---- INSERTIONS ----

counts <- lapply(1:2, function(a){
  res <- c()
  for (i in 1:(mx/interval)){
    res[i] <- sum(all.mid[[a]] > (i-1)*interval & all.mid[[a]] < i*interval)
  }
  names(res) <- as.character(seq(interval,mx,interval))
  return(res / 400 / 26)
})


# adjust the means for the number of patients
means <- lapply(1:2, function(x){
    mapply(function(bin, mean){
    # this is a calculation of how many data sets are still active, decreasing as less data is available
    adj.factor <- sum(maxes > (bin - interval))/ length(maxes)
    print(adj.factor)
    mean / adj.factor
  }, as.numeric(names(counts[[x]])), counts[[x]])
})


#cairo_pdf("~/vindels/Figures/within-host/finalized/ins-timings.pdf",height=8, width=12)
par(xpd=NA, mar=c(0,6,6.5,1), mfrow=c(2,1))
barplot(means[[1]], col="dodgerblue", space=0, xaxt = "n",
        ylab="",
        yaxt="n",
        cex.lab=1.3,
        cex.axis=1.2,
        cex.main=1.7,
        las=1,
        ylim=c(0,4))
axis(2, 0:4, labels=c("",1:4), tick=T,cex.axis=1.2)
title(ylab=" Average Insertion Count \nPer Patient", cex.lab=1.4, line=2.2)
        #ylim=c(0,20))
#axis(1, seq(0,15), labels=F, tick=T, line=0.5)
#text(0:15,rep(-0.7,16), labels=seq(0,7500,500), srt=25, cex=1.1)
#title(xlab="Days After Estimated Start of Infection \n(Branch Midpoints)", line=5, cex.lab=1.4)
# ---- DELETIONS ----

#cairo_pdf("~/vindels/Figures/within-host/finalized/del-timings.pdf",height=8, width=12)
par(xpd=NA, mar=c(6.5,6,0,1))
barplot(means[[2]], col="red", space=0, xaxt = "n",
        ylab="",
        #main="Deletion Timings",
        cex.lab=1.3,
        cex.axis=1.1,
        cex.main=1.7,
        las=1, ylim=c(4,0))
axis(1, seq(0,25,2), labels=F, tick=T, line=0.5)
text(seq(0,40,2),rep(4.5,21), labels=seq(0,4800,400), srt=25, cex=1.2)
title(xlab="Days After Estimated Start of Infection \n(Branch Midpoints)", line=4.75, cex.lab=1.4)
title(ylab="Average Deletion Count \nPer Patient", cex.lab=1.4, line=2.2)
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


