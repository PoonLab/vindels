require(ape)
require(stringr)
require(phangorn)
require(data.table)
require(bbmle)
source("~/vindels/2_within-host/utils.r")

# INSERTION PARSING ----------
path <- "~/work/"
ifolder <- Sys.glob(paste0(path,"9Indels/new-final/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"9Indels/new-final/del/*.tsv"))
sep <- "\t"
trefolder <- paste0(path,"7SampleTrees/new-final/prelim/")

# CASE: Removed patients 49641 and 56549 because they are SUBTYPE B
# reg <- "49641|56549|108869"
# 
# ifolder <- ifolder[!grepl(reg,ifolder)]
# dfolder <- dfolder[!grepl(reg,dfolder)]

# first20 <- "_1{0,1}[0-9]_|_20_"
# 
# ifolder <- ifolder[grepl(first20, ifolder)]
# dfolder <- dfolder[grepl(first20, dfolder)]

tally <- function(infolder){
  name <- basename(infolder)
  name <- gsub("-.+","",name)
  return (table(name))
}

all.ins <- c()
all.del <- c()

maxes <- c()
max_names <- c()

iint <- list()
itip <- list()
dint <- list()
dtip <- list()

finalCols <- c('full.id',"pat", "rep", "vloop", "vlen", "count", "length", "rtt.mid")


for (file in 1:length(ifolder)){
  print(file)
  filename <- basename(ifolder[file])
  full.id <- gsub("_\\d+\\.tsv$","",filename)

  iCSV <- read.csv(ifolder[file], stringsAsFactors = F, sep=sep)
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F, sep=sep)
  
  iCSV$count <- unname(sapply(iCSV$indel, csvcount, delim=":"))
  dCSV$count <- unname(sapply(dCSV$indel, csvcount, delim=":"))
  
  iCSV$full.id <- rep(full.id, nrow(iCSV))
  dCSV$full.id <- rep(full.id, nrow(dCSV))
  
  iCSV$header <- unname(mapply(labels, iCSV$header, iCSV$full.id))
  dCSV$header <- unname(mapply(labels, dCSV$header, dCSV$full.id))
  
  # reads in the tree
  treename <- strsplit(filename, "\\.")[[1]][1]
  tre <- read.tree(paste0(paste0(trefolder,treename,".tree.sample")))
  
  res <- unname(sapply(iCSV$header, findAncestor, tree=tre)) 
  
  rttlens <- node.depth.edgelength(tre)
  
  # midpoint = (rtt length of tip) + (rtt length of ancestor) / 2
  midpoints <- (rttlens[res] + rttlens[tre$edge[match(res, tre$edge[,2]),1]]) / 2
  
  iCSV$length <- tre$edge.length[match(res, tre$edge[,2])]   # midpoints (if we want full depth)
  dCSV$length <- iCSV$length
  
  iCSV$rtt.mid <- midpoints
  dCSV$rtt.mid <- midpoints
  
  res <- as.data.frame(t(unname(sapply(iCSV$full.id, extractPat, removeLetter=T))))
  colnames(res) <- c("pat", "rep")
  
  iCSV <- cbind(iCSV, res)
  dCSV <- cbind(dCSV, res)
  
  # use the full tree lengths as the maximum cutoff
  maxes[file] <- max(rttlens,na.rm=T)
  max_names[file] <- full.id
  
  # ----- TIP + INTERIOR INDEL COUNTS ---
  # Cumulative data frame split by interior vs tip
  id <- full.id #gsub("[ab]_","",full.id)   # 16362-100
  
  # this regexp matches TIP SEQUENCES (NOT CONTAINING LEFT AND RIGHT BRACKETS)
  tips <-  which(grepl("^[^\\(\\):\n]+$", iCSV$header))
  nodes <- which(!grepl("^[^\\(\\):\n]+$", iCSV$header))
  
  iCSV <- iCSV[,-c(1,2,5,6)]
  dCSV <- dCSV[,-c(1,2,5,6)]
  
  iCSV <- iCSV[,finalCols]
  dCSV <- dCSV[,finalCols]

  # if (is.null(iint[[id]])){
  iint[[id]] <- iCSV[nodes,]
  itip[[id]]  <- iCSV[tips,]
  dint[[id]]  <- dCSV[nodes,]
  dtip[[id]]  <- dCSV[tips,]
  # }else{
  #   iint[[id]] <- rbind(iint[[id]], iCSV[nodes,])
  #   itip[[id]]  <- rbind(itip[[id]], iCSV[tips,])
  #   dint[[id]]  <- rbind(dint[[id]] , dCSV[nodes,])
  #   dtip[[id]]  <- rbind(dtip[[id]], dCSV[tips,])
  # }
}
setwd("~/vindels/2_within-host/")
names(maxes) <- max_names
## ========================================
## RSESSION CHECKPOINT : 9_6_all_patients
## ========================================


patnames <- unname(sapply(names(iint), function(x){strsplit(x, "-")[[1]][1]}))
pat.idx <- table(sapply(ifolder, function(x){
   strsplit(basename(x), "-")[[1]][1]
}))

numpats = length(unique(patnames))

all.data <- list(itip, iint, dtip, dint)

rm(iint)
rm(itip)
rm(dtip)
rm(dint)
rm(iCSV)
rm(dCSV)
rm(ifolder)
rm(dfolder)

## ------ SAMPLE ------
# Permanently choose a subset of 100 trees to examine 
require(data.table)
set.seed(0)
idx <- as.vector(sapply(1:numpats, function(i){
  sample(((i-1)*400+1):(i*400), 100, replace=FALSE)
}))

## --- DATA SAMPLING -----
sub.data <- lapply(all.data, function(x) as.data.frame(rbindlist(x[idx])))

# SUBDATA IS USED FOR BOTH SUBSEQUENT ANALYSES

## INDEL RATES ANALYSIS
# -----------------------------

new.data <- lapply(sub.data, function(x) {
  tmp <- x[x$length > 0,]
  split(tmp, tmp$vloop)
})
sizes <- lapply(new.data, function(x){
  sapply(split(x[[2]], x[[2]]$pat), nrow)
})


lower <- matrix(nrow=5, ncol=4)
upper <- matrix(nrow=5, ncol=4)

est.rate <- function(counts, times, vlen, bs=FALSE){
  r <- NA
  tryCatch({
    fit <- glm(counts ~ 1, offset=log(times), family='poisson')
    # if (bs){
    #   ln = sample(vlen,1)
    # }else{
    #   ln = median(vlen)
    # }
    
    ln = median(vlen)
    
    if(ln > 250){
      ln <- 127
    }
    
    r <- exp(coef(fit)[[1]]) / ln * 365 * 1000
  },warning=function(cond){
  })
  return(r)
}

one.boot <- function(counts, times, vlen){
  n <- length(counts)
  sam <- sample(1:n, n, replace=T)
  est.rate(counts[sam], times[sam], vlen[sam], T)
}

ci.rate <- function(vdata){
  mean.rate <- est.rate(vdata$count, vdata$length, vdata$vlen)
  
  # split by full.id and select 
  df <- split(vdata, vdata$full.id)
  
  bs.vals <- sapply(1:100, function(bs){
    idx <- bs + seq(0,(numpats-1)*100,100)
    count <- unlist(sapply(idx, function(i) df[[i]][,"count"]))
    time <- unlist(sapply(idx, function(i) df[[i]][,"length"]))
    vlen <- unlist(sapply(idx, function(i) df[[i]][,"vlen"]))
    est.rate(count,time,vlen)
  })
  
  centered <- quantile(bs.vals - mean.rate, c(0.025,0.975), na.rm=T)
  ci <- c(mean.rate - centered[2], mean.rate - centered[1])
  return(list(rate=mean.rate, ci=ci))
}

get.pats <- function(df){
  counts <- sapply(df, function(x) sum(x$count))
  which(counts != 0)
}

# ---- STANDARD INDEL RATE CALCULATION ----
rates <- sapply(1:4, function(x){
  print(paste0("main ",x))
  z <- sapply(c(1,2,3,4,5), function(v){
    print(paste0("vloop ",v))
    # df <- split(df, df$pat)
    # idx <- get.pats(df)
    # df <- do.call(rbind, df[idx])
    res <- ci.rate(new.data[[x]][[v]])
    lower[v,x] <<- res$ci[1]
    upper[v,x] <<- res$ci[2]
    res$rate
  })
  names(z) <- c('V1','V2','V3','V4','V5')
  z
})


## --- INDEL RATES HORIZONTAL BAR PLOT ---- 
# PLOT PREP
cols <- c('vloop', 'id', 'rate','lower','upper')
f <- 365*1000
irates <- reshape2::melt(rates[,1:2] )
drates <- reshape2::melt(rates[,3:4] )
irates <- cbind(irates, reshape2::melt(lower[,1:2])[,3] , reshape2::melt(upper[,1:2])[,3] )
drates <- cbind(drates, reshape2::melt(lower[,3:4])[,3] , reshape2::melt(upper[,3:4])[,3] )

colnames(irates) <- cols
colnames(drates) <- cols

irates$id <- as.factor(irates$id)
drates$id <- as.factor(drates$id)

levels(irates$id) <- c('terminal', 'internal')
levels(drates$id) <- c('terminal', 'internal')

## Manual fixing for vals bordering 0
irates[8,3:5] <- 0
drates[8,3:5] <- 0

# Plotting 
require(Rmisc)
require(ggplot2)
par(mar=c(3,3,3,2))
iplot <- ggplot() + 
  geom_bar(aes(vloop, rate, fill=id), data=irates, stat='identity', position="dodge") + 
  scale_fill_manual(values=c( "dodgerblue", "red"), name="Type", labels=c("Terminal Branches","Internal Branches")) +
  coord_flip() + geom_errorbar(aes(x=irates$vloop, fill=irates$id, ymax = irates$upper, ymin = irates$lower),
                               width = 0.25, size=0.8,
                               position = position_dodge(0.9)) +
  scale_y_continuous(lim=c(0,9), expand=c(0,0)) + 
  scale_x_discrete(limits = rev(levels(irates$vloop))) +
  labs(#x="Variable Loop", 
    y=substitute(paste(p1, 10^-3,p2), list(p1="         Insertion Rate\n    (events/nt/year x ", p2=")"))) + 
  #title="Indel Rates",
  #color="Subset") +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),#element_line(color="black",size=0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 4, r = 10, b = 5, l = 0, unit = "mm"),
        panel.border=element_rect(fill=NA, size=1),
        #axis.line.x = element_line(colour = "black"), 
        #axis.line.y.left=element_line(colour="black"),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=22,margin=margin(t = 22, r = 3, b = 0, l = 22)),
        axis.text.y = element_text(size=20, colour="black", margin=margin(t = 0, r = 10, b = 2, l = 0)),
        axis.text.x=element_text(size=20, colour="black",margin=margin(t = 0, r = 20, b = 0, l = 0)),
        #plot.title = element_text(size=22, hjust = 0.5),
        legend.position=c(0.65,0.59),
        legend.text=element_text(size=18), 
        legend.background=element_rect(colour="black"),
        legend.title=element_text(size=20),
        legend.spacing.y = unit(2, "mm")
  ) + geom_text(aes(y=c(0.7),x=c(3.25)),label="N/A", size=7)
#iplot
dplot <- ggplot() + 
  geom_bar(aes(x=vloop, y=rate, fill=id), data=drates, stat='identity', position="dodge") + 
  coord_flip() + scale_fill_manual(values=c( "dodgerblue", "red"))+
  geom_errorbar(aes(x=drates$vloop, fill=drates$id, ymax = drates$upper, ymin = drates$lower),
                width = 0.25, size=0.8,
                position = position_dodge(0.9)) +
  labs(x="Variable Loop", 
       y="Deletion Rate") +
  scale_y_reverse(lim=c(9,0), expand=c(0,0)) + 
  scale_x_discrete(limits = rev(levels(drates$vloop))) +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),#element_line(color="black",size=0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        panel.border=element_rect(fill=NA, size=1),
        #axis.line.x=element_line(colour = "black"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 4, r = 5, b = 13.5, l = 12, unit = "mm"),
        #axis.line.y = element_line(colour = "black"), 
        axis.text.y = element_blank(),
        axis.title.y=element_text(size=22),
        axis.title.x=element_text(size=22,margin=margin(t = 7, r = 3, b = 0, l = 12)),
        axis.text = element_text(size=20, colour="black"),
        plot.title = element_text(size=28, hjust = 0.5),
        legend.position="none") + geom_text(aes(y=0.7,x=3.25),label="N/A", size=7)
multiplot(dplot,iplot, cols=2)



## INDEL TIMINGS ANALYSIS
# -----------------------------
# Samples one bootstrap from the all.data list structure
bs.timing <- function(){
  bs.idx <- sample(idx, replace=T)
  #print(paste(bs.idx, sep=' '))
  lapply(all.data, function(x) {
    as.data.frame(rbindlist(x[bs.idx]))   # arbitrarily choose to analyze the first 10 patients
  })
}


max.times <- maxes[idx]
interval <- 100
mx <- 1000
bins <- mx/interval
bin.names <- as.character(seq(interval,mx,interval))


timing.calc <- function(data, max.times=max.times){
  
  # find the mean counts across bins and across patients 
  counts <- lapply(0:1, function(a){
    x1 = data[[(a*2+1)]]
    x2 = data[[(a*2+2)]]
    rtt.mid <- c(x1[x1$count > 0,'rtt.mid'], x2[x2$count > 0,'rtt.mid'])
    #print(max(rtt.mid))
    res <- sapply(1:bins, function(b){
      sum(rtt.mid > (b-1)*interval & rtt.mid < b*interval)
    })
    names(res) <- bin.names
    return(as.data.frame( t(res) / 100 / 2 / numpats))
  })
  
  # adjust the means for the number of patients
  adj.factor <- sapply(as.numeric(names(counts[[1]])), function(bin) mean(max.times > (bin - interval)))
  lapply(counts, function(c) c / adj.factor)
}

# MOST UPDATE TO DATE METHOD 
conf <- function(data, max.times=max.times){
  xsample <- timing.calc(data, max.times)
  
  counts <- lapply(0:1, function(i){
    x1 = data[[(i*2+1)]]
    x2 = data[[(i*2+2)]]
    
    # different per-patient approach ** (splitting)
    df1 <- split(x1, x1$full.id)
    df2 <- split(x2, x2$full.id)
  
    sapply(1:100, function(bs){
      idx <- bs + seq(0,(numpats-1)*100,100)
      rtt.mid <- unlist(lapply(idx, function(i){
        c(df1[[i]][df1[[i]]$count>0,'rtt.mid'], df2[[i]][df2[[i]]$count>0,'rtt.mid'])
      }))
      res <- sapply(1:bins, function(b){
        sum(rtt.mid > (b-1)*interval & rtt.mid < b*interval)
      })
      names(res) <- bin.names
      return( res / 2 / numpats)
    })
  })

  adj.factor <- sapply(as.numeric(rownames(counts[[1]])), function(bin) mean(max.times > (bin - interval)))
  bs.vals <- lapply(1:2, function(i){
    df <- counts[[i]] / adj.factor
    res <- as.data.frame(sapply(1:nrow(df), function(row){
      m <- xsample[[i]][,row]
      q <- quantile(df[row,] - m, c(0.025,0.975))
      c(m - q[[2]], m - q[[1]])
    }))
    colnames(res) <- bin.names
    res
  })
  lapply(1:2, function(i){rbind(xsample[[i]], bs.vals[[i]])})
  
}


timing.conf <- function(n=100){
  # Create n bootstrap replicates and 
  res <- t(replicate(n, timing.calc(bs.timing(), max.times)))
  bs.vals <- apply(res, 2, function(x) as.data.frame(rbindlist(x)))
  
  lapply(1:2, function(y){
    tmp <- as.data.frame(sapply(1:ncol(bs.vals[[y]]), function(z){
      m <- xsample[[y]][1,z]
      bs.quant <- quantile(bs.vals[[y]][,z] - m, c(0.025,0.975))
      ci <- c(m - bs.quant[[2]], m - bs.quant[[1]])
      #print(ci)
      ci
    }))
    names(tmp) <- bin.names
    tmp
  })
}

# New method (patient-wise bootstrapping)
final <- conf(sub.data, max.times)


# ---- DATA  EXPORT ----

# data <- sub.vdata[[1]][[1]]
# patsize <- sizes[[1]]
# dump("data", "~/rate-modeling/v1.data")
# dump("patsize", "~/rate-modeling/patsize.data")
# rm(patsize )
# rm(data)


# ----- INDEL TIMINGS PLOT -----
ymax <- 2.1
par(xpd=NA, mar=c(0,6,6,1), mfrow=c(2,1))
data <- t(final[[1]])
barplot(data[,1], col="dodgerblue", space=0, xaxt = "n",
        ylab="",
        yaxt="n",
        cex.lab=1.3,
        cex.axis=1.2,
        cex.main=1.7,
        las=1,
        ylim=c(0,ymax))

arrows(1:nrow(data)-0.5, data[,2], 1:nrow(data)-0.5, data[,3], length=0.05, angle=90, code=3,lwd=1.5)
axis(2, 0:ymax, labels=c(0:ymax), tick=T,cex.axis=1.2, las=1)
title(ylab="Insertions", cex.lab=1.4, line=2.2)
par(xpd=F)
abline(h=0,lwd=1)

# ---- DELETIONS ----

#cairo_pdf("~/vindels/Figures/within-host/finalized/del-timings.pdf",height=8, width=12)
par(xpd=NA, mar=c(6.5,6,0,1))
data <- t(final[[2]])
barplot(data[,1], col="red", space=0, xaxt = "n",
        ylab="",
        yaxt='n',
        #main="Deletion Timings",
        cex.lab=1.3,
        cex.axis=1.1,
        cex.main=1.7,
        las=1, ylim=c(ymax,0))
arrows(1:nrow(data)-0.5, data[,2], 1:nrow(data)-0.5, data[,3], length=0.05, angle=90, code=3,lwd=1.5)
axis(1, seq(0,bins,2), labels=F, tick=T, line=0.5)
axis(2, 0:ymax, labels=c("",1:ymax), tick=T,cex.axis=1.2, las=1)
text(seq(0,bins,2),rep(2.35,bins/2), labels=seq(0,mx,200), srt=0, cex=1.2)
title(xlab="Days After Estimated Start of Infection \n(Branch Midpoints)", line=4.75, cex.lab=1.4)
title(ylab="Deletions", cex.lab=1.4, line=2.2)
mtext("\t\t\t\t\tAverage Count Per Patient", side=2, line=4,cex=1.4)
par(xpd=F)
abline(h=0,lwd=1)
#dev.off()





## ---- SUPPLEMENTARY; NOT MAIN ------

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



# ---- DATA FLATTENING ----

# COMPARISON OF INSERTION RATES INTERIOR VS TIP
#patnames <- unname(sapply(names(all.data[[1]]), function(x){strsplit(x,"-")[[1]][1]}))
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




