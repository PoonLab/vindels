library("rjags")
library("data.table")
# ---- Simulate Data---- 

rate <- 0.00001
time <- round(rlnorm(50000,4.5,0.3))

lambda <- rate * time

y <- rpois(50000,lambda=lambda)

data <- data.frame(time, y)

# ---- STANDARD GLM ----
fit <- glm(y ~ 1, offset=log(time), data=df, family='poisson')
coef(fit)[[1]]

fit2 <- glm(y ~ offset(log(time)), data=df, family='poisson')
coef(fit2)[[1]]


# ---- REAL DATA -----

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

ci.rate <- function(counts, times, vlen, n=100){
  mean.rate <- est.rate(counts, times, vlen)
  
  # split by full.id and select 
  
  boot <- replicate(n, one.boot(counts, times, vlen))
  centered <- quantile(boot - mean.rate, c(0.025,0.975), na.rm=T)
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
    df <- new.data[[x]][[v]]
    df <- split(df, df$pat)
    idx <- get.pats(df)
    df <- do.call(rbind, df[idx])
    res <- ci.rate(df$count, df$length, df$vlen)
    lower[v,x] <<- res$ci[1]
    upper[v,x] <<- res$ci[2]
    res$rate
  })
  names(z) <- c('V1','V2','V3','V4','V5')
  z
})

# ---- SPLIT BY PATIENT (BOXPLOT) ---- 

new.data <- lapply(1:4, function(a){
  lapply(1:5, function(b){
    n <- new.data[[a]][[b]]
    split(n, n$pat)
  })
})

# ---- PATIENT-WISE INDEL RATES (BOXPLOT) ----
rates <- lapply(1:4, function(x){
  sapply(c(1,2,3,4,5), function(v){
    sapply(1:24, function(p){
      df <- new.data[[x]][[v]][[p]]
      est.rate(df$count, df$length, df$vlen)
    })
  })
})

# Load final data frame 
# (terminal, internal)
final <- list(cbind(rates[[1]], rates[[2]]),
              cbind(rates[[3]], rates[[4]]))
ind <- c()
ind[seq(1,10, by=2)] <- 1:5
ind[seq(2,10, by=2)] <- 6:10


# Graphical parameters
par(mfrow=c(1,2), mar=c(6,2,2,1))
lims <- list(c(20,0), c(0,20))

cols <- c('dodgerblue','red')

data <- final[[2]]
data <- data[,ind]
pos <- c(1.1,1.9,3.1,3.9, 5.1,5.9, 7.1, 7.9, 9.1,9.9)
boxplot(data, xlab="Deletion Rate", horizontal=TRUE, at=pos, cex.axis=1.3,
        cex.lab=1.6, ylim=c(20,0), yaxt='n', col=rep(cols, 5))
axis(4, at=seq(0.5, 10.5, by=2), labels=F)
axis(4, at=seq(1.5,9.5,by=2),tick=F, line=-0.3,
     labels=c("V1","V2","V3","V4","V5"), las=1, cex.axis=1.5)
data <- final[[1]]
data <- data[,ind]
xlabel <- expression(paste("      Insertion Rate\n(events/nt/year x 10"^"-3"*")"))
boxplot(data, xlab="", at=pos,  cex.axis=1.3,
        horizontal=TRUE, col=rep(c('dodgerblue','red')), cex.lab=1.3, ylim=c(0,20), yaxt='n')
axis(2, at=seq(0.5, 10.5, by=2), labels=F)
title(xlab=xlabel, cex.lab=1.6, line=4.3)
legend(x=10, y=6.5, legend=c("Terminal Branches", "Internal Branches"), 
       pch=c(22,22), pt.bg=cols, pt.cex=3, cex=1.4)

# --- PLOT PREPARATION ---- 

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

## Manual fixing
irates[8,4] <- 0
drates[8,4] <- 0


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
  scale_y_continuous(lim=c(0,11), expand=c(0,0)) + 
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
        legend.position=c(0.70,0.59),
        legend.text=element_text(size=18), 
        legend.background=element_rect(colour="black"),
        legend.title=element_text(size=20),
        legend.spacing.y = unit(2, "mm")
  ) + geom_text(aes(y=c(1),x=c(3.25)),label="N/A", size=7)
#iplot
dplot <- ggplot() + 
  geom_bar(aes(x=vloop, y=rate, fill=id), data=drates, stat='identity', position="dodge") + 
  coord_flip() + scale_fill_manual(values=c( "dodgerblue", "red"))+
  geom_errorbar(aes(x=drates$vloop, fill=drates$id, ymax = drates$upper, ymin = drates$lower),
                width = 0.25, size=0.8,
                position = position_dodge(0.9)) +
  labs(x="Variable Loop", 
       y="Deletion Rate") +
  scale_y_reverse(lim=c(11,0), expand=c(0,0)) + 
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
        legend.position="none") #+ geom_text(aes(y=1,x= 3.25),label="N/A", size=7)
multiplot(dplot,iplot, cols=2)

# new.data <- lapply(1:2, function(a){
#   lapply(1:5, function(b){
#     n <- new.data[[a]][[b]]
#     split(n, n$pat)
#   })
# })
# 
# new.data <- 


# ---- Simulate Data RSTAN ---- 

# Data import
data.stan <- list(N=50000,
                  counts=data$y,
                  branches=data$time)

# Stan modeling 
start <- proc.time()
stan.fit <- stan("~/vindels/2_within-host/stan_modeling/rate-new.stan",
                 data= data.stan, 
                 chains=1,
                 iter=1000,
                 control=list(adapt_delta=0.90))
end <- proc.time() - start


final.data <- list(itip, iint, dtip, dint)

final.data <- lapply(final.data, function(x){
  
  # Filter for length 0 
  idx <- x$length == 0
  x[idx,'length'] <- 0.001
  
  # Extract patient strings
  res <- t(sapply(x$pat, function(str){
    
    parsed <- strsplit(str, "_")[[1]]
    # if(grepl("-b_", str)){
    #   parsed[2] = paste(as.numeric(parsed[2]) + 200)
    # }
    # parsed[1] <- strsplit(parsed[1], "-")[[1]][1]
    parsed[1] <- substr(parsed[1], 1, nchar(parsed[1])-2)
    parsed
  }))
  x[, 'pat'] <- res[,1]
  x[, 'rep'] <- res[,2]
  
  split(x, x$vloop)
})

# --- Preprocessing of Data ---

final.data <- lapply(1:4, function(x){
  lvl1 <- final.data[[x]]
  
  lapply(1:5, function(vloop){
    vloop <- final.data[[x]][[vloop]]
    split(vloop , vloop$rep)
  })
})

# CHECK FOR LENGTH 0 ENTRIES
test <- lapply(1:4, function(x){
  lapply(1:5, function(v){
    lapply(1:200, function(rep){
      final.data[[x]][[v]][[rep]]$length
    })
  })
})
sum(unlist(x) == 0)


# REMOVE PATIENTS WITH NO INDELS 
# (getting data for Roux-Cil)
x <- final.data[[1]]
x <- x[x$vloop==1, ]
treeno <- sapply(x$pat, function(a)strsplit(a,"_")[[1]][2])
names <-sapply(x$pat, function(a)strsplit(a,"-")[[1]][1])
x$name <- names 
x$treeno <- treeno
y <- split(x, x$name)
counts <- sapply(y, function(a) sum(a$count))
zeros <- names(which(counts==0))
x = x[!x$name %in% zeros,]
final <- x[,-c(1,2,4,6)]
colnames(final)[2] <- 'times'

y <- split(final, final$name)

# Create lists of matrices 
clist<- lapply(y, function(a){matrix(a$count, ncol=200)})
tlist <- lapply(y, function(a){matrix(a$times, ncol=200)})

# Concatenate matrices together
counts <- do.call(rbind, clist)
times <- do.call(rbind, tlist)




data <- final.data[[1]][[1]]
options(mc.cores = parallel::detectCores()-2)
# Stan modeling 
start <- proc.time()

vmean <- list()
vmed <- list()
vneff <- list()

for (v in 1:5){
  data <- final.data[[1]][[v]]
  mean <- c()
  median <- c()
  neff <- c()
  
  for (i in 1:200){
    # Data import
    df <- data[[i]]
    stan.df <- list(N=nrow(df),
               counts=df$count,
               branches=df$length)
    
    stan.fit <- stan("~/vindels/2_within-host/stan_modeling/rate-new.stan",
                     data=stan.df, 
                     chains=1,
                     iter=10000,
                     control=list(adapt_delta=0.90))
    sum <- summary(stan.fit)$summary
    mean[i] <- exp(sum[1,1]) * (365 / median(df$vlen))
    median[i] <- exp(sum[1,6]) * (365 / median(df$vlen))
    neff[i] <- sum[1,'n_eff']
    print(paste0("Finished: ", i))
  }
  vmean[[v]] <- mean
  vmed[[v]] <- median
  vneff[[v]] <- neff
}


end <- proc.time() - start

lens <- lapply(1:5, function(v){
  sapply(1:200, function(x){
    median(final.data[[1]][[v]][[x]]$vlen)
  })
})

final <- lapply(1:5, function(v){
  sapply(1:200, function(x){
    1000* exp(vmean[[v]][x]) * (365 / median(lens[[v]][x]))
  })
})







# DEPRECATED ----------------
# ----- REAL DATA -----
# ------------------------

# ---- STANDARD GLM ----
df <- final.data[[1]][[1]]
fit <- glm(count ~ 1, offset=log(length), data=df, family='poisson')
coef(fit)[[1]]

fit2 <- glm(count ~ offset(log(length)), data=df, family='poisson')
coef(fit2)[[1]]



# ---- REAL DATA RSTAN ITERATED ----
# POOLING ACROSS ALL PATIENTS 
# NOT MEANINGFUL
data <- final.data[[1]][[1]]


options(mc.cores = parallel::detectCores()-2)
# Stan modeling 
start <- proc.time()

mean <- c()
median <- c()
neff <- c()


for (v in 1:5){
  # Load data 
  data <- final.data[[1]][[v]]
  data$rep <- as.numeric(data$rep)
  
  # Dimensions of matrices 
  npat = length(unique(data$pat))
  ntree = 50 # length(unique(data$rep))
  
  # Count number of branches per patient 
  sizes <- c(table(data[data$rep=="1",'pat']))
  
  # Split by patient 
  byPat <- split(data, data$pat)
  byPat <- lapply(byPat, function(x){
    x[order(x$rep),]
  })
  
  # Search for empty patients 
  pos = 1 
  sums <- c()
  for(i in 1:npat){
    start = pos
    end = pos + sizes[i]-1
    
    sums[i] <- sum(byPat[[i]]$count)
    
    pos = pos + sizes[i]
  }
  
  # Modify sizes to remove the empty patients
  sizes <- sizes[-which(sums==0)]
  npat <- npat - sum(sums==0)
  byPat <- byPat[-which(sums==0)]

  # Initialize matrices
  all.counts <- matrix(nrow=sum(sizes),
                       ncol=ntree)
  all.times <- matrix(nrow=sum(sizes),
                      ncol=ntree)

  # Load matrices
  pos = 1
  for(i in 1:npat){
    start = pos
    end = pos + sizes[i]-1
    print(paste0(start," ",end))
    all.counts[start:end,] = byPat[[i]]$count[sample(ntree)]
    all.times[start:end,] = log(byPat[[i]]$length[sample(ntree)])

    pos = pos + sizes[i]
  }

  stan.df <- list(npat=npat,
                  ntree=ntree,
                  sizes=sizes,
                  counts=all.counts,
                  l_branches=all.times)

  stan.fit <- stan("~/vindels/2_within-host/stan_modeling/rate-perpat-dc.stan",
                   data=stan.df,
                   chains=1,
                   iter=10000,
                   #control=list(adapt_delta=0.90)
                   )
  end <- proc.time() - start

  # Save summary
  sum <- summary(stan.fit)$summary

  # Load results
  mean[i] <- exp(sum[1,1]) * (365 / median(df$vlen))
  median[i] <- exp(sum[1,6]) * (365 / median(df$vlen))
  neff[i] <- sum[1,'n_eff']

  print(paste0("Finished: ", i))

}

