# this script is solely dedicated to statistical analysis of indel rates
# First part simulates data to run tests in RStan
# Second part relies on existing data being loaded from 9_6 

# stan model for indel rates 
library(rstan)
setwd("~/vindels/2_within-host/")


# ---- Simulated Data ----
num.pat <- 20
num.trees <- 200

tru.rate <- 0.003
tru.sd <- 0.1
pat.rate <- rlnorm(num.pat, mean=log(tru.rate), sd=tru.sd)
pat.sd <- 0.1



hist(pat.rate, breaks=5, col="grey")
vlens <- c(75, 126, 105,  90 , 33)
lns <- sample(70:200, num.pat, replace=T)   # round(rlnorm(50000,4.5,0.3))
counts <- matrix(nrow=sum(lns), ncol=num.trees)
lens <- matrix(nrow=sum(lns), ncol=num.trees)
pos <- 1
toCheck <- c()

for (i in 1:num.pat){
  
  idx <- pos:(pos+lns[i]-1)
  print(range(idx))
  tre.rate <- rlnorm(num.trees, mean=log(pat.rate[i]), sd=pat.sd)
  toCheck[((i-1)*30+1):(i*30)] <- tre.rate
  
  for (j in 1:num.trees){
    num <- diff(range(idx))+1
    len.temp <- rlnorm(num,4.5,0.3)   #rexp(diff(range(idx))+1, 0.1)
    lens[idx, j] <- len.temp
    c <- rpois(num, lambda=tre.rate[j]*len.temp)
    if(any(is.na(c))){
      print(tre.rate[j])
      #print(len.temp)
    }
    counts[idx,j] <- c
  }
  pos <- pos + lns[i]
}

rstan_options(auto_write=T)
options(mc.cores = parallel::detectCores()-1)

# Data import
data.stan <- list(npat= num.pat,
                  ntree = num.trees,
                  sizes = lns,
                  counts = counts,
                  branches = lens)

# Stan modeling 
stan.fit <- stan("~/vindels/2_within-host/stan_modeling/rate-perpat-dc.stan",
                 data= data.stan, 
                 chains=1,
                 iter=10000)
# ---- Figure for simulated data ----
par(mar=c(5,6.5,3,2))
hist(w$sub_rate, freq=F, 
     col="dodgerblue",
     xlab="Posterior of Indel Rate",
     las=1,
     cex.axis=1.3,
     cex.lab=1.7, 
     main="",
     ylab="")
abline(h=0.2000004, col="red", lwd=3)
abline(v=0.03, lwd=3)
title(ylab="Density", line=4.2, cex.lab=1.7)



# ---- Real Data ----


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



type <- c("tip","node")
rstan_options(auto_write=T)
setwd("~/vindels/2_within-host/")
set.seed(3124)
options(mc.cores = parallel::detectCores()-8)
irates <- data.frame(stringsAsFactors = F)
drates <- data.frame(stringsAsFactors =F)
ifits <- list()
dfits <- list()
loops <- c(1,2,3,4,5)
for (t in 1:2){ 
  for (i in 1:5){
    cmat <- all.counts[[t]][[loops[i]]]
    tmat <- all.times[[t]][[loops[i]]]
    
    s <- sizes[[t]]
    
    num.pat <- length(s)
    num.tree <- ncol(cmat)

      
    # Data import
    data.stan <- list(npat=num.pat,
                      ntree = num.tree,
                      sizes= s,
                      counts = cmat,
                      branches = tmat)

    # Stan modeling 
    start <- proc.time()
    stan.fit <- stan("stan_modeling/rate-perpat-dc.stan",
                     data= data.stan, 
                     chains=1,
                     iter=100000)
    end <- proc.time() - start
    # export results to a data frame 
    irates <- rbind(irates, data.frame(rate=summary(stan.fit)$summary[1,6],
                                       vloop=paste0("V",as.character(loops[i])),
                                       id=type[t],
                                       lower=summary(stan.fit)$summary[1,4],
                                       upper=summary(stan.fit)$summary[1,8]
    ))     
    idx <- (t-1) * 4 + i
    ifits[[idx]] <- stan.fit
    
    # ----- Deletions ---- 

    cmat <- all.counts[[t+2]][[loops[i]]]
    tmat <- all.times[[t+2]][[loops[i]]]
    
    #cmat <- cmat[,1:100]
    #tmat <- tmat[,1:100]
    
    s <- sizes[[t+2]]
    
    num.pat <- length(s)
    num.tree <- ncol(cmat)
    
    
    # Data import
    data.stan <- list(npat=num.pat,
                      ntree = num.tree,
                      sizes= s,
                      counts = cmat,
                      branches = tmat)
    
    # Stan modeling 
    stan.fit <- stan("rate-perpat-dc.stan",
                     data= data.stan, 
                     chains=1,
                    #control=list(max_treedepth=12),
                     #   adapt_delta=0.7),
                     iter=10000)
    drates<- rbind(drates, data.frame(rate=summary(stan.fit)$summary[1,6],
                                      vloop=paste0("V",as.character(loops[i])),
                                      id=type[t],
                                      lower=summary(stan.fit)$summary[1,4],
                                      upper=summary(stan.fit)$summary[1,8]
    ))
    dfits[[idx]] <- stan.fit
  }
}


# ------ Adjustment factor for indel rates ----- 
med.len <- sapply(all.lens, function(x){
  sapply(x, median)
})
med.len <- c(med.len[,1])
adj <- rep(1e3 * 365 / med.len, 2)

irates[,c(1,4,5)] <- irates[,c(1,4,5)] *  adj
drates[,c(1,4,5)] <- drates[,c(1,4,5)] * adj

irates[,c(1,4,5,6,7)] <- irates[,c(1,4,5,6,7)] *  adj
drates[,c(1,4,5,6,7)] <- drates[,c(1,4,5,6,7)] * adj

names <- c("rate", "vloop","id","lower.old","upper.old","lower","upper")
colnames(irates) <- names
colnames(drates) <- names

var <- c("V1","V2","V3","V4","V5")


# ---- plot for unfixed comparison ---- 
par(mar=c(7,6,3,2),las=1)
cl <- 1.7
ca <- 1.4
pcex <- 3.3
wid <- 3.4
plot(irates$rate, pch=1, bg="white", cex=pcex, lwd=wid,col="blue",
     cex.axis=ca, cex.lab=cl, xaxt="n",xlab="",
     ylab=expression(paste("Insertion Rate (Events/Nt/Year x",10^-3,")",sep="")),
     yaxs="i",ylim=c(0,6))
title(xlab="Group", line=5,cex.lab=cl)
axis(1, labels=rep(c("Terminal","Internal"),5), at= 1:10, line=0, cex.axis=ca-0.2)
axis(1, labels=var, at= seq(1.5,9.5, 2), line=1.5, lwd=0, cex.axis=ca)
arrows(c(1:10), irates$lower, c(1:10), irates$upper, code=3, length=0.07, angle=90,lwd=1.5)
points(x=1:10, y= iunfix$rate, bg="white", cex=pcex, lwd=wid,pch=0, col="red")
arrows(c(1:10), iunfix$lower, c(1:10), iunfix$upper, code=3, length=0.07, angle=90,lwd=1.5)
legend(x=1,y=5.7,legend=c("Fixed","Unfixed"), pch=c(1,0), pt.bg=c("white","white"),col=c("blue","red"),cex=1.5, pt.cex=pcex,pt.lwd=wid,y.intersp = 1.2)

par(mar=c(7,6,3,2),las=1)
cl <- 1.7
ca <- 1.4
pcex <- 3.3
wid <- 3.4
plot(drates$rate, pch=1, bg="white", cex=pcex, lwd=wid,col="blue",
     cex.axis=ca, cex.lab=cl, xaxt="n",xlab="",
     ylab=expression(paste("Deletion Rate (Events/Nt/Year x",10^-3,")",sep="")),
     yaxs="i",ylim=c(0,10))
title(xlab="Group", line=5,cex.lab=cl)
axis(1, labels=rep(c("Terminal","Internal"),5), at= 1:10, line=0, cex.axis=ca-0.2)
axis(1, labels=var, at= seq(1.5,9.5, 2), line=1.5, lwd=0, cex.axis=ca)
arrows(c(1:10), drates$lower, c(1:10), drates$upper, code=3, length=0.08, angle=90,lwd=1.5)
points(x=1:10, y= dunfix$rate, bg="white", cex=pcex, lwd=wid,pch=0, col="red")
arrows(c(1:10), dunfix$lower, c(1:10), dunfix$upper, code=3, length=0.08, angle=90,lwd=1.5)
legend(x=1,y=9.5,legend=c("Fixed","Unfixed"), pch=c(1,0), pt.bg=c("white","white"),col=c("blue","red"),cex=1.5, pt.cex=pcex,pt.lwd=wid,y.intersp = 1.2)





# ---- main indel rates plot ----- 

require(Rmisc)
require(ggplot2)
par(mar=c(3,3,3,2))
iplot <- ggplot() + 
  geom_bar(aes(vloop, rate, fill=id), data=irates, stat='identity', position="dodge") + 
  scale_fill_manual(values=c("red","dodgerblue"), name="Subset", labels=c("Terminal Branches","Internal Branches")) +
  coord_flip() + geom_errorbar(aes(x=irates$vloop, fill=irates$id, ymax = irates$upper, ymin = irates$lower),
                width = 0.25, size=0.8,
                position = position_dodge(0.9)) +
  scale_y_continuous(lim=c(0,10), expand=c(0,0)) + 
  scale_x_discrete(limits = rev(levels(irates$vloop))) +
  labs(#x="Variable Loop", 
       y=expression(paste("        Insertion Rate \n(Events/Nt/Year x",10^-3 ,")", sep = ""))) + 
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
        legend.position=c(0.68,0.85),
        legend.text=element_text(size=18), 
        legend.background=element_rect(colour="black"),
        legend.title=element_text(size=20),
        legend.spacing.y = unit(2, "mm")
        ) #+ geom_text(aes(y=c(1.7,1.7),x=c(3.25,2.75)),label="N/A", size=8)
#iplot
dplot <- ggplot() + 
  geom_bar(aes(x=vloop, y=rate, fill=id), data=drates, stat='identity', position="dodge") + 
  coord_flip() + scale_fill_manual(values=c("red","dodgerblue"))+
  geom_errorbar(aes(x=drates$vloop, fill=drates$id, ymax = drates$upper, ymin = drates$lower),
                width = 0.25, size=0.8,
                position = position_dodge(0.9)) +
  labs(x="Variable Loop", 
       y="Deletion Rate") +
  scale_y_reverse(lim=c(10,0), expand=c(0,0)) + 
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
        legend.position="none") #+ geom_text(aes(y=2.2,x= 3.25),label="N/A", size=8)
#dplot
multiplot(dplot,iplot, cols=2)

png(filename="~/vindels/Figures/within-host/finalized/indel-rates-v2", width=1000,height=700)
multiplot(dplot,iplot, cols=2)
dev.off()
pdf()
traceplot()
dev.off()
rstan::summary(stan.fit)

traceplot(stan.fit, alpha=0.7)
x <- extract(stan.fit)

par(mfrow=c(1,2))
plot(density(x[[1]]))
abline(v=log(t.rate), col="red")
plot(density(x[[2]]))
#abline(v=median(mat), col="red")


plot(density(w$lambda))
