# stan model for indel rates 
library(rstan)
setwd("~/vindels/2_within-host/")
num.pat <- 20
num.trees <- 30

# ---- Simulated Data ----
tru.rate <- 0.03
tru.sd <- 0.002
pat.rate <- rnorm(num.pat, mean=tru.rate, sd=tru.sd)
pat.sd <- 0.008



hist(pat.rate, breaks=5, col="grey")
vlens <- c(75, 126, 105,  90 , 33)
lns <- sample(70:200, num.pat, replace=T)
counts <- matrix(nrow=sum(lns), ncol=num.trees)
lens <- matrix(nrow=sum(lns), ncol=num.trees)
pos <- 1
toCheck <- c()
for (i in 1:num.pat){
  idx <- pos:(pos+lns[i]-1)
  print(range(idx))
  tre.rate <- rnorm(num.trees, mean=pat.rate[i], sd=pat.sd)
  toCheck[((i-1)*30+1):(i*30)] <- tre.rate
  for (j in 1:num.trees){
    len.temp <- rexp(diff(range(idx))+1, 0.1)
    lens[idx, j] <- len.temp
    c <- rpois(diff(range(idx))+1, lambda=tre.rate[j]*len.temp)
    if(any(is.na(c))){
      print(tre.rate[j])
      #print(len.temp)
    }
    counts[idx,j] <- c
  }
  pos <- pos + lns[i]
}

rstan_options(auto_write=T)
options(mc.cores = parallel::detectCores()-2)

# Data import
data.stan <- list(npat= num.pat,
                  ntree = num.trees,
                  sizes = lns,
                  counts = counts,
                  lengths = lens)

# Stan modeling 
stan.fit <- stan("rate-perpat-dc.stan",
                 data= data.stan, 
                 chains=1,
                 iter=100000)
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

rm(all.data)
rm(mat)
rm(tre)
rm(ifolder)
rm(dfolder)
rm(lens)
rm(maxes)
rm(patnames)
rm(tips)
rm(selectRows)
rm(pat.idx)


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
                     iter=40000)
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
# 
# 
# # ---- just tips ----
# mean / median
# irate <- c(2.382861e-05, 2.338753e-05, 4.196071e-05, 4.180047e-05, 4.431492e-05, 4.498296e-05, 8.006266e-05, 8.021932e-05)
# isd <- c(1.483429e-06 ,1.241316e-06 ,1.496752e-06 ,1.273513e-06 ,1.525796e-06 ,1.264572e-06 ,1.522098e-06 ,1.261372e-06)
# 
# drate <- c(7.238429e-05, 7.232098e-05, 5.606446e-05, 5.563663e-05, 6.847116e-05, 6.826069e-05, 1.051018e-04, 1.061269e-04)
# dsd<-  c(1.481917e-06, 1.232391e-06, 1.492279e-06, 1.256041e-06, 1.507554e-06, 1.268980e-06, 1.535336e-06, 1.278238e-06)
# 


isd.mean <- c()
isd.med <- c()
dsd.mean <- c()
dsd.med <- c()

for (j in 1:10){
  if (j == 5){
    isd.mean[j] <- 1.522098e-06 
    isd.med[j] <- 1.261372e-06
    dsd.mean[j] <- 1.535336e-06 
    dsd.med[j] <- 1.278238e-06
  }else{
    if(j < 5){
      idx <- j
    }else{
      idx <- j - 1
    }
    isd.mean[j] <- rstan::summary(ifits[[idx]])$summary[2,1]
    isd.med[j] <- rstan::summary(ifits[[idx]])$summary[2,6]
    dsd.mean[j] <- rstan::summary(dfits[[idx]])$summary[2,1]
    dsd.med[j] <- rstan::summary(dfits[[idx]])$summary[2,6]
  }
  
}


# --- factor of standard deviation estimated from subset --- 
isd <- c(1.428940, 1.497948, 1.468787, 1.476039, 1.457945, 1.554781 ,1.319389 ,1.448355, 1.450687, 1.476578)
dsd <- c(1.472111, 1.502671, 1.484049, 1.469913, 1.477119, 1.383553, 1.499856, 1.506177, 1.488022, 1.487717)

# use the factor to get the SD of the true data set 
isd <- isd * 1e-5
dsd <- dsd * 1e-5

irates$lower.sd <- irates$rate - isd.mean*1.96
irates$upper.sd <- irates$rate + isd.mean*1.96
drates2$lower.sd <- drates2$rate - dsd.med*1.96
drates2$upper.sd <- drates2$rate + dsd.med*1.96

irates[c(3,8),c(6,7)] <- irates[c(3,8),c(4,5)] 
drates[c(3,8),c(6,7)] <- drates[c(3,8),c(4,5)] 

# --- save unfixed data to combine with subset ----

rate <- c(0.120081126915412,0.122842191549769,0.0650445666093929,0.189200446798185,1.03103002736818,0.173485720631059,0.204177862078437,0.120095020342278,0.150454526809474,1.64025170330846) 
lower <- c(0.0818062745261698,0.104004964011112,0.040600743707683,0.156652536273797,0.853595423785998,0.125695208689291,0.16103064244673,0.083219085725915,0.112024276849568,1.37427701063287)
upper <- c(0.160401151109888,0.150362872181409,0.0879200762269469,0.246031483063917,1.19453080501549,0.214411233850892,0.244126807621243,0.165792223855134,0.214776685102334,1.93844426646181)

rate <- c(0.336917532403483,0.172791614066862,0.0519541242881352,0.264269772591424,1.23725613363763,0.0912730560723649,0.172873897935259,0.155197076546366,0.213240930648582,4.34175824513442)
lower <- c(0.274793180553731,0.147491316890742,0.0352698170527366,0.223330008781467,1.05809490155911,0.0476396918003235,0.139853798267994,0.130917696268314,0.17445588434802,3.87255641362116)
upper <- c(0.37667685335889,0.205680539066375,0.0682047654275175,0.320520787004128,1.54391037458851,0.149976880144637,0.212014163604103,0.188816629201296,0.269128696342952,4.63668210730703)





isd <- isd/1e-6
dsd <- dsd/1e-6

med.len <- c(78, 123, 105,  93,  30)


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
