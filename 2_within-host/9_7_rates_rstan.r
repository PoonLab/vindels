# stan model for indel rates 
library(rstan)

num.pat <- 20
num.trees <- 30

# ---- Simulated Data ----
tru.rate <- 0.3
tru.sd <- 0.02
pat.rate <- rnorm(num.pat, mean=tru.rate, sd=tru.sd)
pat.sd <- 0.05



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
  toCheck[idx] <- tre.rate
  for (j in 1:num.trees){
    len.temp <- rexp(diff(range(idx))+1, 0.1)
    lens[idx, j] <- len.temp
    c <- rpois(diff(range(idx))+1,sizes=85, mu=tre.rate[j]*len.temp)
    if(any(is.na(c))){
      print(tre.rate[j])
      print(len.temp)
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
stan.fit <- stan("~/Desktop/rate-perpat-dc.stan",
                 data= data.stan, 
                 chains=1,
                 iter=10000)


# TESTING
nrows <- 100
tre.rate <- rnorm(num.trees, mean=0.02, sd=0.005)
counts <- matrix(nrow=nrows, ncol=num.trees)
lens <- matrix(nrow=nrows, ncol=num.trees)
for (j in 1:num.trees){
  len.temp <- rexp(nrows, 1.4)
  lens[, j] <- len.temp
  c <- rpois(nrows, lambda=tre.rate[j]*len.temp)
  counts[,j] <- c
}

data.stan <- list(ntree = num.trees,
                  nrow=nrows,
                  counts = counts,
                  lengths = lens)

 
stan.fit <- stan("~/vindels/2_within-host/rate-test.stan",
                 data= data.stan, 
                 chains=1,
                 iter=100000)    



# ----- Real Data ----
type <- c("tip","node")
rstan_options(auto_write=T)
options(mc.cores = parallel::detectCores()-2)

for (t in 1:2){ 
  for (i in 1:5){
    cmat <- all.counts[[t]][[i]]
    tmat <- all.times[[t]][[i]]
    
    s <- sizes[[t]]
    
    num.pat <- length(s)
    num.tree <- ncol(cmat)

      
    # Data import
    data.stan <- list(npat=num.pat,
                      ntree = num.tree,
                      sizes= s,
                      counts = cmat,
                      lengths = tmat)
      
    # Stan modeling 
    stan.fit <- stan("~/Desktop/rate-perpat-dc.stan",
                     data= data.stan, 
                     chains=1,
                     iter=10000,
                     )
    #control = list(adapt_delta = 0.99))
    
    # export results to a data frame 
    irates <- rbind(irates, data.frame(rate=summary(stan.fit)$summary[1,6],
                                       vloop=paste0("V",as.character(i)),
                                       id=type[t],
                                       lower=summary(stan.fit)$summary[1,4],
                                       upper=summary(stan.fit)$summary[1,8]
    ))      # ----- Stan Model ----
    # }else{
    #   irates <- rbind(irates, data.frame(rate=NA,
    #                                      vloop=paste0("V",as.character(i)),
    #                                      id=type[t],
    #                                      lower=NA,
    #                                      upper=NA))
    # }

    
    data <- vloop.mat[[i]][[t+2]]
    num.pat <- ncol(data)
    num.data <- nrow(data)
    if (!is.na(data)){
      # ----- Stan Model ----
      data.stan <- list(npat = num.pat,
                        ndata = num.data,
                        mat = data)
      
      stan.fit <- stan("~/vindels/2_within-host/rates.stan",
                       data= data.stan, 
                       chains=1,
                       iter=1000000)
                       #control = list(adapt_delta = 0.99))
      
      drates<- rbind(drates, data.frame(rate=summary(stan.fit)$summary[1,6],
                                        vloop=paste0("V",as.character(i)),
                                        id=type[t],
                                        lower=summary(stan.fit)$summary[1,4],
                                        upper=summary(stan.fit)$summary[1,8]
      ))
    }else{
      drates <- rbind(drates, data.frame(rate=NA,
                                         vloop=paste0("V",as.character(i)),
                                         id=type[t],
                                         lower=NA,
                                         upper=NA))
    }
  }
}
# # --- bootstraps --- 
# not useful; give CIs that are far too narrow
# bs <- c()
# for (i in 1:1000){
#   cols <- ncol(V1.mat[[1]])
#   bs[i] <- mean(V1.mat[[1]][sample(1:200, cols),1:cols])
# }


require(Rmisc)
require(ggplot2)
iplot <- ggplot() + 
  geom_bar(aes(vloop, rate, fill=id), data=irates, stat='identity', position="dodge") + 
  scale_fill_manual(values=c("red","dodgerblue"), name="Subset", labels=c("Tip Sequences","Internal Nodes")) +
  coord_flip() + geom_errorbar(aes(x=irates$vloop, fill=irates$id, ymax = irates$upper, ymin = irates$lower),
                width = 0.25, size=0.8,
                position = position_dodge(0.9)) +
  scale_y_continuous(lim=c(0,35), expand=c(0,0)) + 
  scale_x_discrete(limits = rev(levels(irates$vloop))) +
  labs(#x="Variable Loop", 
       y=expression(paste("      Insertion Rate \n (Events/Nt/Year x",10^-3 ,")", sep = ""))) + 
       #title="Indel Rates",
       #color="Subset") +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),#element_line(color="black",size=0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 4, r = 4, b = 5, l = 0, unit = "mm"),
        panel.border=element_rect(fill=NA, size=1),
        #axis.line.x = element_line(colour = "black"), 
        #axis.line.y.left=element_line(colour="black"),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=24,margin=margin(t = 22, r = 3, b = 0, l = 22)),
        axis.text.y = element_text(size=22, colour="black", margin=margin(t = 0, r = 6, b = 2, l = 2)),
        axis.text.x=element_text(size=22, colour="black"),
        #plot.title = element_text(size=22, hjust = 0.5),
        legend.position=c(0.73,0.85),
        legend.text=element_text(size=22), 
        legend.background=element_rect(colour="black"),
        legend.title=element_text(size=24),
        legend.spacing.y = unit(2, "mm")
        ) + geom_text(aes(y=c(1.7,1.7),x=c(3.25,2.75)),label="N/A", size=8)
iplot
dplot <- ggplot() + 
  geom_bar(aes(x=vloop, y=rate, fill=id), data=drates, stat='identity', position="dodge") + 
  coord_flip() + scale_fill_manual(values=c("red","dodgerblue"))+
  geom_errorbar(aes(x=drates$vloop, fill=drates$id, ymax = drates$upper, ymin = drates$lower),
                width = 0.25, size=0.8,
                position = position_dodge(0.9)) +
  labs(x="Variable Loop", 
       y="Deletion Rate") +
  scale_y_reverse(lim=c(35,0), expand=c(0,0)) + 
  scale_x_discrete(limits = rev(levels(drates$vloop))) +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),#element_line(color="black",size=0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        panel.border=element_rect(fill=NA, size=1),
        #axis.line.x=element_line(colour = "black"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 4, r = 2, b = 13.5, l = 12, unit = "mm"),
        #axis.line.y = element_line(colour = "black"), 
        axis.text.y = element_blank(),
        axis.title.y=element_text(size=22),
        axis.title.x=element_text(size=24,margin=margin(t = 7, r = 3, b = 0, l = 12)),
        axis.text = element_text(size=22, colour="black"),
        plot.title = element_text(size=28, hjust = 0.5),
        legend.position="none") + geom_text(aes(y=2.2,x= 3.25),label="N/A", size=8)
dplot
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
