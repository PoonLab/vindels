# stan model for indel rates 
library(rstan)


# will receive data in from 9_5 indel rates 



# ---- Simulated Data ----
t.rate <- 2.5
sub.rate <- rlnorm(num.pat, meanlog=log(t.rate), sdlog=0.3)

hist(sub.rate, breaks=5, col="grey")

mat <- matrix(nrow=num.data, ncol=num.pat)

for (i in 1:num.pat){
  mat[,i] <- rnorm(num.data, mean=sub.rate[i], sd=0.5)
}

# ----- Real Data ----
for (i in 1:5){
  data <- vloop.mat[[i]][[4]]
  num.pat <- ncol(data)
  num.data <- 200
  rstan_options(auto_write=T)
  options(mc.cores = parallel::detectCores()-2)
  
  
  
  # ----- Stan Model ----
  data.stan <- list(npat = num.pat,
                    ndata = num.data,
                    mat = data)
  
  stan.fit <- stan("~/vindels/2_within-host/rates.stan",
                   data= data.stan, 
                   chains=1,
                   iter=1000000)
  
  #stan_dens(stan.fit)
  summary(stan.fit)$summary
  #irates <- data.frame()
  # irates <- rbind(irates, data.frame(rate=summary(stan.fit)$summary[1,6],
  #                                    vloop=i,id="node",
  #                                    lower=summary(stan.fit)$summary[1,4],
  #                                    upper=summary(stan.fit)$summary[1,8]
  #                                    ))
  drates<- rbind(drates, data.frame(rate=summary(stan.fit)$summary[1,6],
                                    vloop=i,id="node",
                                    lower=summary(stan.fit)$summary[1,4],
                                    upper=summary(stan.fit)$summary[1,8]
  ))
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
  scale_fill_manual(values=c("red","blue"), name="Subset", labels=c("Tip Sequences","Internal Nodes")) +
  geom_errorbar(aes(x=irates$vloop, fill=irates$id, ymax = irates$upper, ymin = irates$lower),
                width = 0.25, size=0.8,
                position = position_dodge(0.9)) +
  labs(x="Variable Loop", 
       y=expression(paste("      Insertion Rate \n (Events/Nt/Year x  ",10^-3 ,")", sep = "")), 
       title="Indel Rates",
       color="Subset") +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 4, r = 1, b = 0, l = 2.5, unit = "mm"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=18,margin=margin(t = 9, r = 3, b = 0, l = 22)),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=16, colour="black"),
        axis.text.x=element_blank(),
        plot.title = element_text(size=22, hjust = 0.5),
        legend.text=element_text(size=16), 
        legend.background=element_rect(colour="black"),
        legend.title=element_text(size=18),
        legend.spacing.y = unit(2, "mm")
        ) + geom_text(aes(y=0.5,x= 3.25),label="N/A", size=6)
iplot
dplot <- ggplot() + 
  geom_bar(aes(x=vloop, y=rate, fill=id), data=drates, stat='identity', position="dodge") + 
  scale_fill_manual(values=c("red","blue"))+
  geom_errorbar(aes(x=drates$vloop, fill=drates$id, ymax = drates$upper, ymin = drates$lower),
                width = 0.25, size=0.8,
                position = position_dodge(0.9)) +
  labs(x="Variable Loop", 
       y="Deletion Rate") +
  scale_y_reverse(lim=c(10,0)) + 
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = -2, r = 55.5, b = 4.5, l = 13.3, unit = "mm"),
        axis.line = element_line(colour = "black"), 
        axis.title=element_text(size=18,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        axis.text = element_text(size=16, colour="black"),
        plot.title = element_text(size=22, hjust = 0.5),
        legend.position="none") + geom_text(aes(y=0.5,x= 3.25),label="N/A", size=6)
dplot

multiplot(iplot,dplot)

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
