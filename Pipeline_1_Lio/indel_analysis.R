require(bbmle)
require(MASS)
require(xtable)
obj.f <- function(rate) -pll(rate, outcomes, times)  # objective function

pll <- function(rate, outcomes, times) {
  # @param rate:  instantaneous rate of indels
  # @param outcomes:  vector of boolean values, TRUE if *no* indel on cherry
  # @param times:  vector of branch lengths associated with cherries
  # @return: log-likelihood of model
  pr <- exp(-rate * times)  # probabilities of no indels for each cherry
  sum( outcomes * log(pr), (1-outcomes) * log(1-pr) )
}


#csvfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/9_2_indels",full.names=TRUE)
csvfolder <- Sys.glob("~/PycharmProjects/hiv-evolution-master/9_2_indels/*.csv")
max.llh <- data.frame(stringsAsFactors = F)
max.llh2 <- data.frame(subtype=character(),stringsAsFactors = FALSE)
con.int <- data.frame(stringsAsFactors = FALSE)
big.df <- data.frame(stringsAsFactors = FALSE)
vlen.tot <- data.frame()
ml.output <- list()
clades <- c()

for (i in 1:length(csvfolder)){
  csv <- read.csv(csvfolder[i])
  
  filename <- strsplit(strsplit(csvfolder[i], "/")[[1]][7], "\\+.")[[1]][1]
  clades[i] <- filename
  # if (filename == "AE"){
  #   filename <- "01_AE"
  # }else if (filename =="AG"){
  #   filename <- "02_AG"
  # }
  max.llh2 [i,1] <- filename

  print(filename)
  
  data.df <- data.frame()
  for (z in 1:5){
    times <- csv[which(!is.na(csv[paste0("VR",z,".indel")])),6]  #retrieves total.length
    outcomes <- csv[which(!is.na(csv[paste0("VR",z,".indel")])),paste0("VR",z,".indel")]
    vlen <- csv[which(!is.na(csv[paste0("VR",z,".indel")])), paste0("VR",z,".len")]
    
    data.df <- rbind(data.df, data.frame(out=outcomes, times=times, vregion=z, vlen=vlen))
    
    vlen.tot <- rbind(vlen.tot, data.frame(clade=filename, vloop=z, vlength=mean(vlen)))
    
    obj.f <- function(rate) -pll(rate, outcomes, times)  # objective function
    mle.result <- bbmle::mle2(obj.f, start=list(rate=1), method = "Brent", lower=1e-12, upper = 1)
    
    rate <- coef(mle.result)[[1]]
    
    max.llh <- rbind(max.llh, data.frame(subtype=filename, vloop=z, rate=rate, rate2=rate/mean(vlen), adj.rate=1000*(rate/mean(vlen)))) 
    max.llh2[i,z+1] <- rate/mean(vlen)
    
    if(!all(outcomes) & !all(!outcomes)){
      int <- confint(mle.result, level = 0.95)
      con.int <- rbind(con.int,data.frame(subtype=filename,vloop=z,lower=1000*(int[[1]]/mean(vlen)),upper=1000*(int[[2]]/mean(vlen))))
      
    }else{
      con.int <- rbind(con.int,data.frame(subtype=filename,vloop=z,lower=0,upper=0))
    }
    
  }
  
  big.df <- rbind(big.df, data.frame(subtype=rep(filename, nrow(data.df)), Time=data.df$times, outcomes=data.df$out, Vregion=as.factor(data.df$vregion), Vlength=data.df$vlen))
  
}


max.llh$subtype <- factor(max.llh$subtype, levels=c("AE", "AG", "A1", "B", "C", "D", "F1"))
max.llh <- max.llh[order(max.llh$subtype),]



max.llh2$subtype <- factor(max.llh2$subtype, levels=c("AE", "AG", "A1", "B", "C", "D", "F1"))
max.llh2 <- max.llh2[order(max.llh2$subtype),]

big.df$subtype <- factor(big.df$subtype, levels=c("AE", "AG", "A1", "B", "C", "D", "F1"))
big.df <- big.df[order(big.df$subtype),]

avg <- c()
for (y in 1:5){
  avg[y] <- median(max.llh[which(max.llh$vloop==y),5])
}


big.df$time.len <- big.df$Time * big.df$Vlength

fit <- glm(!outcomes ~ subtype + Vregion + time.len, family= "binomial", data=big.df)
fit2 <- glm(!outcomes ~ subtype * Vregion + time.len, family= "binomial", data=big.df)
fit.aic <- stepAIC(fit, scope=list(upper=fit2, lower=~1), direction='both', trace=TRUE)


data.df <- data.frame(summary.glm(fit.aic)$coefficients)
data.df$Std..Error <- NULL
data.df$z.value <- NULL
data.df[which(data.df$Pr...z.. < 0.01 & data.df$Pr...z.. > 0.00143),'Signif.'] <- "*"
data.df[which(data.df$Pr...z.. < 0.00143),'Signif.'] <- "***"
data.df$Estimate <- as.numeric(as.character(format(round(data.df$Estimate, 3), nsmall = 3)))
data.df$Pr...z.. <- as.numeric(as.character(format(round(data.df$Pr...z.., 3), nsmall = 3)))

xt <- xtable(data.df, digits=c(0,3,3,0))
write.csv(max.llh, "indel_rates2.csv")
write.csv(con.int, "conf_ints2.csv")




require(ggplot2)
#indel rate plot 
require(RColorBrewer)
colors <- brewer.pal(7, 'Set2')
colors2 <- brewer.pal(7, 'Dark2')
con.int$subtype <- factor(con.int$subtype, levels=c("AE", "AG", "A1", "B", "C", "D", "F1"))
con.int <- con.int[order(con.int$subtype),]

max.llh$subtype <- factor(max.llh$subtype, levels= c("AE", "AG", "A1", "B", "C", "D", "F1"))

# MARKERS

subline <- data.frame(subtype=c("B"), vloop=c(3),adj.rate=c(1.3))
noest <- data.frame(subtype=c("F1"), vloop=c(3), adj.rate=c(0.2))
vline <- data.frame(subtype=c("AE","AG","A1","B","C","D","F1"), vloop=rep(3,7),adj.rate=rep(0.4,7))
vline <- rbind(vline, data.frame(subtype=c("B", "C"), vloop=c(2,2), adj.rate=c(0.4,0.4)))

# vline[3,3] <- 0.016



plot <- ggplot(max.llh, aes(x=vloop,    #*************************************************************************************************************************
                            y=adj.rate, 
                            fill=subtype, 
                            width=1)) + geom_bar(colour="black",
                                                 stat="identity", 
                                                 position="dodge", 
                                                 show.legend=F) +facet_wrap(~subtype,
                                                                             ncol=7,
                                                                             nrow=1)  + geom_text(data=subline,
                                                                                                  aes(x=vloop, y=adj.rate+0.03, label="*"),
                                                                                                  size=9) + geom_segment(data=subline, 
                                                                                                                         aes(x=vloop-2.25,xend=vloop+2.25,y=adj.rate,yend=adj.rate),size=0.8)  + geom_segment(data=subline, aes(x=vloop,xend=vloop,y=adj.rate-0.03,yend=adj.rate-0.11),
                                                                                                                                                                                                                arrow=arrow(length=unit(4,"mm")),
                                                                                                                                                                                                                size=0.8)
plot <- plot + labs(x="Variable Region", 
            y=expression(paste("Indel Rate (Events/Nt/Year x ", 10^-3, ")", sep = "")))+scale_fill_manual(values=colors2)+scale_y_continuous(expand = c(0, 0),
                                                                                                       limits = c(0, 2))+theme(panel.grid.major.y = element_line(color="black",size=0.3),
                                                                                                                     panel.grid.major.x = element_blank(),
                                                                                                                     panel.grid.minor.y = element_blank(),
                                                                                                                     panel.grid.minor.x = element_blank(),
                                                                                                                     panel.spacing=unit(5, "mm"),
                                                                                                                     panel.background=element_rect(fill="gray88",colour="white",size=0),
                                                                                                                     plot.margin=margin(t = 10, r = 18, b = 18, l = 20, unit = "pt"),
                                                                                                                     axis.line = element_line(colour = "black"), 
                                                                                                                     axis.title.y=element_text(size=28,margin=margin(t = 0, r = 22, b = 0, l = 0), colour="black"),
                                                                                                                     axis.title.x=element_text(size=30,margin=margin(t = 22, r = 0, b = 0, l = 0), colour="black"),
                                                                                                                     strip.text.x = element_text(size=26),
                                                                                                                     axis.text=element_text(size=24,colour="black"),
                                                                                                                     legend.position="none")+ geom_errorbar(aes(ymax = con.int$upper, ymin = con.int$lower), 
                                                                                                                                                                     width = 0.25) + geom_segment(data=vline,
                                                                                                                                                                                                  aes(x=vloop,
                                                                                                                                                                                                      y=adj.rate-0.01,
                                                                                                                                                                                                      xend=vloop,
                                                                                                                                                                                                      yend=adj.rate-0.06),
                                                                                                                                                                                                  arrow=arrow(length=unit(3,"mm")),
                                                                                                                                                                                                  size=0.9)+ geom_text(data=vline,
                                                                                                                                                                                                                       aes(x=vloop, y=adj.rate+0.02, label=c("†","†","†","†","‡","†","‡","†","†")), 
                                                                                                                                                                                                                       size=10)  + geom_text(aes(y=adj.rate-0.02),
                                                                                                                                                                                                                                            data=noest, 
                                                                                                                                                                                                                                            label="N/A", 
                                                                                                                                                                                                                                            size=8, 
                                                                                                                                                                                                                                            angle=90)
plot  #**********************************************************************************************************************************************************************************************************************

figure <- ggplot_build(plot)
figure$layout$clip[figure$layout$name=="panel"] <- "off"
figure2 <- figure 
segments((n*7)+0.5,8.5,(n*7)+0.5,-1.4,lwd=3)
#plot +  geom_errorbar(limits, position = position_dodge(0.9),width = 0.25)

pos = c(1:42)
  
axis(1,at=c(3.5,9.5,15.5,21.5,27.5,33.5,39.5), lab=c("AE","AG","A1","B","C","D","F1"),line=1,tick=F,cex.axis=1.3)

barplot(max.llh$rate, col=rep(colors,7), space=rep(c(1,0,0,0,0),7), ylim=c(0,0.08),xaxt='n')
rect(7,0.00,8,0.02,col=rgb(211/255,211/255,211/255))
text(7.5,0.002,adj=0,labels="NO ESTIMATE GIVEN",cex=1,srt=90)
rect(39,0.00,40,0.02,col=rgb(211/255,211/255,211/255))
text(39.5,0.002,adj=0,labels="NO ESTIMATE GIVEN",cex=1,srt=90)


axis(1,at=c(3.5,9.5,15.5,21.5,27.5,33.5,39.5),lab=c("AE","AG","A1","B","C","D","F1"), line=2, tick=F)
axis(1,at=pos2+0.5, lab=rep(c("1","2","3","4","5"),7),line=-1,tick=F,cex.axis=1)



# handle.warning <- function(csv, x) {
#   fit <- tryCatch(
#     {fit <- glm(!csv[,as.integer(2*x+4)] ~ csv$length, family= 'binomial')},
#     warning = function(c){
#       message("A warning was thrown, still running")
#       return (NULL)
#     })
#   return(fit)
# }
