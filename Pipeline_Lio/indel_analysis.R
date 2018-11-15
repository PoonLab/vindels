setwd("~/vindels/Indel_Analysis/")
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


csvfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/9_2_indels",full.names=TRUE)
max.llh <- data.frame(stringsAsFactors = F)
max.llh2 <- data.frame(subtype=character(),stringsAsFactors = FALSE)
con.int <- data.frame(stringsAsFactors = FALSE)
big.df <- data.frame(stringsAsFactors = FALSE)
rates <- c()

for (i in 1:length(csvfolder)){
  csv <- read.csv(csvfolder[i])
  
  filename <- strsplit(strsplit(csvfolder[i], "/")[[1]][7], "\\+.")[[1]][1]
  
  # if (filename == "AE"){
  #   filename <- "01_AE"
  # }else if (filename =="AG"){
  #   filename <- "02_AG"
  # }
  max.llh2 [i,1] <- filename

  print(filename)
  
  
  times.df <- c()
  outcomes.df <- data.frame()
  for (z in 1:5){
    times <- csv$total.length[which(!is.na(csv[paste0("VR",z,".indel")]))]
    times.df <- c(times.df, times)
    outcomes <- csv[which(!is.na(csv[paste0("VR",z,".indel")])),paste0("VR",z,".indel")]
    outcomes.df <- rbind(outcomes.df, data.frame(out=outcomes, vregion=z))
    
    obj.f <- function(rate) -pll(rate, outcomes, times)  # objective function
    mle.result <- bbmle::mle2(obj.f, start=list(rate=1), method = "Brent", lower=1e-12, upper = 1)
    
    max.llh <- rbind(max.llh, data.frame(subtype=filename, vloop=z, rate=coef(mle.result)[[1]])) 
    max.llh2 [i,z+1] <- coef(mle.result)[[1]]
    rates <- c(rates,coef(mle.result)[[1]]) # this is the rate
    #res$value #this is the likelihood
    if(!all(outcomes)){
      int <- confint(mle.result, level = 0.95)
      con.int <- rbind(con.int,data.frame(subtype=filename,vloop=z,lower=int[[1]],upper=int[[2]]))
      
    }else{
      con.int <- rbind(con.int,data.frame(subtype=filename,vloop=z,lower=0,upper=0))
    }
  }
  
  big.df <- rbind(big.df, data.frame(subtype=rep(filename, length(times.df)), Time=times.df, outcomes=outcomes.df$out, Vregion=as.factor(outcomes.df$vregion)))
  
}



#v.order <- c(5,1,2,3,4)
#big.df <- big.df[order(match(big.df$Vregion,v.order)),]

fit <- glm(!outcomes ~ subtype + Vregion + Time, family= "binomial", data=big.df)
fit2 <- glm(!outcomes ~ subtype * Vregion + Time, family= "binomial", data=big.df)


fit.aic <- stepAIC(fit, scope=list(upper=fit2, lower=~1), direction='both', trace=TRUE)
data.df <- data.frame(summary.glm(fit.aic)$coefficients)
data.df$Std..Error <- NULL
data.df$z.value <- NULL
data.df[which(data.df$Pr...z.. < 0.01 & data.df$Pr...z.. > 0.00143),'Signif.'] <- "**"
data.df[which(data.df$Pr...z.. < 0.00143),'Signif.'] <- "***"
data.df$Estimate <- as.numeric(as.character(format(round(data.df$Estimate, 3), nsmall = 3)))
data.df$Pr...z.. <- as.numeric(as.character(format(round(data.df$Pr...z.., 3), nsmall = 3)))



xt <- xtable(data.df, digits=c(0,3,3,0))
write.csv(max.llh, "indel_rates2.csv")
write.csv(con.int, "conf_ints2.csv")

# outcomes and times are vectors that have to come from your data frame

# + geom_rect(data=vbox, 
#             mapping=aes(xmin=vloop-1, xmax=vloop+2.25,ymin=rate, ymax=rate+0.036),
#             color="black", fill="gray88") + geom_text(data=vloops,aes(x=vloop-0.2),
#                                                       label=c("V1","V2","V3","V4","V5"),
#                                                       size=4.5) +geom_text(data=vbox, 
#                                                                            aes(x=vloop+0.3,y=rate+0.0162),
#                                                                            label="*",
#                                                                            size=10) + geom_segment(data=vbox, aes(x=vloop-0.5,xend=vloop-0.5,y=rate+0.020,yend=rate+0.016),
#                                                                                                    arrow=arrow(length=unit(2,"mm")),
#                                                                                                    size=0.8)
max.llh$subtype <- factor(max.llh$subtype, levels=c("AE", "AG", "A1", "B", "C", "D", "F1"))
max.llh <- max.llh[order(max.llh$subtype),]

max.llh2$subtype <- factor(max.llh2$subtype, levels=c("AE", "AG", "A1", "B", "C", "D", "F1"))
max.llh2 <- max.llh2[order(max.llh2$subtype),]

require(ggplot2)
#indel rate plot 
require(RColorBrewer)
colors <- brewer.pal(7, 'Set2')
colors2 <- brewer.pal(7, 'Dark2')
con.int$subtype <- factor(con.int$subtype, levels=c("AE", "AG", "A1", "B", "C", "D", "F1"))
con.int <- con.int[order(con.int$subtype),]

# MARKERS

submark <- data.frame(subtype=max.llh$subtype,vloop=max.llh$vloop,rate=rep(NaN,35))
submark[18,3] <- 0.12


rates <- split(max.llh[,2:3], max.llh[,1])
subline <- data.frame(subtype=c("B"), vloop=c(3),rate=c(0.118))
noest <- data.frame(subtype=c("F1"), vloop=c(3), rate=c(0.02))
vbox <- data.frame(subtype=c("F1"), vloop=c(3),rate=c(0.103))
vloops <- data.frame(subtype=c(rep("F1",5)), vloop=c(rep(4.5,5)),rate=c(0.100+1:5*0.007))
vline <- data.frame(subtype=c(rep(c("AE","AG","A1","B","C","D","F1"),2)), vloop=c(rep(3,7), rep(5,7)),rate=rep(0.012,7))
vline <- vline[-7,]
vline[6,3] <- 0.021
vline[2,3] <- 0.016
vline[3,3] <- 0.016



plot <- ggplot(max.llh, aes(x=vloop, 
                            y=rate, 
                            fill=subtype, 
                            width=1)) + geom_bar(colour="black",
                                                 stat="identity", 
                                                 position="dodge", 
                                                 show.legend=F) +facet_wrap(~subtype,
                                                                             ncol=7,
                                                                             nrow=1)  + geom_text(data=subline,
                                                                                                  aes(x=vloop, y=rate+0.002, label="*"),
                                                                                                  size=10) + geom_segment(data=subline, 
                                                                                                                         aes(x=vloop-2.25,xend=vloop+2.25,y=rate,yend=rate),size=0.8)  + geom_segment(data=subline, aes(x=vloop,xend=vloop,y=rate-0.003,yend=rate-0.01),
                                                                                                                                                                                                                arrow=arrow(length=unit(3,"mm")),
                                                                                                                                                                                                                size=0.8)
plot <- plot + labs(x="Variable Loop", 
            y="Indel Rate (Events/Lineage/Year)")+scale_fill_manual(values=
                                                              colors2)+scale_y_continuous(expand = c(0, 0),
                                                                                         limits = c(0, 0.155))+theme(panel.grid.major.y = element_line(color="black",size=0.3),
                                                                                                                     panel.grid.major.x = element_blank(),
                                                                                                                     panel.grid.minor.y = element_blank(),
                                                                                                                     panel.grid.minor.x = element_blank(),
                                                                                                                     panel.spacing=unit(1, "mm"),
                                                                                                                     panel.background=element_rect(fill="gray88",colour="white",size=0),
                                                                                                                     plot.margin =margin(t = 10, r = 10, b = 20, l = 8, unit = "pt"),
                                                                                                                     axis.line = element_line(colour = "black"), 
                                                                                                                     axis.title.y=element_text(size=20,margin=margin(t = 0, r = 15, b = 0, l = 0)),
                                                                                                                     axis.title.x=element_text(size=20,margin=margin(t = 15, r = 0, b = 0, l = 0)),
                                                                                                                     strip.text.x = element_text(size=16),
                                                                                                                     axis.text=element_text(size=14),
                                                                                                                     legend.position="none")+ geom_errorbar(aes(ymax = con.int$upper, ymin = con.int$lower), 
                                                                                                                                                                     width = 0.25) + geom_segment(data=vline,
                                                                                                                                                                                                  aes(x=vloop,
                                                                                                                                                                                                      y=rate-0.001,
                                                                                                                                                                                                      xend=vloop,
                                                                                                                                                                                                      yend=rate-0.005),
                                                                                                                                                                                                  arrow=arrow(length=unit(2,"mm")),
                                                                                                                                                                                                  size=0.7)+ geom_text(data=vline,
                                                                                                                                                                                                                       label="*", 
                                                                                                                                                                                                                       size=8)  + geom_text(aes(y=rate-0.005),
                                                                                                                                                                                                                                            data=noest, 
                                                                                                                                                                                                                                            label="no estimate", 
                                                                                                                                                                                                                                            size=6, 
                                                                                                                                                                                                                                            angle=90)

plot

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
