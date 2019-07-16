


acute <- read.csv("~/PycharmProjects/glyc-analysis/9_glycs/acute.csv")
chronic <- read.csv("~/PycharmProjects/glyc-analysis/9_glycs/chronic.csv")

rownames(acute) <- acute$position
rownames(chronic) <- chronic$position

acute$prop <- acute$count / sum(acute$count)
chronic$prop <- chronic$count / sum(chronic$count)

combined <- c(acute$position, chronic$position)
dupl <- combined[duplicated(combined)]

toPlot <- data.frame(Acute=acute[acute$position %in% dupl, "prop"], Chronic=chronic[chronic$position %in% dupl, "prop"])
par(pty="s", mar=c(5,8,4,1))
lim <- c(0,0.04)
plot(toPlot, xlim=lim, ylim=lim, cex.lab=1.4, cex.main=1.6, cex.axis=1.1, main="N-Glyc site prevalence")
abline(0,1)

acute$unique <- unname(sapply(acute$position, function(x){x %in% dupl}))



# INDEL ABOVE/BELOW MULTIPLOTS 
# -------------------------------------

require(Rmisc)
require(ggplot2)

g1 <- ggplot(toPlot2, aes(x=Position, y=Acute)) + 
  geom_bar(colour="black", stat="identity",fill="dodgerblue",position="dodge",show.legend=F)+
  labs(x="Position", 
       y="")+
  ylim(0,250)+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 30, r = 10, b = 0, l = 18, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=20,margin=margin(t = 0, r = 3, b = 5, l = 12)),
        axis.title.x=element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        legend.position="none") + geom_text(aes(y=225,x=30 ),
                                          label="Acute", 
                                          size=8)
#g1



g2 <- ggplot(toPlot2, aes(x=Position, y=Chronic)) + 
  geom_bar(colour="black", stat="identity",fill="dodgerblue",position="dodge",show.legend=F)+
  labs(x="Position", 
       y="                                         Count")+
  scale_y_reverse(lim=c(250,0))+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = -10, r = 10, b = 15, l = 18, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=20,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        axis.title.x=element_text(size=20,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        strip.text.x = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.position="none")+ geom_text(aes(y=175,x=35),
                                           label="Chronic", 
                                           size=8)
#g2

multiplot(g1,g2)

g2 <- ggplot(delrates, aes(x=vloop, y=rate,width=0.8)) + 
  geom_bar(colour="black", stat="identity",fill="firebrick1",position="dodge",show.legend=F) + 
  geom_errorbar(aes(ymax = delrates$upper, ymin = delrates$lower), 
                width = 0.25, size=1.1) +
  labs(x="Variable Loop", 
       y="Deletion Rate")+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_y_reverse(lim=c(5,0))+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = -10, r = 10, b = 8, l = 24, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=18,margin=margin(t = 0, r = 11, b = 0, l = 6)),
        axis.title.x=element_text(size=18,margin=margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size=16),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=14),
        legend.position="none") + geom_text(aes(y=0.3,x=3 ),
                                            label="N/A", 
                                            size=6)


