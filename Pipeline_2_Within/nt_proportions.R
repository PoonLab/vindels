# INSERTIONS nucleotide proportions

require(stringr)
ifolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/10_nucleotides/ins/*.csv")
dfolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/10_nucleotides/del/*.csv")

charCount <- function(str){
  if (str == ""){
    print(str)
  }
  else{
    print(str)
  }
}

removeNA <- function(input){
  if (is.na(input)){
    input <- ""
  }
  input
}

ntcount <- c()
all.ins <- data.frame(stringsAsFactors = F)
all.del <- data.frame(stringsAsFactors = F)
for (file in 1:length(ifolder)){
  iCSV <- read.csv(ifolder[file], na.strings=c("NA"),row.names = 1, stringsAsFactors = F)
  dCSV <- read.csv(dfolder[file], na.strings=c("NA"),row.names = 1, stringsAsFactors = F)
  
  if (all(!is.na(iCSV$Seq))){
    iCSV$Seq <- sapply(iCSV$Seq, removeNA)
    dCSV$Seq <- sapply(dCSV$Seq, removeNA)
    
    colnames(iCSV) <- c("Accno", "Vloop", "Seq", "VSeq")
    colnames(dCSV) <- c("Accno", "Vloop", "Seq", "VSeq")
    
    all.ins <- rbind(all.ins, iCSV)
    all.del <- rbind(all.del, dCSV)
  }
}

all.ins <- all.ins[all.ins$Seq!="",]
all.del <- all.del[all.del$Seq!="",]




ins <- list()
del <- list()
nucleotides <- c("A","C","G","T")

ins <- data.frame(nucl=nucleotides)
del <- data.frame(nucl=nucleotides)

# a vector of two totals
# iTotals[1] = total number of nucleotides in insertion sequences
# iTotals[2] = total nuimebr of nucleotides in the vloops associated with insertion sequences 
iTotals <- c(sum(unname(sapply(all.ins$Seq, nchar))), sum(unname(sapply(all.ins$VSeq, nchar))))
dTotals <- c(sum(unname(sapply(all.del$Seq, nchar))), sum(unname(sapply(all.del$VSeq, nchar))))

iProps <- c()
dProps <- c()

iVProps <- c()
dVProps <- c()
for (nuc in nucleotides){
  iProps <- c(iProps, sum(str_count(all.ins$Seq, nuc)) / iTotals[1])
  dProps <- c(dProps, sum(str_count(all.del$Seq, nuc)) / dTotals[1])
  
  iVProps <- c(iVProps, sum(str_count(all.ins$VSeq, nuc)) / iTotals[2])
  dVProps <- c(dVProps, sum(str_count(all.del$VSeq, nuc)) / dTotals[2])
}



ins.props <- data.frame(nt=nucleotides, iprops=iProps, vloop=iVProps)
del.props <- data.frame(nt=nucleotides, dprops=dProps, vloop=dVProps)



require(RColorBrewer)
colors <- brewer.pal(4, "Set1")

#A
cex=2
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0.10,0.45)
plot(ins.props[,c(3,2)], pch=21, bg=colors,xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Insertions")
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Inside Insertions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.38,0.22,legend=nucleotides, pch=21,cex=1.9, pt.bg=colors,x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=F)
abline(0,1)






cex=2
par(pty="s", xpd=NA, mar=c(6,8,4,1),las=0)

lim = c(0.10,0.45)
plot(del.props[,c(3,2)], pch=21, bg=colors,xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2, ylab='', xlab='',cex=3.5, main="Deletions")
#text(0.187,0.475,labels="a)", cex=1.5)
#text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Inside Deletions", line=3.5,cex.lab=1.75)
title(xlab="Proportion in Variable Loops", line=3.5,cex.lab=1.75)
legend(0.38,0.22,legend=nucleotides, pch=21,cex=1.9, pt.bg=colors,x.intersp = 1.0,y.intersp=1.0, pt.cex=3)
par(xpd=F)
abline(0,1)

require(ggplot2)

plot <- ggplot(ntcount, aes(x=nt, 
                               y=props,
                               width=1)) + geom_bar(colour="black",
                                                    stat="identity", 
                                                    fill="dodgerblue",
                                                    position="dodge", 
                                                    show.legend=F) 

plot <- plot + labs(x="Nucleotide", 
                    y="Proportion in Insertions")+scale_y_continuous(expand = c(0, 0),
                                                                                                                                                     limits = c(0, 1))+theme(panel.grid.major.y = element_line(color="black",size=0.3),
                                                                                                                                                                               panel.grid.major.x = element_blank(),
                                                                                                                                                                               panel.grid.minor.y = element_blank(),
                                                                                                                                                                               panel.grid.minor.x = element_blank(),
                                                                                                                                                                               panel.spacing=unit(1, "mm"),
                                                                                                                                                                               #panel.background=element_rect(fill="gray88",colour="white",size=0),
                                                                                                                                                                               plot.margin =margin(t = 20, r = 20, b = 20, l = 8, unit = "pt"),
                                                                                                                                                                               axis.line = element_line(colour = "black"), 
                                                                                                                                                                               axis.title.y=element_text(size=20,margin=margin(t = 0, r = 15, b = 0, l = 0)),
                                                                                                                                                                               axis.title.x=element_text(size=20,margin=margin(t = 15, r = 0, b = 0, l = 0)),
                                                                                                                                                                               strip.text.x = element_text(size=16),
                                                                                                                                                                               axis.text=element_text(size=14),
                                                                                                                                                                               legend.position="none")



#A
cex=2
par(pty="s", mfrow=c(2,2), xpd=NA, mar=c(3,8,4,1),las=0)

lim = c(0.24,0.46)
plot(list.df[[1]][,1:2], cex=sizes.v, pch=(list.df[[1]]$vregion+20), bg=rep(colors, 7),xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=0.9,cex.main=1.5, ylab='', xlab='')
text(0.187,0.475,labels="a)", cex=1.5)
text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Within Indels", line=2.5,cex.lab=1.15)
title(xlab="Proportion Outside Indels", line=2.1,cex.lab=1.15)
par(xpd=F)
abline(0,1)
