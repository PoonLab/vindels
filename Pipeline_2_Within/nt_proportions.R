# INSERTIONS nucleotide proportions

require(stringr)
csvfolder <- Sys.glob("~/Lio/*.csv")

charCount <- function(str){
  if (str == ""){
    print(str)
  }
  else{
    print(str)
  }
}

ntcount <- c()
total <- data.frame()
for (file in 1:length(csvfolder)){
  csv <- read.csv(csvfolder[file], na.strings=c("NA"))
  
  if (all(!is.na(csv$Seq))){
    total <- rbind(total, csv)
    
   
    
  }
  
  
}

allseqs <- total$Seq[total$Seq != ""]
allseqs<-as.character(allseqs)
charcount <- sum(unname(sapply(allseqs, nchar)))

a <- sum(str_count(allseqs, "A")) / charcount
c <- sum(str_count(allseqs, "C")) / charcount
g <- sum(str_count(allseqs, "G")) / charcount
t <- sum(str_count(allseqs, "T")) / charcount




ntcount <- data.frame(nt=c("A","C","G","T"), props=c(a,c,g,t))



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