ifolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/10_10_Indels_i", full.names=TRUE)
tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/10_9_Indels_total", full.names=TRUE)

totals.df <- data.frame(subtype=character(), stringsAsFactors = F)
intf.df <- data.frame(subtype=character(), stringsAsFactors = F)
subtypes <- c()
for (n in 1:length(ifolder)){
  subtype <- strsplit(strsplit(tfolder[n], "/")[[1]][7],"\\+")[[1]][1]
  subtypes <- c(subtypes, subtype)
  
  
  totals <- read.csv(tfolder[n], sep=":", header=F,  stringsAsFactors = F, colClasses = "character")
  interfered <- read.csv(ifolder[n], sep=":", header=F, stringsAsFactors = F, colClasses="character")
  interfered[is.na(interfered)]<-""
  totals$vloop <- sapply(totals$V1, function(x) strsplit(x,"\\.")[[1]][3])
  totals$subtype <- rep(subtype, nrow(totals))
  totals$count <- mapply(function(x, y) { length(strsplit(x, ",")[[1]]) + 
      length(strsplit(y, ",")[[1]])}, totals$V2, totals$V3)
  
  
  interfered$vloop <- sapply(interfered$V1, function(x) strsplit(x,"\\.")[[1]][3])
  interfered$subtype <- rep(subtype, nrow(interfered))
  interfered$count <- mapply(function(x, y) { length(strsplit(x, ",")[[1]]) + 
      length(strsplit(y, ",")[[1]])}, interfered$V2, interfered$V3)
  
  

  
  totals.df <- rbind(totals.df, totals)
  intf.df <- rbind(intf.df, interfered)
}

#MAIN 3 DATA FRAMES
disrupt <- data.frame(stringsAsFactors = F)
prop.df <- data.frame(stringsAsFactors = F)

#REGRESSION FIGURE --- Proportion of NGS sites vs Proportion of Disrupted NGS Sites
#contains the proportions of amino acids belonging to NGSites
ngProps <- read.csv("~/PycharmProjects/hiv-evolution-master/ngCounts.csv") #generated using 10_9a Nglycs
ngProps.df <- data.frame()


for (i in subtypes){
  for (m in 1:5){
    tCount <- sum(totals.df[which(totals.df$subtype==i & totals.df$vloop==m),6])
    iCount <- sum(intf.df[which(intf.df$subtype==i & intf.df$vloop==m),6])
    
    ngProps.df <- rbind(ngProps.df, data.frame(subtype=i, vloop=m, prop=ngProps[which(ngProps$subtype==i), paste0('V',m)]))
    
    
    if (tCount != 0){
      disrupt <- rbind(disrupt, data.frame(subtype=i, vloop=m, count=iCount,total=tCount, prop=iCount/tCount))
    
    }else{
      disrupt <- rbind(disrupt, data.frame(subtype=i, vloop=m, count=iCount,total=tCount, prop=0))
    }
    
  }
}


#ngProps.df <- ngProps.df[-33,]
#disrupt <- disrupt[-33,]


#FIGURE-----
require(RColorBrewer)
colors <- rep(brewer.pal(5, 'Set1'),7)
#colors <- colors[-33]
lim = c(0,0.75)
par(pty="s",xpd=F,mar=c(8,9,4,4))
plot(x=ngProps.df$prop, y=disrupt$prop,  las=1,pch=(disrupt$vloop+20),bg=colors,cex=3, cex.axis=1.5,main=NULL, xlab="",ylab="",xlim=lim,ylim=lim) 
par(xpd=NA)
legend(0.85,0.5,legend=c('V1 ','V2 ','V3 ','V4 ','V5'), pch=c(21,22,23,24,25), cex=1.6, pt.bg=colors,x.intersp = 1.5,y.intersp=1.2, pt.cex=3)
par(xpd=F)
#abline(0,1)
title(ylab="Proportion of Indels affecting PNGS", line=4.5, cex.lab=1.7)
title(xlab="Proportion of PNGS Content", line=3.5, cex.lab=1.7)
text(0,0.04,label="F1",cex=1.5)
text(0.09,0.18,label="AE",cex=1.5)
abline(0,1,lty=5,lwd=1.5)

#abline(lm(disrupt$prop ~ ngProps.df$prop))



#HEAT MAP FIGURE 
#contains the total counts of NG sites in all sequences 
# format ---- 
ngTotal <- read.csv("~/PycharmProjects/hiv-evolution-master/ngTotal.csv", stringsAsFactors = F)
ngTotal$accno <- NULL
v1 <- ngTotal[which(ngTotal$vloop==1),c(1,3)]
v2 <- ngTotal[which(ngTotal$vloop==2),c(1,3)]
v3 <- ngTotal[which(ngTotal$vloop==3),c(1,3)]
v4 <- ngTotal[which(ngTotal$vloop==4),c(1,3)]
v5 <- ngTotal[which(ngTotal$vloop==5),c(1,3)]

#stats analysis between subtypes
summary(glm(count ~ subtype, data = v5, family=poisson))

v1.df <- as.data.frame.matrix(table(v1))
v2.df <- as.data.frame.matrix(table(v2))
v3.df <- as.data.frame.matrix(table(v3))
v4.df <- as.data.frame.matrix(table(v4))
v5.df <- as.data.frame.matrix(table(v5))


v1.prop <- v1.df/rowSums(v1.df)
v2.prop <- v2.df/rowSums(v2.df)
v3.prop <- v3.df/rowSums(v3.df)
v4.prop <- v4.df/rowSums(v4.df)
v5.prop <- v5.df/rowSums(v5.df)


clades <- c("AE","AG","A1","B","C","D","F1")
sizes <- c(nrow(ngTotal[which(ngTotal$subtype=="01_AE" & ngTotal$vloop==1),]),
           nrow(ngTotal[which(ngTotal$subtype=="02_AG" & ngTotal$vloop==1),]),
           nrow(ngTotal[which(ngTotal$subtype=="A1" & ngTotal$vloop==1),]),
           nrow(ngTotal[which(ngTotal$subtype=="B" & ngTotal$vloop==1),]),
           nrow(ngTotal[which(ngTotal$subtype=="C" & ngTotal$vloop==1),]),
           nrow(ngTotal[which(ngTotal$subtype=="D" & ngTotal$vloop==1),]),
           nrow(ngTotal[which(ngTotal$subtype=="F1" & ngTotal$vloop==1),]))


#NGS Plot #2
par(mar=c(8,4,2,2), pty="m")
plot(1:8,c(0,2,4,5,7,8,8,8),cex=0, ylim=c(-0.5,8.5),xlim=c(1,35), xaxt='n',yaxt='n', 
     main=NULL,xlab="", ylab="", cex.lab=1.7,cex.main=2, bty='n')
title(xlab="V-Loop / Clade", line=4, cex.lab=1.5)
title(ylab="Number of N-linked Glycosylation Sites", line=1.7, cex.lab=1.5)

par(las=1)
axis(1,at=c(1:35),lab=rep(clades,5), line=-1.8, tick=F,cex.axis=0.7)
axis(1,at=c(4,11,18,25,32), lab=c("V1","V2","V3","V4","V5"),line=0,tick=F,cex.axis=1.2)
par(las=2)
axis(2, at=0:8,pos=0.5,cex.axis=1.2)
par(xpd=NA)

#x axis vertical line 
segments(0.5,8.5,0.5,-1.4,lwd=2.5)


#bright color scheme
#colors <- c(1,0,0,0,0.5,1,0,1,0,1,0,0.8,1,0.5,0)
#moderate color scheme
require(RColorBrewer)
colors <- brewer.pal(5, 'Set1')
require(plotrix)
require(scales)
full.df <- data.frame(stringsAsFactors = F)
for (n in 1:5){
  temp <- as.data.frame.matrix(table(ngTotal[which(ngTotal$vloop==n),c(1,3)]))
  prop <- temp/rowSums(temp)
  bool <- temp < 0.10*rowSums(temp)
  col <- colors[((n*3)-2):(n*3)]
  r <- col[1]
  g <- col[2]
  b <- col[3]
  
  #vertical line segments
  segments((n*7)+0.5,8.5,(n*7)+0.5,-1.4,lwd=3)
  
  
  #iterate through all possible positions in the 7x9 grid, check for the proportion and draw the color density with it
  for (x in 1:7){
    for (y in 1:9){
      if (!is.null(prop[x,y])){
        alv <- prop[x,y]
        
      }else{
        alv <- 0.00
      }
      #boxes for the whole plot
      rect(((n-1)*7)+x-0.5, y-1.5,((n-1)*7)+x+0.5,y-0.5, col=alpha(colors[n],alv),lwd=0.8)
      
      #for low density boxes      
      if(!is.null(prop[x,y]) && isTRUE(bool[x,y]) && temp[x,y] > 0){
        text(((n-1)*7)+x,y-1,labels=temp[x,y],col=colors[n], cex=0.8) #colors[n] #rgb(r,g,b)
      }
      
    }
  }
}

#significance
higher <- data.frame(x=c(1,11,18,19,23,24,25),y=c(3,2,1,1,3,3,3))
for (i in 1:7){
  text(higher[i,1],(higher[i,2]-0.2),labels="*",cex=1.7)
  arrows(higher[i,1],higher[i,2],higher[i,1],(higher[i,2]+0.25), length=0.07, lwd=1.5)
}
lower <- data.frame(x=c(33,28),y=c(1,3))
for (i in 1:2){
  text(lower[i,1],(lower[i,2]+0.2),labels="*",cex=1.7)
  arrows(lower[i,1],lower[i,2],lower[i,1],(lower[i,2]-0.25), length=0.07, lwd=1.5)
}

#color legend bottom left hand corner
rect(0.5,-2.8,1.5,-1.8, col=alpha(colors[1],0.25),lwd=0.8)
rect(2.5,-2.8,3.5,-1.8, col=alpha(colors[1],0.5),lwd=0.8)
rect(4.5,-2.8,5.5,-1.8, col=alpha(colors[1],0.75),lwd=0.8)
rect(6.5,-2.8,7.5,-1.8, col=alpha(colors[1],1),lwd=0.8)
text(1,-3, labels="25%")
text(3,-3, labels="50%")
text(5,-3, labels="75%")
text(7,-3, labels="100%")
text(30, -1.9, labels="Available Sequences")
for (x in 1:7){
  text(23+(x*1.8),-2.3, labels=clades[x])
  text(23+(x*1.8),-2.8, labels=sizes[x])
}







