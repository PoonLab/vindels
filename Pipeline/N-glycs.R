ifolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/10_8_NGlycs_i", full.names=TRUE)
tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/10_7_NGlycs_total", full.names=TRUE)

count.df <- data.frame(subtype=character(), stringsAsFactors = F)
intf.df <- data.frame(subtype=character(), stringsAsFactors = F)
subtypes <- c()
for (n in 1:length(ifolder)){
  subtype <- strsplit(strsplit(tfolder[n], "/")[[1]][7],"\\+")[[1]][1]
  subtypes <- c(subtypes, subtype)
  
  
  totals <- read.csv(tfolder[n], sep=":", header=F,stringsAsFactors = F)
  interfered <- read.csv(ifolder[n], sep=":", header=F, stringsAsFactors = F)
  
  totals$vloop <- sapply(totals$V1, function(x) strsplit(x,"\\.")[[1]][3])
  totals$subtype <- rep(subtype, nrow(totals))
  totals$count <- mapply(function(x, y) { length(strsplit(x, ",")[[1]]) + 
      length(strsplit(y, ",")[[1]])}, totals$V2, totals$V3)
  
  
  interfered$vloop <- sapply(interfered$V1, function(x) strsplit(x,"\\.")[[1]][3])
  interfered$subtype <- rep(subtype, nrow(interfered))
  interfered$count <- mapply(function(x, y) { length(strsplit(x, ",")[[1]]) + 
      length(strsplit(y, ",")[[1]])}, interfered$V2, interfered$V3)
  
  

  
  count.df <- rbind(count.df, totals)
  intf.df <- rbind(intf.df, interfered)
}

#MAIN 3 DATA FRAMES
totals.df <- data.frame(stringsAsFactors = F)
disrupt <- data.frame(stringsAsFactors = F)
prop.df <- data.frame(stringsAsFactors = F)

#REGRESSION FIGURE --- Proportion of NGS sites vs Proportion of Disrupted NGS Sites
#contains the proportions of amino acids belonging to NGSites
ngProps <- read.csv("~/PycharmProjects/hiv-evolution-master/ngCounts.csv")
ngProps.df <- data.frame()


for (i in subtypes){
  for (m in 1:5){
    tCount <- sum(count.df[which(count.df$subtype==i & count.df$vloop==m),6])
    iCount <- sum(intf.df[which(intf.df$subtype==i & intf.df$vloop==m),6])
    
    ngProps.df <- rbind(ngProps.df, data.frame(subtype=i, vloop=m, prop=ngProps[which(ngProps$subtype==i), paste0('V',m)]))
    
    
    if (tCount != 0){
      disrupt <- rbind(disrupt, data.frame(subtype=i, vloop=m, count=iCount,total=tCount, prop=iCount/tCount))
    
    }else{
      disrupt <- rbind(disrupt, data.frame(subtype=i, vloop=m, count=iCount,total=tCount, prop=0))
    }
    
  }
}

#FIGURE-----
require(RColorBrewer)
colors <- brewer.pal(5, 'Set1')
lim = c(0,0.55)
par(pty="m",xpd=F,mar=c(6,7,2,2))
plot(x=ngProps.df$prop, y=disrupt$prop, panel.first = grid(NULL,NULL,lty=1), las=1,pch=(disrupt$vloop+20),bg=rep(colors, 7),cex=2.8, cex.lab=1.3,main=NULL, xlab="",ylab="") 
legend(0.003,0.285,legend=c('V1  ','V2  ','V3  ','V4  ','V5  '), pch=c(21,22,23,24,25), cex=1.6, pt.bg=colors,x.intersp = 1.5,y.intersp=1.6, pt.cex=2.4)
#abline(0,1)
title(ylab="Disrupted PNLGS / Total PNLGS", line=4, cex.lab=1.5)
title(xlab="PNLGS Amino Acids / Total Amino Acids", line=3, cex.lab=1.5)
text(0,0.01,label="F1")
text(0.09,0.125,label="01_AE")


#abline(lm(disrupt$prop ~ ngProps.df$prop))



#HEAT MAP FIGURE 
#contains the total counts of NG sites in all sequences 
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

#NGS Plot #2
par(mar=c(6,4,2,2), pty="m")
plot(1:8,c(0,2,4,5,7,8,8,9),cex=0, ylim=c(-0.5,8.5),xlim=c(1,35), xaxt='n',yaxt='n', 
     main=NULL,xlab="", ylab="", cex.lab=1.7,cex.main=2, bty='n')
title(xlab="V-Loop / Subtype", line=4, cex.lab=1.5)
title(ylab="Number of N-linked Glycosylation Sites", line=1.7, cex.lab=1.5)

par(las=1)
axis(1,at=c(1:35),lab=rep(c("AE","AG","A1","B","C","D","F1"),5), line=-1.8, tick=F,cex.axis=0.7)
axis(1,at=c(4,11,18,25,32), lab=c("V1","V2","V3","V4","V5"),line=1,tick=F,cex.axis=1.2)
par(las=2)
axis(2, at=0:7,pos=0.5,cex.axis=1.2)
par(xpd=NA)

#first vertical segment
segments(0.5,8.5,0.5,-1.4,lwd=2.5)

#bright color scheme
#colors <- c(1,0,0,0,0.5,1,0,1,0,1,0,0.8,1,0.5,0)
#moderate color scheme
require(RColorBrewer)
colors <- brewer.pal(5, 'Set1')
require(plotrix)
require(scales)
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
      if(!is.null(prop[x,y]) && isTRUE(bool[x,y]) && temp[x,y] > 1){
        text(((n-1)*7)+x,y-1,labels=temp[x,y],col=colors[n], cex=0.8) #colors[n] #rgb(r,g,b)
      }
      
    }
  }
}
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







