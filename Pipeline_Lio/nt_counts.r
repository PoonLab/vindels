#Two csv files containing the numbers of indel and non-indel nucleotides in all 35 variable loops 
noindel <- read.csv("~/PycharmProjects/hiv-evolution-master/noindel.csv")
indel <- read.csv("~/PycharmProjects/hiv-evolution-master/indel.csv",stringsAsFactors = F)

require(MASS)
require(RColorBrewer)

#data frame to hold the chi square statistical test values
data.df <- data.frame(filename=indel$filename,stringsAsFactors = F)

#data frame to hold nucleotide counts + proportions for all subtypes + vloops
big.df <- data.frame()
toTest <- data.frame()
for (i in 1:nrow(indel)){
  
  subtype <- toString(indel[i,1])
  
  for (x in 1:5){
    
    in.df <- indel[i,c(1+x, 6+x, 11+x, 16+x)]
    out.df <- noindel[i, c(1+x, 6+x, 11+x, 16+x)]
    
    #prep for chisq test 
    
    toTest <- rbind(in.df,out.df)
    toTest <- matrix(unlist(toTest), nrow=2,ncol=4)
    
    #build the big.df
    combined <- rbind(in.df,in.df,in.df,in.df,out.df,out.df,out.df,out.df)
    names(combined) <- c('A', 'C', 'G', 'T')
    combined$total <- rowSums(combined)
    combined$nt <- rep(c('A','C','G','T'), 2)
    combined$subtype <- c(rep(subtype,8))
    combined$vregion <- c(rep(x,8))
    combined$indel <- c(rep("inside",4), rep("outside",4))
    
    big.df <- rbind(big.df, combined)
    if (x == 3 && subtype == "F1"){
      next
    }
    #df$indel <- c('indel', 'no_indel')
    rownames(toTest) <- c('inside', 'outside')
    colnames(toTest) <- c('A', 'C', 'G', 'T')
    chi <- chisq.test(toTest)
    data.df[i, paste0('V',x)] <- chi$p.value
    
    
    
  }
}

ntides <- strsplit('ACGT',"")[[1]]

for (n in ntides){
  big.df[which(big.df$nt == n),'freq'] <- big.df[which(big.df$nt == n) , n] 
  big.df[which(big.df$nt == n),'prop'] = big.df[which(big.df$nt == n), 'freq'] / big.df[which(big.df$nt == n), 'total'] 
}


#DATA ANALYSIS PART 1 ------------------------------------
#for pooling the subtypes together to test if the vregions are significantly different   
#nucleotide counts are summed across all the subtypes
require(MASS)
vr.df <- cbind(big.df[1:40, 'freq'],big.df[41:80, 'freq'],big.df[81:120, 'freq'],big.df[121:160, 'freq'],big.df[161:200, 'freq'],big.df[201:240, 'freq'],big.df[241:280, 'freq'] )
vr.df <- data.frame(vr.df)
vr.df$sum <- rowSums(vr.df)

vrcount <- array(data=vr.df$sum, dim=c(4,2,5), dimnames= list("nucleotide"=c(ntides), "indel"=c('inside','outside'), "vregion"=c(1,2,3,4,5)))
vrcount2 <- as.data.frame(as.table(vrcount))
vrFit <- glm(Freq ~ nucleotide + indel + vregion, data = vrcount2, family=poisson)

vrFit2 <- glm(Freq ~ nucleotide * indel * vregion, data = vrcount2, family=poisson)
#step AIC function -- vregions 
vrFit.aic <- stepAIC(vrFit, scope=list(upper=vrFit2, lower=~1), direction='both', trace=TRUE)


#DATA ANALYSIS PART 2 -------------------------------------
#for pooling the vregions together to test if the subtypes are significantly different
st.df <- data.frame()
for (i in 1:7){
  for (j in 1:5){
    for (k in 1:8){
      st.df[(i-1)*8 + k,j] <- big.df[((i-1)*40)+((j-1)*8)+k, 'freq']
    }
  }
}
st.df$sum <- rowSums(st.df)

stcount <- array(data=st.df$sum, dim=c(4,2,7), dimnames=list("nucleotide"=c(ntides), "indel"= c('inside','outside'), "subtype" =c("01_AE","02_AG", "A1","B","C","D","F1")))
stcount2 <- as.data.frame(as.table(stcount))
stFit <- glm(Freq ~ nucleotide + indel + subtype, data = stcount2, family=poisson)

stFit2 <- glm(Freq ~ nucleotide * indel * subtype, data = stcount2, family=poisson)
stFit.aic <- stepAIC(stFit, scope=list(upper=stFit2, lower=~1), direction='both', trace=TRUE)

#4 way table
ntcount <- array(data=big.df$freq, dim=c(4,2,5,7), dimnames=list("nucleotide"=c(ntides), "indel"= c('inside','outside'), "vregion"=c(1,2,3,4,5), "subtype" =c("01_AE","02_AG", "A1","B","C","D","F1")))

ntcount2 <- as.data.frame(as.table(ntcount))
ntfit <- glm(Freq ~ nucleotide + indel + vregion + subtype, data=ntcount2, family=poisson)
ntfit2 <- glm(Freq ~ nucleotide * indel * vregion * subtype, data=ntcount2, family=poisson)

ntfit.aic <- stepAIC(ntfit, scope=list(upper=ntfit2, lower=~1, direction='both', trace=TRUE))


#SCALE SIZES -----------------------------------------
sizes.df <- data.frame()
sizes.v <- c()
pvalues.df <- data.frame()
for (n in 1:nrow(data.df)){
  sizes <- c()
  pvalues <- c()
  for (i in 1:5){
    toAdd <- data.df[n,(i+1)]
    
    if (is.na(toAdd) || is.nan(toAdd) || toAdd >= (0.05/35)){
      #toAdd <- 1
      sizes <- c(sizes, 1.1)
    }else{
      #toAdd <- -log(toAdd, base=10)
      sizes <- c(sizes, 2.2)
    }
    pvalues <- c(pvalues, toAdd)
  }
  print(sizes)
  sizes.df <- rbind(sizes.df, sizes)
  pvalues.df <- rbind(pvalues.df, pvalues)
  sizes.v <- c(sizes.v,sizes)
}



colors <- brewer.pal(5, 'Set1')
nt <- "G"
nt <- "T"
nt <- "C"
nt <- "A"

inside <- data.frame(big.df[which(big.df$indel == 'inside' & big.df$nt==nt),])
outside <- data.frame(big.df[which(big.df$indel == 'outside' & big.df$nt==nt) ,])

list.df <- list()

for (nt in ntides){
  inside <- data.frame(big.df[which(big.df$indel == 'inside' & big.df$nt==nt),])
  outside <- data.frame(big.df[which(big.df$indel == 'outside' & big.df$nt==nt) ,])
  
  df <- data.frame(outside=outside$prop, inside=inside$prop, nt=inside$nt, vregion=inside$vregion)
  list.df[[nt]] <- df
}


#A
cex=2
par(pty="s", mfrow=c(2,2), xpd=NA, mar=c(3,8,4,1),las=0)

lim = c(0.24,0.46)
plot(list.df[[1]][,1:2], cex=sizes.v, pch=(list.df[[1]]$vregion+20), bg=rep(colors, 7),xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=0.9,cex.main=1.5, ylab='', xlab='')
text(0.187,0.475,labels="a)", cex=1.5)
text(0.245,0.452,labels="A", cex=1.5)
title(ylab="Proportion Within Indels", line=2.2,cex.lab=1)
title(xlab="Proportion Outside Indels", line=2.2,cex.lab=1)
par(xpd=F)
abline(0,1)

#C
lim=c(0.05,0.30)
par(mar=c(3,1,4,8))
plot(list.df[[2]][,1:2], cex=sizes.v, pch=(list.df[[2]]$vregion+20), bg=rep(colors, 7),xlim=lim,ylim=lim,
     cex.lab=1.2, cex.axis=0.9,cex.main=1.5, ylab='', xlab='')
abline(0,1)
title(ylab="Proportion Within Indels", line=2.2,cex.lab=1)
title(xlab="Proportion Outside Indels", line=2.2,cex.lab=1)

par(xpd=NA)
text(-0.010,0.317,labels="b)", cex=1.5)
text(0.056,0.29,labels="C", cex=1.5)
par(xpd=F)
#G
lim=c(0.12,0.40)
par(mar=c(4,8,3,1))
plot(list.df[[3]][,1:2], cex=sizes.v, pch=(list.df[[3]]$vregion+20), bg=rep(colors, 7),xlim=lim,ylim=lim,
     cex.lab=1, cex.axis=0.9,cex.main=1.5, ylab='', xlab='')
abline(0,1)
title(ylab="Proportion Within Indels", line=2.2,cex.lab=1)
title(xlab="Proportion Outside Indels", line=2.2,cex.lab=1)

par(xpd=NA)
text(0.053,0.42,labels="c)", cex=1.5)
text(0.127,0.387,labels="G", cex=1.5)
par(xpd=F)
legend(0.33,0.26,legend=c('V1  ','V2  ','V3  ','V4  ','V5  '), pch=c(21,22,23,24,25),cex=1.2, pt.bg=colors,x.intersp = 1.4,y.intersp=1.2, pt.cex=2.2)

#T
lim=c(0.11,0.33)
par(mar=c(4,1,3,8))
plot(list.df[[4]][,1:2], cex=sizes.v, pch=(list.df[[4]]$vregion+20), bg=rep(colors, 7),xlim=lim,ylim=lim,
     cex.lab=1.2, cex.axis=0.9,cex.main=1.5, ylab='', xlab='')
abline(0,1)
title(ylab="Proportion Within Indels", line=2.2,cex.lab=1)
title(xlab="Proportion Outside Indels", line=2.2,cex.lab=1)
par(xpd=NA)
text(0.058,0.345,labels="d)", cex=1.5)
text(0.115,0.32,labels="T", cex=1.5)
par(xpd=F)


