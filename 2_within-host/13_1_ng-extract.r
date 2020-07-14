require(ape)
require(stringr)
require(Biostrings)

source("~/vindels/2_within-host/utils.r")



insRandTest <- function(seq, indel, start){
  # this will be the start point of the test insertion 
  smpl <- sample(nchar(seq)+1,500, replace=T)
  
  # for every number in this random sample 
  seq <- sapply(smpl, function(x){insert(seq,indel,x)})
  if(sum(is.na(seq))>0){
    print(seq[which(is.na(seq))])
  }
  
    # generate the result sequence ; use substring to add the insertion sequence into the seqestor 
  aa.seq <- unname(sapply(seq, translate))
  seq.glycs <- unname(sapply(aa.seq,extractGlycs))
  
  png.count <- unname(sapply(seq.glycs,csvcount))
  return(png.count - start)
}

delRandTest <- function(seq, indel, start){
  # this will be the start point of the test insertion 
  smpl <- sample(nchar(indel):nchar(seq),500, replace=T)
  
  # for every number in this random sample 
  seq <- sapply(smpl, function(x){delete(seq,indel,x)})
  if(sum(is.na(seq))>0){
    print(seq[which(is.na(seq))])
  }

  # generate the result sequence ; use substring to add the insertion sequence into the seqestor 
  aa.seq <- unname(sapply(seq, translate))
  seq.glycs <- unname(sapply(aa.seq,extractGlycs))
  
  png.count <- unname(sapply(seq.glycs,csvcount))

  return(png.count - start)
}

# GLYC SITE RANDOMIZATION TEST 

observedGlycChange <- function(anc, indel, pos, option="i"){
  aa.seq <- unname(sapply(anc, translate))
  glycs <- unname(sapply(aa.seq,extractGlycs))
  before <- unname(sapply(glycs,csvcount))
  
  # generate the new sequence by applying the appropriate insertion or deletion
  if (option == "i"){
    newanc <- insert(anc, indel, pos)
  }else{
    newanc <- delete(anc, indel, pos)
  }
  if(is.na(newanc)){
    print(newanc)
    print(indel)
    print(pos)
  }
  
  # recalculate the number of N-glyc sites 
  new.aa <- unname(sapply(newanc, translate))
  glycs <- unname(sapply(new.aa,extractGlycs))
  after <- unname(sapply(glycs,csvcount))

  # report this as a net change in the number of Nglyc sites
  return(after - before)
} 

glycCount <- function(seq){
  # determine the locations of all N-glyc sites in the ancestral sequence 
  aa.seq <- translate(seq)
  # if (is.na(aa.seq)){
  #   return(0)
  # }
  seq.glycs <- extractGlycs(aa.seq)
  csvcount(seq.glycs)
}

#PycharmProjects/hiv-withinhost/

ins <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/all/ins-sep.csv",  sep="\t", stringsAsFactors = F)
del <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/all/del-sep.csv", sep="\t", stringsAsFactors = F)


ins <- ins[,-c(3,4,5)]
del <- del[,-c(3,4,5)]


# apply an adjust to the deletion locations to make them the same as insertions 
#del$pos <- as.numeric(del$pos) + nchar(del$indel)

# Insertions : fill in gaps found in the tip sequences 
ins$tip <- unname(mapply(restoreOtherSeq,ins$tip, ins$anc))

# Deletions : fill in gaps found in the ancestral sequences 
# this is to include any insertions in the ancestor 
del$anc <- unname(mapply(restoreOtherSeq,del$anc, del$tip))

del <- del[-which(nchar(del$indel) > nchar(del$anc)),]
del <- del[-which(del$pos > nchar(del$anc)),]

# Insertions : 
# adds all the other insertions into the ancestor
ins$anc <- unname(mapply(restoreInsAnc, ins$anc, ins$tip, ins$indel, ins$pos))
ins$anc <- gsub("-","",ins$anc)
ins$pos <- ins$pos - nchar(ins$indel) + 1
# not needed for deletions because no sequences contain more than 1 deletion
#del$tip <- unname(mapply(restoreOtherIndels, del$tip, del$anc, del$indel, del$pos))


ins.v <- split(ins, ins$vloop)
del.v <- split(del, del$vloop)


# ----- Randomization Test ---- 
ins.data <- data.frame()
del.data <- data.frame()

for (n in 1:5){
  print(paste0("V-loop ",n))
  iTemp <- ins.v[[n]]
  dTemp <- del.v[[n]]
  
  icounts <- nrow(ins.v[[n]])
  dcounts <- nrow(del.v[[n]])
  
  iTemp$glycs <- unname(sapply(iTemp$anc, glycCount))
  dTemp$glycs <- unname(sapply(dTemp$anc, glycCount))
  
  
  # EXPECTED GLYC CHANGES (RANDOMIZATION TEST)
  # ---------------
  # Insertions
  ires <- t(unname(mapply(insRandTest, iTemp$anc,iTemp$indel, iTemp$glycs)))
  ires <- split(ires, rep(1:nrow(ires), each=ncol(ires)))
  
  iedist <- unname(unlist(lapply(ires, mean)))
  
  iemean <- mean(iedist)
  # Boostraps for expected insertions
  bs.means <- c()
  for (i in 1:100){
    sam <- sample(length(iedist), length(iedist), replace=T)
    iexp.bs <- iedist[sam]
    bs.means[i] <- mean(iexp.bs)
  }
  iequantiles <- quantile(bs.means, c(0.025,0.975))
  
  # Deletions
  dres <- t(unname(mapply(delRandTest, dTemp$anc,dTemp$indel, dTemp$glycs)))
  dres <- split(dres, rep(1:nrow(dres), each=ncol(dres)))

  dedist <- unname(unlist(lapply(dres, mean)))
  
  demean <- mean(dedist)
  # Boostraps for expected deletions
  bs.means <- c()
  for (i in 1:100){
    sam <- sample(length(dedist), length(dedist), replace=T)
    dexp.bs <- dedist[sam]
    bs.means[i] <- mean(dexp.bs)
  }
  dequantiles <- quantile(bs.means, c(0.025,0.975))
  
  
  # OBSERVED GLYCOSYLATION SITE CHANGES (from the data)
  # ----------------------
  iobs <- unname(mapply(observedGlycChange, iTemp$anc, iTemp$indel, iTemp$pos, "i"))
  
  iomean <- mean(iobs)
  # Boostraps for observed insertions
  bs.means <- c()
  for (i in 1:100){
    sam <- sample(length(iobs), length(iobs), replace=T)
    iobs.bs <- iobs[sam]
    bs.means[i] <- mean(iobs.bs)
  }
  ioquantiles <- quantile(bs.means, c(0.025,0.975))
  
  dobs <- unname(mapply(observedGlycChange, dTemp$anc, dTemp$indel, dTemp$pos, "d"))
  
  domean <- mean(dobs)
  # Boostraps for observed deletions
  bs.means <- c()
  for (i in 1:100){
    sam <- sample(length(dobs), length(dobs), replace=T)
    dobs.bs <- dobs[sam]
    bs.means[i] <- mean(dobs.bs)
  }
  doquantiles <- quantile(bs.means, c(0.025,0.975))
  
  ins.data <- rbind(ins.data, data.frame(exp=iemean, 
                                         elower=iequantiles[[1]],
                                         eupper=iequantiles[[2]],
                                         obs=iomean, 
                                         olower=ioquantiles[[1]],
                                         oupper=ioquantiles[[2]],
                                         counts=icounts))
  del.data <- rbind(del.data, data.frame(exp=demean, 
                                         elower=dequantiles[[1]],
                                         eupper=dequantiles[[2]],
                                         obs=domean, 
                                         olower=doquantiles[[1]],
                                         oupper=doquantiles[[2]],
                                         counts=dcounts))
}
del.data[2,1:3] <- del.data[2,1:3] + 0.005
del.data[3,1:3] <- del.data[3,1:3] - 0.005


require(RColorBrewer)
colors <- brewer.pal(5, "Set1")
vloops <- c("V1","V2","V3","V4","V5")
cex=1
par(pty="s", xpd=F, mar=c(7,7,3,1),las=0)
#as.numeric(row.names(data))+20 
# this take in data either as ins.data or del.data

# Deletion data points 
data <- del.data

sizes <- 0.42*sqrt(data$counts)
sizes[3] <- 2.4

lim = c(-0.85,0.5)

plot(data[,c('exp','obs')], pch=1, cex=sizes, lwd=c(10,10,5,10,10), col=colors,xlim=lim,ylim=lim,
     cex.lab=1.85, cex.axis=1.4,cex.main=1.8,las=1, ylab='', xlab='')#main="N-linked Glycosylation Site Changes")
abline(0,1)
title(ylab="Observed Change in PNGS Count\nPer Indel", line=4,cex.lab=1.7)
title(xlab="Expected Change in PNGS Count\nPer Indel", line=5,cex.lab=1.7)
arrows(data[,1], data[,5], data[,1], data[,6], length=0.035, angle=90, code=3)
arrows(data[,2], data[,4], data[,3], data[,4], length=0.035, angle=90, code=3)

# Insertion data points 
data <- ins.data[-3,]
newcol <- colors[-3]
sizes <- 0.42*sqrt(data$counts)

points(data[,c("exp","obs")], pch=21, col=newcol, cex=sizes,lwd=1, bg=newcol )
arrows(data[,1], data[,5], data[,1], data[,6], length=0.035, angle=90, code=3)
arrows(data[,2], data[,4], data[,3], data[,4], length=0.035, angle=90, code=3)

legend(0.25,-0.3,legend=vloops, pch=21,cex=1.5, pt.bg=colors,x.intersp = 1.0,y.intersp=1.0, pt.cex=2.8)
#legend(0.45,0.2,legend=c("Ins", "Del"), pch=c(21,1),cex=1.3, pt.bg=colors[1],col=colors[1], x.intersp = 1.0,y.intersp=1.3, pt.cex=3)
rect(0.25,-0.23,0.46,-0.03)
text(0.4, -0.09, labels="Ins", cex=1.5)
text(0.4, -0.17, labels="Del", cex=1.5)
points(c(0.30,0.30), c(-0.09, -0.17), pch=c(21,1), cex=2.5, lwd=7, col='black', bg='black')
# positions for V1,V2,V3,V4,V5 (top set first )
xposi <- c(-0.42, -0.18,-0.65, -0.03)
yposi <- c(0.48, 0.13, 0.33, 0.07)

xposd <- c(-0.57, -0.12, -0.28, -0.81, -0.31)
yposd <- c(-0.40, -0.26, -0.33,-0.68, -0.06)

text(xposi, yposi, labels=c("V1","V2","V4","V5"),cex=1.2)
text(xposd, yposd, labels=c("V1","V2","V3","V4","V5"),cex=1.2)




# --- out of date ---- 
require(RColorBrewer)
colors <- brewer.pal(5, "Set1")
vloops <- c("V1","V2","V3","V4","V5")
cex=1
par(pty="s", xpd=F, mar=c(6,8,4,1),las=0)
#as.numeric(row.names(data))+20
# this take in data either as ins.data or del.data
sizes <- 0.45*sqrt(data$counts)
sizes[3] <- 1.3
data <- del.data
lim = c(-0.8,0.8)
plot(data[,c('exp','obs')], pch=21, cex=sizes, bg=colors,xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.2,cex.main=1.8, ylab='', xlab='', main="Deletions - PNLGS")
abline(0,1)
title(ylab="Observed Net Change in N-Glyc Sites", line=3,cex.lab=1.3)
title(xlab="Expected Net Change in N-Glyc Sites", line=3,cex.lab=1.3)
arrows(data[,1], data[,5], data[,1], data[,6], length=0.05, angle=90, code=3)
arrows(data[,2], data[,4], data[,3], data[,4], length=0.05, angle=90, code=3)
legend(0.5,-0.2,legend=vloops, pch=21,cex=1.2, pt.bg=colors,x.intersp = 1.0,y.intersp=1.0, pt.cex=2.5)





adjust <- list()
for (i in 1:length(ires)){
  adjust[[i]] <- ires[[i]] - iobs[i]
}

isign <- c()
for (n in 1:length(iobs)){
  idist <- ires[[n]]
  ins.p <- iobs[n]
  
  iQT <- quantile(idist, probs=c(0.025,0.975))
  
  # highlight significant differences 
  if (ins.p < iQT[[1]]){
    isign <- c(isign, "lower")
  }else if(ins.p > iQT[[2]]){
    isign <- c(isign, "higher")
  }else{
    isign <- c(isign, "")
  }
}

cols <- brewer.pal(5,"Set1")
cex=2
par(pty="s", mar=c(6,5,4,1),las=0)

lim = c(0,0.8)

plot(x=i)


vloops <- vector(mode = "list", length = 5)
for (v in 1:5){
  idx <- which(ins$vloop==v)
  for (i in idx){
    vloops[[v]] <- c(vloops[[v]], adjust[[i]])
  }
}


# PLOTTING 

par(mar=c(3,5,2,0))
caxis=1.3
clab=1.4
cmain=1.5

par(mfrow=c(5,1),las=1)
hist(vloops[[1]], col="red", prob=T, xaxt="n",breaks=seq(min(vloops[[4]])-0.5,max(vloops[[4]])+0.5),main="V1")
lines(density(vloops[[1]],bw=0.6),lwd=2)

hist(vloops[[2]], col="red", prob=T, xaxt="n", breaks=seq(min(vloops[[4]])-0.5,max(vloops[[4]])+0.5),main="V2")
lines(density(vloops[[2]],bw=0.6),lwd=2)

hist(vloops[[3]], col="red", prob=T,xaxt="n", breaks=seq(min(vloops[[4]])-0.5,max(vloops[[4]])+0.5),main="V3")
lines(density(vloops[[3]],bw=0.6),lwd=2)

hist(vloops[[4]], col="red", prob=T,xaxt="n",breaks=seq(min(vloops[[4]])-0.5,max(vloops[[4]])+0.5),main= "V4")
lines(density(vloops[[4]],bw=0.6),lwd=2)

hist(vloops[[5]], col="red", prob=T, breaks=seq(min(vloops[[4]])-0.5,max(vloops[[4]])+0.5),main="V5")
lines(density(vloops[[5]],bw=0.6),lwd=2)




ins[which(isign=="higher"),]


dres <- t(unname(mapply(randomizationTest, del$Tip,del$indel)))
dres <- split(dres, rep(1:nrow(dres), each=ncol(dres)))

iobs <- unname(sapply(ins$anc, glycCount))

isign <- c()
for (n in 1:length(iobs)){
  idist <- res[[n]]
  
  iQT <- quantile(idist, probs=c(0.025,0.975))
  
  ins.p <- observed[n]
  
  # highlight significant differences 
  if (ins.p < iQT[[1]]){
    isign <- c(isign, "lower")
  }else if(ins.p > iQT[[2]]){
    isign <- c(isign, "higher")
  }else{
    isign <- c(isign, "")
  }
  
}




#headers <- c("accno", "vloop", "indel", "pos", "tip","anc", "patient")
#colnames(ins) <- headers
#colnames(del) <- headers

new.ins <- ins[,c(1,3,4,5)]
new.del <- del[,c(1,3,4,5)]

# new.ins$anc <- unname(mapply(insAlign, ins$seq, ins$pos, ins$anc, ins$tip))
# new.ins$tip <- ins$tip
# 
# new.del$anc <- del$anc
# new.del$tip <- unname(mapply(delAlign, del$seq, del$pos, del$anc, del$tip))

new.ins$tip <- sapply(ins$tip, translate)
new.del$tip <- sapply(del$tip, translate)

new.ins$anc <- sapply(ins$anc, translate)
new.del$anc <- sapply(del$anc, translate)

new.ins$tip.glycs <- unlist(sapply(new.ins$tip, extractGlycs))
new.del$tip.glycs <- unlist(sapply(new.del$tip, extractGlycs))

new.ins$anc.glycs <- unlist(sapply(new.ins$anc, extractGlycs))
new.del$anc.glycs <- unlist(sapply(new.del$anc, extractGlycs))




write.table(new.ins, paste0(path,"13_nglycs/ins-edit.csv"), sep="\t", quote=F, row.names=F)
write.table(new.del, paste0(path,"13_nglycs/del-edit.csv"), sep="\t", quote=F, row.names=F)



ins$aaseq <- NULL
del$aaseq <- NULL

ins$original <- mapply(insOriginal, indel=ins$indel, pos=ins$pos, vseq=ins$tip)
del$original <- mapply(delOriginal, indel=ins$indel, pos=ins$pos, vseq=ins$tip)
# for both: 
# determine the start and stop of all nglycs 
# collect them in one column separated by "-", comma separated


# insertions 
# determine what the original sequence is 
# start = pos - nchar(insertion sequence)
# end = pos 

# original seq: insertion removed  
# original <- paste0(substring(vloop, 0, start), substring(vloop, end, nchar(vloop)))


# deletions 
# start = pos - nchar(insertion sequence)
# end = pos 

# original seq: deletion added back in 
# original <- paste0(substring(vloop,0,start), del, substring(vloop,end,nchar(vloop)))

# --------------- OUTDATED

# insAlign <- function(indels, pos, anc, seq){
#   i.list <- str_split(indels, ",")[[1]]
#   p.list <- str_split(pos, ",")[[1]]
#   
#   for (idx in 1:length(i.list)){
#     
#     len <- nchar(i.list[idx])
#     ix <- i.list[idx]
#     px <- as.numeric(p.list[idx])
#     
#     anc <- paste0(substr(anc, 0, px-len), paste(rep("-", len),collapse=""), substr(anc,px-len+1, nchar(anc)))
#   }
#   
#   anc
# }


# 
# delAlign <- function(indels, pos, anc, seq){
#   i.list <- str_split(indels, ",")[[1]]
#   p.list <- str_split(pos, ",")[[1]]
#   p.list <- as.numeric(p.list)
#   
#   for (idx in 1:length(i.list)){
#     len <- nchar(i.list[idx])
#     ix <- i.list[idx]
#     px <- p.list[idx]
#     
#     seq <- paste0(substr(seq, 0, px), paste(rep("-", len),collapse=""), substr(seq,px+1, nchar(seq)))
#     p.list[(idx+1):length(p.list)] <- p.list[(idx+1):length(p.list)] + len
#   }
#   
#   seq
# }

