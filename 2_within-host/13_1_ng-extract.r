require(ape)
require(stringr)
require(Biostrings)

source("~/vindels/2_within-host/utils.r")



insRandTest <- function(seq, indel, start){
  # this will be the start point of the test insertion 
  smpl <- sample(nchar(seq)+1,300, replace=T)
  
  # for every number in this random sample 
  seq <- sapply(smpl, function(x){insert(seq,indel,x)})
  
    # generate the result sequence ; use substring to add the insertion sequence into the seqestor 
  aa.seq <- unname(sapply(seq, translate))
  seq.glycs <- unname(sapply(aa.seq,extractGlycs))
  
  png.count <- unname(sapply(seq.glycs,csvcount))
  return(png.count - start)
}

delRandTest <- function(seq, indel, start){
  # this will be the start point of the test insertion 
  smpl <- sample((nchar(seq) - nchar(indel) + 1),300, replace=T)
  
  # for every number in this random sample 
  seq <- sapply(smpl, function(x){delete(seq,indel,x)})

  # generate the result sequence ; use substring to add the insertion sequence into the seqestor 
  aa.seq <- unname(sapply(seq, translate))
  seq.glycs <- unname(sapply(aa.seq,extractGlycs))
  
  png.count <- unname(sapply(seq.glycs,csvcount))

  return(png.count - start)
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
path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
ins <- read.csv(paste0(path, "13_nglycs/ins-sep.csv"),  sep="\t", stringsAsFactors = F)
del <- read.csv(paste0(path,"13_nglycs/del-sep.csv"), sep="\t", stringsAsFactors = F)

ins <- ins[-c(which(ins$Pos ==0)),]

ins <- ins[,-c(3,4)]
del <- del[,-c(3,4)]

ins$Vpos <- NULL
del$Vpos <- NULL

# apply an adjust to the deletion locations to make them the same as insertions 
del$Pos <- as.numeric(del$Pos) + nchar(del$Seq)

# Insertions : fill in gaps found in the tip sequences 
res <- as.data.frame(t(unname(mapply(restoreTipDel,ins$Vseq, ins$Anc, ins$Seq,ins$Pos))))
ins$Vseq <- as.character(res[,1])
ins$Pos <- as.numeric(as.character(res[,2]))


# Deletions : fill in gaps found in the ancestral sequences 
del$Anc <- unname(mapply(restoreAncIns, del$Anc, del$Vseq, del$Seq))

# Insertions : 
ins$Anc <- unname(mapply(removeOtherGaps, ins$Anc,ins$Vseq, ins$Seq, ins$Pos))
ins$Anc <- gsub("-","",ins$Anc)
# not needed for deletions because no sequences contain more than 1 deletion
#ins$Anc <- unname(mapply(removeOtherGaps, del$Vseq,del$Anc, del$Seq, del$Pos))


ins.v <- split(ins, ins$Vloop)
del.v <- split(del, del$Vloop)


# GLYC SITE RANDOMIZATION TEST 

observedGlycChange <- function(anc, indel, pos, option="i"){
  anc <- gsub("-","",anc)
  aa.seq <- unname(sapply(anc, translate))
  glycs <- unname(sapply(aa.seq,extractGlycs))
  before <- unname(sapply(glycs,csvcount))
  
  if (option == "i"){
    newanc <- insert(anc, indel, (pos - nchar(indel) + 1))
  }else{
    newanc <- delete(anc, indel, (pos - nchar(indel) + 1))
  }
  
  #print(newanc)
  # generate the result sequence ; use substring to add the insertion sequence into the seqestor 
  
  new.aa <- unname(sapply(newanc, translate))
  glycs <- unname(sapply(new.aa,extractGlycs))
  after <- unname(sapply(glycs,csvcount))
  # recalculate the number of N-glyc sites 
  # adjust the location of all N-glyc locations falling AFTER the random position to check their similarity 
  
  # report this as a positive or a negative result zzzz
  return(after - before)
} 

ins.data <- data.frame()
del.data <- data.frame()

for (n in 1:5){
  print(paste0("V-loop ",n))
  iTemp <- ins.v[[n]]
  dTemp <- del.v[[n]]
  
  icounts <- nrow(ins.v[[n]])
  dcounts <- nrow(del.v[[n]])
  
  iTemp$glycs <- unname(sapply(iTemp$Anc, glycCount))
  dTemp$glycs <- unname(sapply(dTemp$Anc, glycCount))
  
  
  # EXPECTED GLYC CHANGES (RANDOMIZATION TEST)
  # ---------------
  # Insertions
  ires <- t(unname(mapply(insRandTest, iTemp$Anc,iTemp$Seq, iTemp$glycs)))
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
  dres <- t(unname(mapply(delRandTest, dTemp$Anc,dTemp$Seq, dTemp$glycs)))
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
  iobs <- unname(mapply(observedGlycChange, iTemp$Anc, iTemp$Seq, iTemp$Pos, "i"))
  
  iomean <- mean(iobs)
  # Boostraps for observed insertions
  bs.means <- c()
  for (i in 1:100){
    sam <- sample(length(iobs), length(iobs), replace=T)
    iobs.bs <- iobs[sam]
    bs.means[i] <- mean(iobs.bs)
  }
  ioquantiles <- quantile(bs.means, c(0.025,0.975))
  
  dobs <- unname(mapply(observedGlycChange, dTemp$Anc, dTemp$Seq, dTemp$Pos, "d"))
  
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


require(RColorBrewer)
colors <- brewer.pal(5, "Set1")
vloops <- c("V1","V2","V3","V4","V5")
cex=1
par(pty="s", xpd=F, mar=c(6,8,4,1),las=0)
#as.numeric(row.names(data))+20
# this take in data either as ins.data or del.data

# Deletion data points 
data <- del.data
sizes <- 0.4*sqrt(data$counts)
sizes[3] <- 1.2

lim = c(-0.8,0.8)

plot(data[,c('exp','obs')], pch=1, cex=sizes, lwd=10, col=colors,xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.2,cex.main=1.8, ylab='', xlab='', main="Insertions - PNLGS")
abline(0,1)
title(ylab="Observed Net Change in N-Glyc Sites", line=3,cex.lab=1.3)
title(xlab="Expected Net Change in N-Glyc Sites", line=3,cex.lab=1.3)
arrows(data[,1], data[,5], data[,1], data[,6], length=0.05, angle=90, code=3)
arrows(data[,2], data[,4], data[,3], data[,4], length=0.05, angle=90, code=3)

# Insertion data points 
data <- ins.data
sizes <- 0.5*sqrt(data$counts)
sizes[3] <- 1.2
points(data[,c("exp","obs")], pch=21, col=colors, cex=sizes,lwd=1, bg=colors )
arrows(data[,1], data[,5], data[,1], data[,6], length=0.05, angle=90, code=3)
arrows(data[,2], data[,4], data[,3], data[,4], length=0.05, angle=90, code=3)

legend(0.45,-0.2,legend=vloops, pch=21,cex=1.3, pt.bg=colors,x.intersp = 1.0,y.intersp=1.0, pt.cex=2.5)
legend(0.45,0.2,legend=c("Ins", "Del"), pch=c(21,1),cex=1.3, pt.bg=colors[1],col=colors[1], x.intersp = 1.0,y.intersp=1.3, pt.cex=3)

# positions for V1,V2,V3,V4,V5 (top set first )
xposi <- c(-0.44, -0.18, -0.2, -0.73, 0)
yposi <- c(0.35, 0.3, 0.05, 0.35, 0.30)

xposd <- c(-0.55, -0.4, -0.1, -0.75, -0.1)
yposd <- c(-0.30, -0.2, -0.22,-0.35, -0.03)

text(xposi, yposi, labels=c("V1","V2","V3","V4","V5"))
text(xposd, yposd, labels=c("V1","V2","V3","V4","V5"))

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
  idx <- which(ins$Vloop==v)
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


dres <- t(unname(mapply(randomizationTest, del$Tip,del$Seq)))
dres <- split(dres, rep(1:nrow(dres), each=ncol(dres)))

iobs <- unname(sapply(ins$Anc, glycCount))

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

new.ins$tip <- sapply(ins$Vseq, translate)
new.del$tip <- sapply(del$Vseq, translate)

new.ins$anc <- sapply(ins$Anc, translate)
new.del$anc <- sapply(del$Anc, translate)

new.ins$tip.glycs <- unlist(sapply(new.ins$tip, extractGlycs))
new.del$tip.glycs <- unlist(sapply(new.del$tip, extractGlycs))

new.ins$anc.glycs <- unlist(sapply(new.ins$anc, extractGlycs))
new.del$anc.glycs <- unlist(sapply(new.del$anc, extractGlycs))




write.table(new.ins, paste0(path,"13_nglycs/ins-edit.csv"), sep="\t", quote=F, row.names=F)
write.table(new.del, paste0(path,"13_nglycs/del-edit.csv"), sep="\t", quote=F, row.names=F)



ins$aaseq <- NULL
del$aaseq <- NULL

ins$original <- mapply(insOriginal, indel=ins$Seq, pos=ins$Pos, vseq=ins$Vseq)
del$original <- mapply(delOriginal, indel=ins$Seq, pos=ins$Pos, vseq=ins$Vseq)
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

