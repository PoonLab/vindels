require(ape)
require(stringr)
require(Biostrings)

source("~/vindels/2_within-host/utils.r")



insRandTest <- function(seq, indel, start){
  # this will be the start point of the test insertion 
  smpl <- sample(nchar(seq)+1,1, replace=T)
  
  # for every number in this random sample 
  aa.seq <- sapply(smpl, function(x){translate(insert(seq,indel,x))})
  if(sum(is.na(aa.seq))>0){
    print(aa.seq[which(is.na(aa.seq))])
  }
  
    # generate the result sequence ; use substring to add the insertion sequence into the seqestor 
  #aa.seq <- unname(sapply(seq, translate))
  #seq.glycs <- unname(sapply(aa.seq,extractGlycs))
  png.count <- unname(sapply(aa.seq,function(x) {csvcount(extractGlycs(x))}))
  
  #png.count <- unname(sapply(seq.glycs,csvcount))
  return(png.count - start)
}

delRandTest <- function(seq, indel, start){
  # this will be the start point of the test insertion 
  smpl <- sample(nchar(indel):nchar(seq),1, replace=T)
  
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

observedGlycChange <- function(anc, indel, pos, before, option="i",n){

  # generate the new sequence by applying the appropriate insertion or deletion
  if (option == "i"){
    newanc <- insert(anc, indel, pos)
  }else{
    newanc <- delete(anc, indel, pos)
  }
  if(is.na(newanc) || newanc == ""){
    print(newanc)
    print(indel)
    print(pos)
  }
  
  # recalculate the number of N-glyc sites 
  new.aa <- translate(newanc)
  after <- csvcount(extractGlycs(new.aa))

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

bootstrap <- function(sampl, reps, n=1000){
  replicate(n, mean(sample(sampl, length(sampl)/reps, replace=T)))
}

#PycharmProjects/hiv-withinhost/

ins <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/all/ins-current.csv",  sep="\t", stringsAsFactors = F)
del <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/all/del-current.csv", sep="\t", stringsAsFactors = F)

ins = ins[!grepl("[^ACGT]", ins$indel),]
del = del[!grepl("[^ACGT]", del$indel),]

ins <- ins[,-c(4,5,6)]
del <- del[,-c(4,5,6)]

del <- del[!(del$pos - nchar(del$indel) < 0),]


# apply an adjust to the deletion locations to make them the same as insertions 
#del$pos <- as.numeric(del$pos) + nchar(del$indel)

# Insertions : fill in gaps found in the tip sequences 
ins$tip <- unname(mapply(restoreAllGaps,ins$tip, ins$anc))

# Deletions : fill in gaps found in the ancestral sequences 
# this is to include any insertions in the ancestor 
del$anc <- unname(mapply(restoreAllGaps,del$anc, del$tip))


# Insertions : 
# adds all the other insertions into the ancestor
ins$anc <- unname(mapply(restoreOtherIndels, ins$anc, ins$tip, ins$indel, ins$pos))
ins$anc <- gsub("-","",ins$anc)
ins$pos <- ins$pos - nchar(ins$indel) + 1
# not needed for deletions because no sequences contain more than 1 deletion
del$tip <- unname(mapply(restoreOtherIndels, del$tip, del$anc, del$indel, del$pos))

ins$glycs <- unname(sapply(ins$anc, glycCount))
del$glycs <- unname(sapply(del$anc, glycCount))

any(grepl("-",ins$anc))
any(grepl("-",del$anc))


ins.v <- split(ins, ins$vloop)
del.v <- split(del, del$vloop)


# ----- Randomization Test ---- 
ins.data <- data.frame()
del.data <- data.frame()

ins.rep <- matrix(nrow=1000, ncol=4)
del.rep <- matrix(nrow=1000, ncol=4)


for (n in 1:5){
  print(n)
  iTemp <- ins.v[[n]]
  dTemp <- del.v[[n]]
  
  icounts <- nrow(ins.v[[n]])
  dcounts <- nrow(del.v[[n]])
  
  
  # EXPECTED GLYC CHANGES 
  # ---------------
  # Insertions
  print("EXPECTED INSERTIONS...")
  ires <- unname(mapply(insRandTest, iTemp$anc,iTemp$indel, iTemp$glycs))
  #ires <- split(ires, rep(1:nrow(ires), each=ncol(ires)))
  
  iedist <- ires
  
  iemean <- mean(iedist)
  # Boostraps for expected insertions
  bs.means <- bootstrap(iedist, 200)
  bs.q <- quantile(bs.means - iemean, c(0.025,0.975))
  
  ieconint <- c()
  ieconint[1] <- iemean - bs.q[[2]]
  ieconint[2] <- iemean - bs.q[[1]]
  
  # Deletions
  print("EXPECTED DELETIONS...")
  dres <- unname(mapply(delRandTest, dTemp$anc,dTemp$indel, dTemp$glycs))
  #dres <- split(dres, rep(1:nrow(dres), each=ncol(dres)))

  dedist <- dres
  
  demean <- mean(dedist)
  
  # Boostraps for expected deletions
  bs.means <- bootstrap(dedist, 200)
  bs.q <- quantile(bs.means - demean, c(0.025,0.975))
  
  deconint <- c()
  deconint[1] <- demean - bs.q[[2]]
  deconint[2] <- demean - bs.q[[1]]
  
  
  # OBSERVED GLYCOSYLATION SITE CHANGES (from the data)
  # ----------------------
  print("OBSERVED INSERTIONS...")
  iobs <- unname(mapply(observedGlycChange, iTemp$anc, iTemp$indel, iTemp$pos, iTemp$glycs, "i", 1:nrow(iTemp)))
  
  iomean <- mean(iobs)
  # Boostraps for observed insertions
  bs.means <- bootstrap(iobs, 200)
  bs.q <- quantile(bs.means - iomean, c(0.025,0.975))
  
  ioconint <- c()
  ioconint[1] <- iomean - bs.q[[2]]
  ioconint[2] <- iomean - bs.q[[1]]
  
  print("OBSERVED DELETIONS...")
  dobs <- unname(mapply(observedGlycChange, dTemp$anc, dTemp$indel, dTemp$pos, dTemp$glycs, "d", 1:nrow(dTemp)))
  
  domean <- mean(dobs)
  # Boostraps for observed deletions
  bs.means <- bootstrap(dobs, 200)
  bs.q <- quantile(bs.means - domean, c(0.025,0.975))
  
  doconint <- c()
  doconint[1] <- domean - bs.q[[2]]
  doconint[2] <- domean - bs.q[[1]]
  
  
  
  iQT <- quantile(ires, c(0.025, 0.975))
  if (iomean < iQT[[1]]){
    isign <- "lower"
  }else if(iomean > iQT[[2]]){
    isign <- "higher"
  }else{
    isign <- ""
  }
  
  #dnull <- apply(dres, 2, mean)
  dQT <- quantile(dres, c(0.025, 0.975))
  if (domean < dQT[[1]]){
    dsign <- "lower"
  }else if(domean > dQT[[2]]){
    dsign <- "higher"
  }else{
    dsign <- ""
  }

  ins.data <- rbind(ins.data, data.frame(exp=iemean, elower=ieconint[1],eupper=ieconint[2],
                                         obs=iomean, olower=ioconint[1], oupper=ioconint[2],
                                         counts=icounts, sign=isign))
  del.data <- rbind(del.data, data.frame(exp=demean, elower=deconint[1],eupper=deconint[2],
                                         obs=domean, olower=doconint[1],oupper=doconint[2],
                                         counts=dcounts, sign=dsign))
  
  ioreps <- unname(sapply(split(iobs, ins[ins$vloop==n, 'rep']), mean))
  iereps <- unname(sapply(split(iedist, ins[ins$vloop==n, 'rep']), mean))
  
  doreps <- unname(sapply(split(dobs, del[del$vloop==n, 'rep']), mean))
  dereps <- unname(sapply(split(dedist, del[del$vloop==n, 'rep']), mean))
  
  idx <- ((n-1)*200+1):(n*200)
  ins.rep[idx, 1] <- rep(n,200)
  ins.rep[idx, 2] <- 1:200
  ins.rep[idx, 3] <- iereps
  ins.rep[idx, 4] <- ioreps
  
  del.rep[idx, 1] <- rep(n,200)
  del.rep[idx, 2] <- 1:200
  del.rep[idx, 3] <- dereps
  del.rep[idx, 4] <- doreps
  
}



require(RColorBrewer)
require(scales)
colors <- brewer.pal(5, "Set1")
vloops <- c("V1","V2","V3","V4","V5")
cex=1
par(pty="s", xpd=F, mar=c(7,7,3,1),las=0)

# this take in data either as ins.data or del.data
data <- del.data
d.rep <- del.rep
v3offset <- 0
#data[3,1:6] <- data[3,1:6] + v3offset
#d.rep[41:60,3:4] <- d.rep[41:60,3:4] + v3offset

sizes <- 0.06*sqrt(data$counts)
sizes[3] <- 2.4

lim = c(-0.85,0.5)

plot(data[,c('exp','obs')], pch=1, cex=sizes, lwd=c(10,10,5,10,10), col=alpha(colors,0.75),xlim=lim,ylim=lim,
     cex.lab=1.85, cex.axis=1.4,cex.main=1.8,las=1, ylab='', xlab='')#main="N-linked Glycosylation Site Changes")

# CLOUDS 
points(d.rep[,3], d.rep[,4], pch=1, col=alpha(colors[d.rep[,1]],0.4), cex=0.7, lwd=1.4)

abline(0,1)
title(ylab="Observed Change in PNGS Count\nPer Indel", line=4,cex.lab=1.7)
title(xlab="Expected Change in PNGS Count\nPer Indel", line=5,cex.lab=1.7)
#arrows(data[,1], data[,5], data[,1], data[,6], length=0.035, angle=90, code=3)
#arrows(data[,2], data[,4], data[,3], data[,4], length=0.035, angle=90, code=3)

# Insertion data points 
i <- c(1,2,0,3,4)
data <- ins.data[-3,]
newcol <- colors[-3]
sizes <- 0.06*sqrt(data$counts)
i.rep <- ins.rep[-c(401:600),]

# MEANS 
points(data[,c("exp","obs")], pch=21, col=newcol, cex=sizes,lwd=1, bg=alpha(newcol,0.55))


# CLOUDS 
points(i.rep[,3], i.rep[,4], pch=21, cex=0.7, lwd=1.4,bg=alpha(newcol[i[i.rep[,1]]],0.4),col=alpha(newcol[i[i.rep[,1]]],0.3))

# CONFIDENCE INTERVALS 
#arrows(data[,1], data[,5], data[,1], data[,6], length=0.035, angle=90, code=3)
#arrows(data[,2], data[,4], data[,3], data[,4], length=0.035, angle=90, code=3)



legend(0.25,-0.26,legend=vloops, pch=21,cex=1.5, pt.bg=colors,x.intersp = 1.5,y.intersp=1.2, pt.cex=3.5)
#legend(0.45,0.2,legend=c("Ins", "Del"), pch=c(21,1),cex=1.3, pt.bg=colors[1],col=colors[1], x.intersp = 1.0,y.intersp=1.3, pt.cex=3)
rect(0.25,-0.23,0.453,-0.03)
text(0.39, -0.09, labels="Ins", cex=1.5)
text(0.39, -0.17, labels="Del", cex=1.5)
points(c(0.30,0.30), c(-0.09, -0.17), pch=c(21,1), cex=2.5, lwd=7, col='black', bg='black')

# positions for V1,V2,V3,V4,V5 (top set first )
xposi <- c(-0.41, -0.28,-0.61, 0.0)
yposi <- c(0.53, 0.11, 0.42, 0.08)

xposd <- c(-0.52, -0.26, -0.05, -0.86, -0.32)
yposd <- c(-0.42, -0.41, -0.32,-0.67, -0.09)

text(xposi, yposi, labels=c("V1","V2","V4","V5"),cex=1.2)
text(xposd, yposd, labels=c("V1","V2","V3","V4","V5"),cex=1.2)







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

