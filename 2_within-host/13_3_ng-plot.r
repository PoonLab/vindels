
path <- "~/Lio/"
path <- "~/PycharmProjects/hiv-withinhost/"
ins <- read.csv(paste0(path,"13_nglycs/interfered/insertions.csv"), sep="\t",stringsAsFactors = F)
del <- read.csv(paste0(path,"13_nglycs/interfered/deletions.csv"),sep="\t", stringsAsFactors = F)

t.ins <- read.csv(paste0(path,"10_nucleotide/ins-sep-only.csv"),row.names=1, stringsAsFactors = F)
t.del <- read.csv(paste0(path,"10_nucleotide/del-sep-only.csv"),row.names=1, stringsAsFactors = F)

iseqs <- read.csv(paste0(path,"13_nglycs/ins-edit.csv"),row.names=1,stringsAsFactors = F)
dseqs <- read.csv(paste0(path,"13_nglycs/del-edit.csv"),row.names=1,stringsAsFactors = F)


ngprops <- read.csv(paste0(path,"13_nglycs/interfered/ngprops.csv"), row.names=1)
require(stringr)
require(RColorBrewer)
splitcsv <- function(str){
  fields <- strsplit(str, ",")[[1]]
  as.numeric(fields)
}

count <- function(str){
  if (grepl(",", str)){
    str_count(str,",") + 1
  }else{
    1
  }
}


# total insertion counts 
itotals <- c()
dtotals <- c()
for (i in 1:5){
  itotals[i] <- nrow(t.ins[t.ins$Vloop==i,])
  dtotals[i] <- nrow(t.del[t.del$Vloop==i,])
}
# retrieved from the nrow totals of 10_nt_proportions
ins$vloop <- sapply(ins$header, function(x){as.numeric(strsplit(x, "_")[[1]][4])})
ins$count <- sapply(ins$pos,count)
ins_overlap <- sapply(c(1:5),function(x){sum(ins[ins$vloop==x,"count"])} )

ins_props <- ins_overlap / itotals

# total deletion counts
del$vloop <- sapply(del$header, function(x){as.numeric(strsplit(x, "_")[[1]][4])})
del$count <- sapply(del$pos,count)

del_overlap <- sapply(c(1:5),function(x){sum(del[del$vloop==x,"count"])} )

del_props <- del_overlap / dtotals




# RANDOMIZATION TEST 
# -----------------------------------------
sampleString <- function(len, vloop){
  len <- len-1
  idx <- sample(1:(nchar(vloop)-len),100, replace=TRUE)
  strings <- sapply(idx, function(x){substr(vloop, x, x+len)})
  #a <- unname(sapply(strings, function(x){str_count(x, "A")/nchar(x)}))
  #c <- unname(sapply(strings, function(x){str_count(x, "C")/nchar(x)}))
  #g <- unname(sapply(strings, function(x){str_count(x, "G")/nchar(x)}))
  #t <- unname(sapply(strings, function(x){str_count(x, "T")/nchar(x)}))
  return(strings)
}




iSample <- list(c(),c(),c(),c())
dSample <- list(c(),c(),c(),c())
# generates the randomly sampled substrings for each indel
for (row in 1:nrow(ins)){
  itemp <- sampleString(total.ins[row,"len"], total.ins[row,"Vseq"])
  for (i in 1:4){
    iSample[[i]] <- c(iSample[[i]], itemp[[i]])
  }
}
for (row in 1:nrow(total.del)){
  dtemp <- sampleString(total.del[row,"len"], total.del[row,"Vseq"])
  for (i in 1:4){
    dSample[[i]] <- c(dSample[[i]], dtemp[[i]])
  }
}

# compares the observed proportion to the overall distribution of each nucleotide 
isign <- c()
dsign <- c()
for (i in 1:4){
  idist <- iSample[[i]]
  ddist <- dSample[[i]]
  
  iQT <- quantile(idist, probs=c(0.025,0.975))
  dQT <- quantile(ddist, probs=c(0.025,0.975))
  
  ins.p <- ins.nt[i,2]
  del.p <- del.nt[i,2]
  
  # highlight significant differences 
  if (ins.p < iQT[[1]]){
    isign <- c(isign, "lower")
  }else if(ins.p > iQT[[2]]){
    isign <- c(isign, "higher")
  }else{
    isign <- c(isign, "")
  }
  
  # highlight significant differences 
  if (del.p < dQT[[1]]){
    dsign <- c(dsign, "lower")
  }else if(del.p > dQT[[2]]){
    dsign <- c(dsign, "higher")
  }else{
    dsign <- c(dsign, "")
  }
  
}

ins.nt$sign <- isign
del.nt$sign <- dsign


cols <- brewer.pal(5,"Set1")
cex=2
par(pty="s", mar=c(6,5,4,1),las=0)

lim = c(0,0.8)
plot(y=ins_props, x=unlist(ngprops[1,]), pch=c(21:25), bg=cols,xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2,cex=3.5, main="Insertions", xlab="N-glyc Content", ylab="Proportion of \nInsertions Generating PNGS")
abline(0,1)
legend(0.6,0.37,legend=c('V1 ','V2 ','V3 ','V4 ','V5'), pch=c(21:25), cex=1.5, pt.bg=cols,x.intersp = 1.3,y.intersp=1.2, pt.cex=3)

plot(y=del_props, x=unlist(ngprops[2,]), pch=c(21:25), bg=cols,xlim=lim,ylim=lim,
     cex.lab=1.3, cex.axis=1.3,cex.main=2.2,cex=3.5, main="Deletions", xlab="N-glyc Content", ylab="Proportion of \nDeletions Removing PNGS")
abline(0,1)
legend(0.6,0.37,legend=c('V1 ','V2 ','V3 ','V4 ','V5'), pch=c(21:25), cex=1.5, pt.bg=cols,x.intersp = 1.3,y.intersp=1.2, pt.cex=3)
