ins <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/interfered/insertions.csv", sep="\t",stringsAsFactors = F)
del <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/interfered/deletions.csv",sep="\t", stringsAsFactors = F)
ngprops <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/interfered/ngprops.csv", row.names=1)
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
itotals <- c(44,23,1,15,29)      # retrieved from the nrow totals of 10_nt_proportions
ins$vloop <- sapply(ins$header, function(x){as.numeric(strsplit(x, "_")[[1]][4])})
ins$count <- sapply(ins$pos,count)
ins_overlap <- sapply(c(1:5),function(x){sum(ins[ins$vloop==x,"count"])} )

ins_props <- ins_overlap / itotals

# total deletion counts
dtotals <- c(84,45,6,56,50)
del$vloop <- sapply(del$header, function(x){as.numeric(strsplit(x, "_")[[1]][4])})
del$count <- sapply(del$pos,count)

del_overlap <- sapply(c(1:5),function(x){sum(del[del$vloop==x,"count"])} )

del_props <- del_overlap / dtotals


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
