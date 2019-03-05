
logfiles <- Sys.glob("~/PycharmProjects/hiv-withinhost/6BEASTout/*.log")
setwd("~/PycharmProjects/hiv-withinhost/6BEASTout/thinned_logs/")

args <- commandArgs(TRUE)


for (file in logfiles){
  log <- read.table(file, header=T,comment.char = '#')
  name <- strsplit(basename(file),"\\.")[[1]][1]
  print("hello")

  log2 <- log[seq(1, nrow(log), length.out=1001), ]

  write.table(log2, file=paste0(name,'-thin.log'), sep='\t', quote=F, row.names=F)

}
