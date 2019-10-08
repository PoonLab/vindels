path = "~/PycharmProjects/hiv-withinhost/6_hm_constant/unfinished/"
logfiles <- Sys.glob(paste0(path,"*.log"))


for (file in logfiles){
  log <- read.table(file, header=T,comment.char = '#')
  name <- strsplit(basename(file),"\\.")[[1]][1]
  
  loglen <- nrow(log)-1
  interval <- c(((loglen-1)*0.1)+1,loglen)
  log2 <- log[seq(interval[1], interval[2], length.out=901), ]
  
  print(name)
  #print(head(log2))
  write.table(log2, file=paste0(path, name,'-thin.log'), sep='\t', quote=F, row.names=F)

}
