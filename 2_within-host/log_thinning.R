
args <- commandArgs(trailingOnly = T)
if (length(args) != 2 ){
  print("USAGE: Rscript log_thinning [input dir] [output dir]")
  quit()
}
#path = "~/PycharmProjects/hiv-withinhost/6_hm_constant/unfinished/"
for (i in 1:length(args)){
  if (!endsWith(args[i], "/")){
    args[i] <-  paste0(args[i], "/")
  }
}

logfiles <- Sys.glob(paste0(args[1],"*.log"))
outpath <- args[2]

for (file in logfiles){
  log <- read.table(file, header=T,comment.char = '#')
  name <- strsplit(basename(file),"\\.")[[1]][1]
  
  loglen <- nrow(log)-1
  interval <- c(((loglen)*0.1)+1,loglen) + 1
  log2 <- log[seq(interval[1], interval[2], length.out=901), ]
  
  print(name)
  #print(head(log2))
  write.table(log2, file=paste0(outpath, name,'-thin.log'), sep='\t', quote=F, row.names=F)

}
