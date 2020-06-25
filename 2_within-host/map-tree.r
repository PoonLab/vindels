require(ape)
args <- commandArgs(trailingOnly = T)
folder <- args[1]
outfolder <- args[2]

logfiles <- Sys.glob(paste0(folder,"*.log"))
treefiles <- Sys.glob(paste0(folder, "trees/*"))

for (k in 1:length(logfiles)){
  filename <- basename(logfiles[k])
  print(filename)
  
  # DIRECTORY NUMBER 1 PROCESSING
  data <- read.table(logfiles[k], sep="\t", skip=4, head=T)
  loglen <- nrow(data) -1
  #print(loglen)
  # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
  interval <- c(loglen*0.1+1,loglen+1)
  sliced <- data[interval[1]:interval[2],]
  sliced <- sliced[as.numeric(sliced$state) %% 100000 == 0, ]

  m.state <- sliced[which.max(sliced$posterior),"state"]
  
  infile <- file(treefiles[k])
  trees <- readLines(infile)
  #print(tail(trees))
  print(paste0("STATE_",m.state))
  write(strsplit(trees[grepl(paste0("STATE_",m.state),trees)]," ")[[1]][6],
             file=paste0(outfolder,filename))
}


