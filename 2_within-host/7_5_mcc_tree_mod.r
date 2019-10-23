# Tree modifications 
# used for modifying different attributes of a tree and returning the result

require(ape)

args <- commandArgs(trailingOnly = T)

# input directory of sampled BEAST trees
# relies on the presence of a "prelim" folder being present
if (length(args) != 2){
  quit(status="USAGE: Rscript 7_5_mcc_tree_mod.r [working directory] [log file directory]")
}
for (i in 1:length(args)){
  if (!endsWith(args[i],"/")){
    args[i] <- paste0(args[i],"/") 
  }
}

path <- args[1]
logpath <- args[2]

if (!dir.exists(paste0(path,"prelim/"))){
  quit(status="USAGE: Rscript 7_5_mcc_tree_mod.r [working directory] [log file directory]")
}

dir.create(paste0(path,"rescaled/"), showWarnings = FALSE)
infolder <- Sys.glob(paste0(path,"prelim/*.tree"))

# for (f in folders){
#   trees <- Sys.glob(paste0(f, "/*.tree.sample"))
#   foldername <- paste0(basename(f),"/") 
#   
#   dir.create(paste0(path,"rescaled_multi/",foldername), showWarnings = FALSE)
#   
count <- 0
r.vec2 <- c()
for (treefile in infolder){
  count <- count + 1
  print(treefile)
  
  intree <- read.nexus(treefile)
  filename <- basename(treefile)

  # edits the tree file name to get the BEAST log file name
  logname <- paste0(strsplit(filename, "\\.")[[1]][1], ".log")

  # uses log file name to find and read BEAST log file
  logfile <- read.csv(paste0(logpath,logname), sep="\t", skip=3)

  print(logname)

  #counts the number of MCMC steps
  loglen <- nrow(logfile) -1
  print(loglen)
  
  if (loglen != 10000){
    print(status="INCOMPLETE FILE")
    next()
  }
  
  # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
  interval <- c(loglen*0.1+1,loglen+1)
  print(interval)

  # calculates the rescale factor using the median of the UCLD.MEAN column (can check that this matches UCLD.MEDIAN on tracer)
  means <- logfile$ucld.mean[interval[1]:interval[2]]

  rescale.factor <- median(means)
  print(rescale.factor)
  
  r.vec[count] <- rescale.factor

  # rescales all the edge lengths of the tree
  intree$edge.length <- (intree$edge.length * rescale.factor)

  # writes the rescaled tree to a new folder called "rescaled"
  #write.tree(intree,paste0(outpath, filename))
}
#}

