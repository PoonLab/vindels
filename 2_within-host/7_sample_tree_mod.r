# Tree modifications 
# used for modifying different attributes of a tree and returning the result

require(ape)

args <- commandArgs(trailingOnly = T)

if (length(args) != 2){
  print("USAGE: Rscript 7_sample_tree_mod.r [working directory] [log directory]")
  quit()
}

for (i in 1:length(args)){
  if (!endsWith(args[i],"/")){
    args[i] <- paste0(args[i],"/") 
  }
}

# input directory of sampled BEAST trees
# relies on the presence of a "prelim" folder being present
infolder <- args[1]
logpath <- args[2]

dir.create(paste0(infolder,"final/"), showWarnings = F)
treefolder <- Sys.glob(paste0(infolder,"prelim/*.tree.sample"))

for (treefile in treefolder){
  
  print(treefile)
  
  intree <- read.tree(treefile)
  filename <- basename(treefile)
  
  state <- as.numeric(paste0(strsplit(strsplit(filename, "\\.")[[1]][1], "_")[[1]][3],"00000"))

  # edits the tree file name to get the BEAST log file name
  logname <- paste0(strsplit(filename, "_")[[1]][1], ".log")

  # uses log file name to find and read BEAST log file
  logfile <- read.csv(paste0(logpath,logname), sep="\t", skip=4)

  print(logname)
  
  rescale.factor <- logfile[which(logfile$state == state),"ucld.mean"]
  #counts the number of MCMC steps
  #loglen <- nrow(logfile) -1
  #print(loglen)
  # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
  #interval <- c(loglen*0.1+1,loglen+1)
  #print(interval)
  # calculates the rescale factor using the median of the UCLD.MEAN column (can check that this matches UCLD.MEDIAN on tracer)
  #rescale.factor <- median(logfile$ucld.mean[interval[1]:interval[2]])
  print(rescale.factor)

  # rescales all the edge lengths of the tree
  intree$edge.length <- (intree$edge.length * rescale.factor)

  # writes the rescaled tree to a new folder called "rescaled"
  write.tree(intree,paste0("~/PycharmProjects/hiv-withinhost/7SampleTrees/final/", gsub(".sample","",filename)))
}
#}

