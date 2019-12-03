# Tree modifications 
# used for modifying different attributes of a tree and returning the result

require(ape)

args <- commandArgs(trailingOnly = T)

# input directory of sampled BEAST trees
# relies on the presence of a "prelim" folder being present
if (length(args) != 1){
  print("USAGE: Rscript 7_5_mcc_tree_mod.r [working directory]")
  quit()
}
for (i in 1:length(args)){
  if (!endsWith(args[i],"/")){
    args[i] <- paste0(args[i],"/") 
  }
}

treefolder <- args[1]

if (!dir.exists(paste0(treefolder,"prelim/"))){
  quit(status="USAGE: Rscript 7_5_mcc_tree_mod.r [working directory]")
}



dir.create(paste0(treefolder,"final/"), showWarnings = F)
treefiles <- Sys.glob(paste0(treefolder,"prelim/*.tree"))


r.vec <- c()
for (i in 1:length(treefiles)){
  tre <- read.nexus(treefiles[i])
  
  print("INPUT")
  print(treefiles[i])
  
  
  inlog <- strsplit(basename(treefiles[i]),"\\.")[[1]][1]
  # uses log file name to find and read BEAST log file
  logfile <- read.csv(paste0("~/PycharmProjects/hiv-withinhost/6_hm/final/",inlog,".log"), sep="\t", skip=4)
  
  print("LOGFILE")
  print(inlog)
  
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
  
  r.vec[i] <- rescale.factor

  # rescales all the edge lengths of the tree
  tre$edge.length <- (tre$edge.length * rescale.factor)

  # writes the rescaled tree to a new folder called "final"
  write.tree(tre,paste0(treefolder,"final/",treefiles[i]))
}
#}

