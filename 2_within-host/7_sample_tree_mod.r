# Tree modifications 
# used for modifying different attributes of a tree and returning the result

require(ape)

args <- commandArgs(trailingOnly = T)

if (length(args) != 3){
  print("USAGE: Rscript 7_sample_tree_mod.r [log dir] [input dir] [output dir]")
  quit()
}

for (i in 1:length(args)){
  if (!endsWith(args[i],"/")){
    args[i] <- paste0(args[i],"/") 
  }
}

# input directory of sampled BEAST trees
# relies on the presence of a "prelim" folder being present
logpath <- args[1]
tfolder <- args[2]
outfolder <- args[3]


dir.create(outfolder, showWarnings = F)
logfolder <- Sys.glob(paste0(logpath,"*.log"))

for (logfile in logfolder){
  
  print(logfile)
  filename <- basename(logfile)
  pat <- strsplit(filename, "\\.")[[1]][1]
  
  # uses log file name to find and read BEAST log file
  logfile <- read.csv(logfile, sep="\t", skip=4)
  
  pat_trees <- Sys.glob(paste0(tfolder,pat,"*"))
  for (tree in pat_trees){
    treename <- basename(tree)
    
    treename <- gsub("-original","",treename)
    
    intree <- read.tree(tree)
    
    state <- as.numeric(paste0(strsplit(strsplit(treename, "\\.")[[1]][1], "_")[[1]][3],"00000"))
    
    rescale.factor <- logfile[which(logfile$state == state),"ucld.mean"]

    print(rescale.factor)
  
    # rescales all the edge lengths of the tree
    intree$edge.length <- (intree$edge.length * rescale.factor)

    # writes the rescaled tree to a new folder called "rescaled"
    write.tree(intree,paste0(outfolder, gsub(".sample","",treename)))
  }
}

