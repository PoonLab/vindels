require(ape)
require(usedist)
args <- commandArgs(trailingOnly = T)
folder <- args[1]
outfolder <- args[2]

logfiles <- Sys.glob(paste0(folder,"*.log"))
treefiles <- Sys.glob(paste0(folder, "trees/*"))

getDict <- function(infile){
  seqList <- list()
  lines <- readLines(infile)
  #print(head(lines))
  for (l in lines){
    search <- grepl("\\d+\\s'.+'", l)
    #print(search)
    # this will load the translate dictionary (when search has been found)
    if (search){
      fields <- strsplit(l, " ")[[1]]
      num <- gsub("\t","",fields[1])
      header <- gsub("[',]","",fields[2])
      seqList[[num]] <- header
    }
  }
  return(seqList)
}

convertTree <- function(str, dict){
  require(ape)
  tree <- strsplit(str, " ")[[1]][6]
  newtree <- gsub("\\[[^]]+\\]","", tree)
  newtree <- read.tree(text=newtree)
  
  newtree$tip.label <- sapply(newtree$tip.label, function(n){dict[[n]]})
  return(newtree)
}

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

  #m.state <- sliced[which.max(sliced$posterior),"state"]
  
  tfile <- file(treefiles[k])
  
  seqDict <- getDict(tfile)
  #print(length(seqDict))
  #print(head(seqDict))
  
  trees <- readLines(tfile)
  trees <- trees[grepl("STATE",trees)]
  
  tlen <- length(trees) - 1
  interval <- c(tlen*0.1+2,tlen+1)
  trees <- trees[interval[1]:interval[2]]
  
  mat <- matrix(ncol=length(seqDict)*2-2, nrow=length(trees))
  for (j in 1:length(trees)){
    newtree <- convertTree(trees[j], seqDict)
    mat[j,] <- newtree$edge.length
  }
  
  # calculate the euclidean distance ; find the centroid
  dmat <- dist(mat, method="euclidean")
  cent <- dist_to_centroids(dmat, rep(1,length(trees)))
  
  choice <- which.min(cent[,3])
  
  finaltree <- convertTree(trees[choice], seqDict)
  write.tree(finaltree, file=paste0(outfolder, strsplit(filename,"\\.")[[1]][1], ".tree"))

}


