args <- commandArgs(trailingOnly = T)
if (!endsWith(args[1],"/")) args[1] <- paste0(args[1],"/")
source("~/vindels/2_within-host/utils.r")

require(ape)
require(stringr)
require(phangorn)
require(data.table)
require(bbmle)


# INSERTION PARSING ----------
path <- "~/work/"
ifolder <- Sys.glob(paste0(path,"indels/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"indels/del/*.tsv"))
sep <- "\t"
trefolder <- paste0(path,"trees/original/")

newreg <- "30647"

ifolder <- ifolder[grepl(newreg,ifolder)]
dfolder <- dfolder[grepl(newreg,dfolder)]

iint <- list()
itip <- list()
dint <- list()
dtip <- list()
dict <- list()

for (file in 1:length(ifolder)){
  filename <- basename(ifolder[file])
  print(file)
  full.id <- gsub("_\\d+\\.tsv$","",filename)
  
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F, sep=sep)
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F, sep=sep)
  
  iCSV$count <- unname(sapply(iCSV$indel, csvcount, delim=":"))
  dCSV$count <- unname(sapply(dCSV$indel, csvcount, delim=":"))
  
  iCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(iCSV))
  dCSV$pat <- rep(strsplit(filename, "\\.")[[1]][1], nrow(dCSV))
  
  iCSV$header <- unname(mapply(labels, iCSV$header, iCSV$pat))
  dCSV$header <- unname(mapply(labels, dCSV$header, dCSV$pat))
  
  # reads in the tree
  treename <- strsplit(filename, "\\.")[[1]][1]
  tre <- read.tree(paste0(paste0(trefolder,treename,".tree.sample")))

  res <- unname(sapply(iCSV$header, findAncestor, tree=tre)) 
  iCSV$index <- res
  if (file == 1){
    idx <- seq(1,1320, 5)
    dict <- lapply(1:length(unique(res)), function(x){
      name <- iCSV[idx[x], "header"]
    })  
  }
  
  iCSV$length <- tre$edge.length[match(res, tre$edge[,2])]
  dCSV$length <- iCSV$length
  
  lens <- node.depth.edgelength(tre)
  # midpoint = (rtt length of tip) + (rtt length of ancestor) / 2
  iCSV$rtt.mid <- (lens[res] + lens[tre$edge[match(res, tre$edge[,2]),1]]) / 2
  dCSV$rtt.mid <- iCSV$rtt.mid
  
  # ----- TIP + INTERIOR INDEL COUNTS ---
  # Cumulative data frame split by interior vs tip
  id <- gsub("[ab]_","",full.id)   # 16362-100
  
  # this regexp matches TIP SEQUENCES (NOT CONTAINING LEFT AND RIGHT BRACKETS)
  tips <-  which(grepl("^[^\\(\\):\n]+$", iCSV$header))
  nodes <- which(!grepl("^[^\\(\\):\n]+$", iCSV$header))
  
  iCSV <- iCSV[,-c(1,2,5,6,8)]
  dCSV <- dCSV[,-c(1,2,5,6,8)]
  
  if (is.null(iint[[id]])){
    iint[[id]] <- iCSV[nodes,]
    itip[[id]]  <- iCSV[tips,]
    dint[[id]]  <- dCSV[nodes,]
    dtip[[id]]  <- dCSV[tips,]
  }else{
    iint[[id]] <- rbind(iint[[id]], iCSV[nodes,])
    itip[[id]]  <- rbind(itip[[id]], iCSV[tips,])
    dint[[id]]  <- rbind(dint[[id]] , dCSV[nodes,])
    dtip[[id]]  <- rbind(dtip[[id]], dCSV[tips,])
  }
}

test <- read.tree("~/work/trees/original/30647-a_100_569.tree")

indexes <- lapply(itip, function(x){
  table(x[x$count > 0, "index"])
})




trefiles<- Sys.glob(paste0(args[1],"*"))
require(ape)
for (t in 1:length(trefiles)){
  filename <- basename(trefiles[t])
  nwk <- read.tree(trefiles[t])

  nwk$edge.length <- rep(1, (Ntip(nwk) + Nnode(nwk) - 1))
  write.tree(nwk, file=paste0("~/work/trees/rescaled/",filename))
}

#new <- read.tree("~/work/trees/rescaled/30647-a_100_569.tree")

