require(ape)
require(stringr)
require(phangorn)
require(data.table)
require(bbmle)
source("~/vindels/2_within-host/utils.r")

# INSERTION PARSING ----------
path <- "~/PycharmProjects/hiv-withinhost/"
ifolder <- Sys.glob(paste0(path,"9Indels/main/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"9Indels/main/del/*.tsv"))
sep <- "\t"
trefolder <- paste0(path,"7SampleTrees/prelim200/")

# CASE: Removed patient 56552 because could not complete Historian runs 
# CASE: Removed patients 49641 and 56549 because they are SUBTYPE B
# CASE: Removed patient 28376 and B because of very bad Rsquared value 
reg <- "56552|49641|56549|28376"
newreg <- "56549|49641"

ifolder <- ifolder[grepl(newreg,ifolder)]
dfolder <- dfolder[grepl(newreg,dfolder)]

iint <- list()
itip <- list()
dint <- list()
dtip <- list()


for (file in 1:length(ifolder)){
  print(file)
  filename <- basename(ifolder[file])
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
  dCSV$index <- res
  idx <- seq(1,length(res), 5)
  if (file == 1){
    dict <- lapply(1:length(unique(res)), function(x){
      iCSV[idx[x], "header"]
    })  
  }
  
  iCSV$length <- tre$edge.length[match(res, tre$edge[,2])]
  dCSV$length <- iCSV$length
  
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
  }
}


# Paul's data extraction 

tipnode <- Ntip(tre) + Nnode(tre)
n <- Ntip(tre)

# insertions 
vec <- iCSV[iCSV$count >0 ,"index"]
tb <- tabulate(unique(vec),nbins=tipnode)
boolean <- tb > 0

vec <- dCSV[dCSV$count >0 ,"index"]
tb <- tabulate(unique(vec),nbins=tipnode)
boolean <- tb > 0

tre$tip.label <- 1:n
tre$node.label <- (n+1):tipnode

tre$tip.label[!boolean[1:n]] <- ""
tre$node.label[!boolean[(n+1):tipnode]] <- ""
  
  
  which(grepl("OC", ifolder))
# clear the values of
for (i in which(grepl("OC", patnames))){
  itip[[i]] <- NA
  dtip[[i]] <- NA
}
  