source("~/vindels/2_within-host/utils.r")

# INSERTION PARSING ----------
path <- "~/Lio/"
path <- "~/PycharmProjects/hiv-withinhost/"
ifolder <- Sys.glob(paste0(path,"9Indels/rep/wholetree/ins/*.tsv"))
dfolder <- Sys.glob(paste0(path,"9Indels/rep/wholetree/del/*.tsv"))



for (file in 1:length(ifolder)){
  print(file)
  filename <- strsplit(basename(ifolder[file]),"\\.")[[1]][1]
  runno <- strsplit(filename, "_")[[1]][2]
  count <- count + 1
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F, sep='\t')
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F, sep='\t')
  
  for (i in 2:6){
    res <- unname(sapply(iCSV[,i], function(x){csvcount(x,":")}))
    iCSV[,i] <- res
    res <- unname(sapply(dCSV[,i], function(x){csvcount(x,":")}))
    dCSV[,i] <- res
  }
  
  # reads in the tree
  tre <- read.tree(paste0(paste0(path,"7SampleTrees/prelim/108869-a_10_368.tree.sample")))#,filename , ".tree.sample")))
  lens <- node.depth.edgelength(tre)    # [(length(tre$tip.label)+1):(length(tre$edge.length)+1)]  #used if you want to only access internal nodes and not tips

  # remove the root from the rtt length vector because it is NOT found in the reconstruction or the indel extraction
  #lens <- lens[-(Ntip(tre)+1)]
  res <- unname(sapply(iCSV$header, function(x){
    tips <- str_match_all(x,"([^\\)\\(,\n:]+):")[[1]][,2]
    if (length(tips) == 0){
      # the index in the tre$tip.label vector is the final result
      index <- match(x, tre$tip.label)
    }else{
      desc <- Descendants(tre)

      # find the numeric labels of all extracted tips 
      matches <- match(tips, tre$tip.label)
      
      # find the SINGLE node in the descendants list that contains the exact same subset of tips
      index <- which(sapply(desc, function(x){ifelse(length(x) == length(matches) && all(x==matches),T,F)}))
    }
    if (length(index)!=1){
      return(NaN)
    }
    return(lens[index])
  }))
  
  iCSV$length <- res
  dCSV$length <- res
  
  
  # names are in the order of the sequences as they appear in the Historian file 
  
  tip.rtt <- node.depth.edgelength(tre)[1:Ntip(tre)]
  names(tip.rtt) <- tre$tip.label
  
  tip.rows <- which(tre$edge[,2] <= Ntip(tre))
  tip.nodes <- tre$edge[tip.rows,1]
  anc.rtt <- node.depth.edgelength(tre)[tip.nodes]
  names(anc.rtt) <- tre$tip.label
  
  iheaders <- unname(sapply(iCSV$header, function(x){gsub("_\\d+$","",x)[[1]]}))
  dheaders <- unname(sapply(dCSV$header, function(x){gsub("_\\d+$","",x)[[1]]}))
  
  iCSV$mid.rtt <- (tip.rtt[iheaders] + anc.rtt[iheaders]) / 2
  dCSV$mid.rtt <- (tip.rtt[dheaders] + anc.rtt[dheaders]) / 2
  
  
  # adjusts the tre tip labels to match the accession numbers
  #tre$tip.label <- unname(sapply(tre$tip.label, function(x){strsplit(x,"_")[[1]][1]}))
  tre <- read.tree(paste0(paste0(path,"7SampleTrees/prelim/",filename , ".tree.sample")))
  # retrieves branch lengths from the tree

  
  # matches the branch length to each of the sequences 
  iCSV$Date <- tre$edge.length[match(sub("_\\d*$","",iCSV$header), tre$tip.label)]
  dCSV$Date <- tre$edge.length[match(sub("_\\d*$","",dCSV$header), tre$tip.label)]
}
iTotal <- split(csv.ins, csv.ins$Run)
dTotal <- split(csv.del, csv.del$Run)