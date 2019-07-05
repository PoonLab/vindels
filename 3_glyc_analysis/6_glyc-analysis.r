require(ape)
require(stringr)
translate <- function(dna) {
  require(ape)
  
  if (nchar(dna) %% 3 != 0) {
    return(NA)
  }
  dnabin <- as.DNAbin.DNAString(dna)
  aabin <- trans(dnabin)[[1]]
  aaseq <- paste(as.character(aabin),collapse="")
  aaseq
}


extractGlycs <- function(dna){
  require(ape)
  
  if (nchar(dna) %% 3 != 0) {
    return(NA)
  }
  dnabin <- as.DNAbin.DNAString(dna)
  aabin <- trans(dnabin)[[1]]
  aaseq <- paste(as.character(aabin),collapse="")
  result <- gregexpr("N[^P][ST][^P]", aaseq)[[1]] * 3 - 2
  paste(result,collapse=",")
}

extractStage <- function(x){
  temp <- strsplit(x, "\\.")[[1]][5]
  if (temp == "AIDS"){
    temp <- "chronic"
  }else if(temp == "acute_infection"){
    temp <- "acute"
  }
  temp 
}
countGlycs <- function(str){
  if (grepl(",", str)){
    return(str_count(str, ",")+1)
  }else if (str == ""){
    return(0)
  }else{
    return(1)
  }
}

tre <- read.tree("~/PycharmProjects/glyc-analysis/6_pruned/pruned.tree")

full.len <- read.csv("~/PycharmProjects/glyc-analysis/3_sequences/full/gp120.csv", stringsAsFactors = F)
# retrieve the sequences from the tree 


idx <- match(tre$tip.label, full.len$header)
# look up and retrieve their FULL LENGTH sequences 

fasta <- full.len[idx,]
subtype <- unname(sapply(fasta$header, function(x){strsplit(x, "\\.")[[1]][1]}))
accno <- unname(sapply(fasta$header, function(x){strsplit(x, "\\.")[[1]][2]}))

nglycs <- data.frame(accno=accno, subtype=subtype, stage=unname(sapply(fasta$header, extractStage)),  pos=unname(sapply(fasta$seq, extractGlycs)), stringsAsFactors = F)
nglycs$count <- unname(sapply(nglycs$pos, countGlycs))


# analyze their sequences for glyc sites 