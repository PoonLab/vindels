require(ape)
require(stringr)
removeNA <- function(input, repl=""){
  if (is.na(input)){
    input <- repl
  }
  input
}

extractGlycs <- function(seq){
  dnabin <- as.DNAbin.DNAString(seq)
  aabin <- trans(dnabin)[[1]]
  aaseq <- paste(as.character(aabin),collapse="")
  
  result <- gregexpr("N[^P][ST][^P]", aaseq)
  result
}



ins <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/insertions.csv", row.names=1, stringsAsFactors = F)
del <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/deletions.csv", row.names=1, stringsAsFactors = F)

ins$Pos <- sapply(ins$Pos, removeNA)
del$Pos <- sapply(del$Pos, removeNA)


ins$glycs <- sapply(ins$Vseq, extractGlycs)
del$glycs <- sapply(del$Vseq, extractGlycs)