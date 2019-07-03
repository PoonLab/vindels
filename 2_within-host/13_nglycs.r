require(ape)
require(stringr)
removeNA <- function(input, repl=""){
  if (is.na(input)){
    input <- repl
  }
  input
}

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

extractGlycs <- function(seq){
  dnabin <- as.DNAbin.DNAString(seq)
  aabin <- trans(dnabin)[[1]]
  aaseq <- paste(as.character(aabin),collapse="")
  
  result <- gregexpr("N[^P][ST][^P]", aaseq)
  result
}



ins <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/insertions.csv", row.names=1, stringsAsFactors = F)
del <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/deletions.csv", row.names=1, stringsAsFactors = F)

#ins$Pos <- ins$Pos + 1
#del$Pos <- del$Pos + 1


ins$Seq <- sapply(ins$Seq, translate)
del$Seq <- sapply(del$Seq, translate)

ins$Vseq <- sapply(ins$Vseq, translate)
del$Vseq <- sapply(del$Vseq, translate)



ins$glycs <- sapply(ins$Vseq, extractGlycs)
del$glycs <- sapply(del$Vseq, extractGlycs)


ins$