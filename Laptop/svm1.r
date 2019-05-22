# code for machine learning in R 
# packages e1017 and penalizedSVM

require(e1071)
require(penalizedSVM)
require(ape)
require(Biostrings)


folder <- Sys.glob("~/Github/vindels/Laptop/VariableLoopSeqs/*.csv")

nglycs <- data.frame(stringsAsFactors = F)

for (file in folder){
  file <- folder[1]
  
  data <- read.csv(file, stringsAsFactors = F)
  
  # remove the unneeded start and stop columns
  data <- data[,-c(3,4,6,7,9,10,12,13,15,16)]
  
  for (seq in 1:nrow(data)){
    accno <- data[seq,1]
    ngCounts <- c()
    for (idx in 1:5){
      vloop <- data[seq,idx+1]

      aaseq <- as.character(trans(as.DNAbin(DNAString(vloop))))[[1]]
      aaseq <- paste(aaseq,collapse="")

      search <- gregexpr("N[^P][ST][^P]", aaseq)[[1]]
      ngCounts <- c(ngCounts, length(search))
    }
    nglycs <- rbind(nglycs, data.frame(accno=accno, v1=ngCounts[1],v2=ngCounts[2],v3=ngCounts[3],v4=ngCounts[4],v5=ngCounts[5]))
    
  }
    
}


