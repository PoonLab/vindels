require(bbmle)
require(stringr)
require(ape)
poisll <- function (rate, counts, times){
  prob <- (exp(-rate*times))*((rate*times)^(counts))/(factorial(counts))    # (e^-rt)*((rt)^k)/(k!)
  sum(log(prob))
}


csvcount <- function(input){
  commas <- str_count(input, ",")
  if (commas > 0){
    result <- commas + 1  
  }else if(input == ""){
    result <- 0
  }else{
    result <- 1
  }
}

extractInfo <- function(input){
  if (length(input)==1 && input == ""){
    return(c("",""))
  }else{
    insertions <- str_split(input, ",")
  }
  seq <- c()
  pos <- c()
  
  for (ins in insertions[[1]]){
     fields <- str_split(ins, "-")
     seq <- c(seq, fields[[1]][1])
     pos <- c(pos, as.numeric(fields[[1]][2]))
  }
  return(c(paste(seq,collapse=","), paste(pos,collapse=",")))
}

cutHeader <- function(header){
  newheader <- str_split(as.character(header),"_")[[1]][1]
  newheader
}


#TESTING ----------------------
csvfile <- read.csv(file, sep="\t", stringsAsFactors = F)
output <- sapply(csvfile$Ins, extractInfo)



# INSERTION PARSING ----------
folder <- Sys.glob("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/insertions/*.csv")
treedir <- "/home/jpalmer/PycharmProjects/hiv-withinhost/7SampleTrees/prelim/"
data.df <- data.frame()
total.df <- list()
count <- 0
for (file in folder){
  print(file)
  basename <- str_split(basename(file),"\\.")[[1]][1]
  count <- count + 1
  csvfile <- read.csv(file,sep="\t", stringsAsFactors = F)
  if (all(is.na(csvfile$Ins))){
    csvfile$Ins <- ""
  }
  csvfile$Count <- sapply(csvfile$Ins, csvcount) 
  
  tre <- read.tree(paste0(treedir, basename, ".tree.sample"))
  branches <- tre$edge.length[tre$edge[,2] <=Ntip(tre)]   #branches will match exactly with the tre$tip.label order
  csvfile$Date <- branches[match(csvfile$Accno, tre$tip.label)]
  csvfile$Accno <- unname(sapply(csvfile$Accno, cutHeader))
  
  output <- sapply(csvfile$Ins, extractInfo)
  output <- unname(output)
  output <- t(output)
  output <- as.data.frame(output)
  output$V1 <- as.character(output$V1)
  output$V2 <- as.character(output$V2)
  csvfile <- cbind(csvfile, output)
  csvfile$Ins <- NULL
  
  colnames(csvfile) <- c("Accno", "Vloop", "Count","Date", "Seq", "Pos")
  vrdf <- split(csvfile, csvfile$Vloop)
  
  for (i in 1:5){
    total.df[i][[1]] <- rbind(total.df[i][[1]], vrdf[i][[1]])
  }
}

# RATE ANALYSIS -------------
rates <- c()
vlengths <- c(84,156,105,90,33)
for (vloop in 1:5){
  #obj.f <- function(rate) -poisll(rate, total.df[[vloop]]$Count, total.df[[vloop]]$Date)
  fit <- glm(total.df[[vloop]]$Count ~ total.df[[vloop]]$Date, family="poisson") 
  rates <- c(rates, coef(fit)[2][[1]]*365/vlengths[vloop])
}

# Get raw insertion counts 
sum(total.df[[1]]$Count)
sum(total.df[[2]]$Count)


