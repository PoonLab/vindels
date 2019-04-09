require(bbmle)
require(stringr)
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


#TESTING ----------------------
csvfile <- read.csv(file, sep="\t", stringsAsFactors = F)
output <- sapply(csvfile$Ins, extractInfo)



# INSERTION PARSING ----------
folder <- Sys.glob("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/insertions/*.csv")
data.df <- data.frame()
total.df <- list()
count <- 0
for (file in folder){
  print(file)
  count <- count + 1
  csvfile <- read.csv(file,sep="\t", stringsAsFactors = F)
  if (all(is.na(csvfile$Ins))){
    csvfile$Ins <- ""
  }
  csvfile$Count <- sapply(csvfile$Ins, csvcount) 
  output <- sapply(csvfile$Ins, extractInfo)
  output <- t(output)
  rownames(output) <- NULL
  output <- as.data.frame(output)
  output$V1 <- as.character(output$V1)
  output$V2 <- as.character(output$V2)
  csvfile <- cbind(csvfile, output)
  csvfile$Ins <- NULL
  colnames(csvfile) <- c("Accno", "Date", "Vloop", "Count", "Seq", "Pos")
  vrdf <- split(csvfile, csvfile$Vloop)
  
  for (i in 1:5){
    total.df[i][[1]] <- rbind(total.df[i][[1]], vrdf[i][[1]])
  }
}

# RATE ANALYSIS -------------
rates <- c()
for (vloop in 1:5){
  obj.f <- function(rate) -poisll(rate, total.df[[vloop]]$Count, total.df[[vloop]]$Date)
  mle.result <- bbmle::mle2(obj.f, start=list(rate=0.1), method = "Brent", lower=1e-12, upper = 0.1)
  rates <- c(rates, coef(mle.result)[[1]])
}

# Get raw insertion counts 
sum(total.df[[1]]$Count)
sum(total.df[[2]]$Count)


