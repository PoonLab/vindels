require(ape)

poissonlogll <- function (rate, outcomes, times){
  prob <- exp(-rate*times)  # e ^ -rate*t
  sum(outcomes*log(prob), (1-outcomes)*log(1-prob))
}

csvcount <- function(input){
  if (input != ""){
    count <- length(strsplit(as.character(input), ","))
  }else{
    count <- 0
  }
}

extractInfo <- function(input){
  headers <- strsplit(as.character(input), "\\.")
  name <- c()
  date <- c()
  for (i in c(1,2)){
    name <- c(name, strsplit(headers[[1]][i],"_")[[1]][1])
    date <- c(date, strsplit(headers[[1]][i],"_")[[1]][2])
  }
  
  list(name1=name[1], name2=name[2], date1=date[1], date2=date[2])
}

#count1 <- sapply(data1, function(x){ length(strsplit(as.character(x),","))})
#count1 + sapply(data2, function(x){ length(strsplit(as.character(x),","))})

reconfolder <- Sys.glob("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/*")

data.df <- data.frame(name=character(), date=numeric(), count=integer())
#for (file in reconfolder){
file <- "~/PycharmProjects/hiv-withinhost/9Indels/30631-a_recon.csv"
recon <- read.csv(file, sep="\t")
  
insertion <- data.frame(recon[,c(1:3)])
insertion$count1 <- sapply(insertion$Ins1, csvcount) 
insertion$count2 <- sapply(insertion$Ins2, csvcount)
  
insertion$name1 <- unlist(insertion$name1)
insertion$name2 <- unlist(insertion$name2)
insertion$date1 <- unlist(insertion$date1)
insertion$date2 <- unlist(insertion$date2)
  
insertion <- cbind( as.data.frame(t(sapply(insertion$AccNo, extractInfo))), insertion)
insertion$AccNo <- NULL
  
data.df <- rbind(data.df, data.frame(name=insertion$name1, date=insertion$date1, ins=insertion$Ins1, count=insertion$count1))
data.df <- rbind(data.df, data.frame(name=insertion$name2, date=insertion$date2, ins=insertion$Ins2, count=insertion$count2))
#} 


