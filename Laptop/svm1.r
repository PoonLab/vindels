# code for machine learning in R 
# packages e1017 and penalizedSVM

require(e1071)
require(penalizedSVM)
require(ape)
require(Biostrings)


csvfile <- read.csv("~/PycharmProjects/glyc-analysis/3_glycs/glycs.csv")

nglycs <- data.frame(stringsAsFactors = F)



csvfile[,1] <- as.character(csvfile[,1])

subtype <- sapply(csvfile[,1], function(x){strsplit(x, "\\.")[[1]][1]})
subtype <-unname(subtype)

accno <- sapply(csvfile[,1], function(x){strsplit(x, "\\.")[[1]][2]})
accno <- unname(accno)

status <- sapply(csvfile[,1], function(x){
  temp <- strsplit(x, "\\.")[[1]][5]
  if (temp == "AIDS"){
    temp <- "chronic"
  }else if(temp == "acute_infection"){
    temp <- "acute"
  }
  temp 
  })
status <- unname(status)

csvfile[,1] <- NULL
csvfile <- cbind(accno, subtype, status, csvfile)

# spl <- sample(nrow(csvfile), 610)
# glyc_train <- csvfile[spl,-c(1,2)]
# glyc_test <- csvfile[-spl,-c(1,2)]
# 
# svmfit <- svm(status ~ ., data=glyc_train, kernel="polynomial", cost=0.1, scale=F)
# print(svmfit)
# 
# tuned <- tune(svm, status~., data=glyc_train, kernel="linear", ranges=list(cost=c(0.001, 0.01, 0.1,1,10,100)))
# 
# p <- predict(svmfit, glyc_test, type="class")


subC <- csvfile[csvfile$subtype=="C",]
subB <- csvfile[csvfile$subtype=="B",]

require(reshape)

newC <- melt(subC, id=c("accno","subtype","status"))
newB <- melt(subB, id=c("accno","subtype","status"))

boxplot(value ~ variable + status, data=newB)
boxplot(value ~ variable + status, data=newC, main="Subtype C")
abline(v=5.5)



