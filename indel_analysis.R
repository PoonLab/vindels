csvfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/9_indels",full.names=TRUE)


setwd("~/PycharmProjects/hiv-evolution-master/10_analysis")
print(csvfolder)

output <- data.frame("subtype"=character(), stringsAsFactors = FALSE)
for (i in 1:length(csvfolder)){
  #binomial analysis
  
  csv <- read.csv(csvfolder[i])
  filename <- strsplit(strsplit(csvfolder[i], "/")[[1]][7], "\\.")[[1]][1]
  print(filename)
  output[i,1] <- filename
  
  if (!is.numeric(csv$tip1.len) | !is.numeric(csv$tip2.len)){
    print("skipped")
    next
  }else{
    csv$length <- csv$tip1.len + csv$tip2.len 
  }
  
  for (x in 1:5){
    fit <- handle.warning(csv, x)
    if (is.null(fit)){
      next
    }else if(!is.numeric(fit$coefficients[[2]])){
      next
    }
    
    print(csv[,(2*x+4)])
    inv.logit <- (exp(fit$coefficients[[2]]))/(exp(fit$coefficients[[2]])+1)
    
    #unsure of the logic behind this 
    inv.logit <- -log(inv.logit)
    
    header <- paste0("VR", as.character(x), ".estimate")
    
    output[i, header] <- inv.logit
  }
  break
}

handle.warning <- function(csv, x) {
  fit <- tryCatch(
    {fit <- glm(!csv[,as.integer(2*x+4)] ~ csv$length, family= 'binomial')},
    warning = function(c){
      message("A warning was thrown, still running")
      return (NULL)
    })
  return(fit)
}
