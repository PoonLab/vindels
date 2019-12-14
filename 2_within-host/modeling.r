

require(bbmle)
#slips <- matrix(rep(rep(0,121),10), nrow=10,ncol=121)
slips <- c(rep(0,121000))
slips[sample(121000,100)] <- sample(3,100,replace = T) * 3

affinell <- function(p.slip, p.stay, slips){
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  #(1 - p.slip)^x * p.slip^y * (1-p.stay)^y * p.stay^z
  2 * y * x * z * log(1-p.slip) * log(p.slip) * log(p.stay)
}

objf <- function(p.slip, p.stay){
  -affinell(p.slip, p.stay, slips)
}

result <- mle2(objf, start=list(p.slip=1,p.stay=1), method = "L-BFGS-B", lower=1e-12, upper = 1)


geomll <- function(p.copy){
  N <- length(counts)
  # log likelihood of the geometric distribution
  N * log(p.copy) + sum(counts) * log(1-p.copy)
}

counts <- c(0,0,0,0,0,0,0,0,0,0,1,0)

likelihood <- function(param){
  sum(dgeom(counts,prob=param, log=T))
}

prior <- function(param){
  prior <- dunif(param, log = T)
  return(prior)
}


posterior <- function(param){
  prior(x) + likelihood(param)
}

x <- rep(0,1000)
x[1] <- 3

for (i in 2:1000){
  currentx <- x[i-1]
  proposedx <- currentx + rnorm(1, mean=0, sd=1)
  prop <- target(proposedx) / target(currentx)
  # if the proportion exceeds the random uniform sample, accept the proposed value
  if (runif(1) < prop) {
    x[i] <- proposedx
  
  # if the proportion is less than the random uniform sample, reject the proposed value stick with current 
  } else {
    x[i] <- currentx
  }
  
}