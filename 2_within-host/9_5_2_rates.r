library("rjags")

# Load 9_6 rsession

iint <- as.data.frame(rbindlist(iint))
itip <- as.data.frame(rbindlist(itip))
dint <- as.data.frame(rbindlist(dint))
dtip <- as.data.frame(rbindlist(dtip))

# type 
  # vloop
    # data frame

final.data <- list(itip, iint, dtip, dint)

final.data <- lapply(final.data, function(x){
  x[,"pat"] <- sapply(x$pat, function(str){
    strsplit(str, "-")[[1]][1]
  })
  split(x, x$vloop)
})



# ---- Simulate Data JAGS ---- 

rate <- 0.0001
time <- rexp(10000,rate=0.05)

lambda <- exp(rate) * time

y <- rpois(10000,lambda=0.1)

data <- data.frame(time, y)

mod_string <- " 
model {
  for (i in 1:N){
    counts[i] ~ dpois(lam[i])
    log(lam[i]) = log(length[i]) + int
  }

  int ~ dlnorm(log(0.0001),0.5)
  
}"

set.seed(113)
params = c("int")
forJags <- list(counts = data$y,
                length = data$time,
                N = nrow(data))

model.fit <- jags.model(textConnection(mod_string), 
                        data = forJags,
                        n.chains=1)
update(model.fit, 1000)

mod.sim = coda.samples(model=model.fit,
                       variable.names=params,
                       n.iter=10000)

mod.csim <- as.mcmc(do.call(rbind, mod.sim))

# ---- Simulate Data RSTAN ---- 

# Data import
data.stan <- list(N=1000,
                  counts=data$y,
                  branches=data$time)

# Stan modeling 
start <- proc.time()
stan.fit <- stan("~/vindels/2_within-host/stan_modeling/rate-new.stan",
                 data= data.stan, 
                 chains=1,
                 iter=50000,
                 control=list(adapt_delta=0.99))
end <- proc.time() - start

