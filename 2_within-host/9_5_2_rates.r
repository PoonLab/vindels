library("rjags")
library("data.table")
# Load 9_6 rsession
# 9_6 loads finalized indel count and timing data 

iint <- as.data.frame(rbindlist(iint))
itip <- as.data.frame(rbindlist(itip))
dint <- as.data.frame(rbindlist(dint))
dtip <- as.data.frame(rbindlist(dtip))

# type 
  # vloop
    # data frame

final.data <- list(itip, iint, dtip, dint)

final.data <- lapply(final.data, function(x){
  
  res <- t(sapply(x$pat, function(str){
    strsplit(str, "_")[[1]]
    
  }))
  x[, 'pat'] <- res[,1]
  x[, 'rep'] <- res[,2]
  split(x, x$vloop)
})


# --- Preprocessing of Data ---

final.data <- lapply(1:4, function(x){
  lvl1 <- final.data[[x]]
  
  lapply(1:5, function(vloop){
    vloop <- final.data[[x]][[vloop]]
    split(vloop , vloop$rep)
  })
})


# ---- Simulate Data JAGS ---- 

rate <- 0.00001
time <- round(rlnorm(50000,4.5,0.3))

lambda <- rate * time

y <- rpois(50000,lambda=lambda)

data <- data.frame(time, y)

# ---- STANDARD GLM ----
fit <- glm(y ~ 1, offset=log(time), data=df, family='poisson')
exp(coef(fit)[[1]])


mod_string <- " 
model {
  for (i in 1:N){
    counts[i] ~ dpois(lam[i])
    log(lam[i]) = int + log(length[i])
  }

  int ~ dunif(0,0.1)
  
}"

set.seed(113)
params = c("int")
forJags <- list(counts = data$y,
                length = data$time,
                N = nrow(data))

model.fit <- jags.model(textConnection(mod_string), 
                        data = forJags,
                        n.chains=3)
update(model.fit, 100)

mod.sim = coda.samples(model=model.fit,
                       variable.names=params,
                       n.iter=1000)

mod.csim <- as.mcmc(do.call(rbind, mod.sim))

# ---- Simulate Data RSTAN ---- 

# Data import
data.stan <- list(N=50000,
                  counts=data$y,
                  branches=data$time)

# Stan modeling 
start <- proc.time()
stan.fit <- stan("~/vindels/2_within-host/stan_modeling/rate-new.stan",
                 data= data.stan, 
                 chains=1,
                 iter=2000,
                 control=list(adapt_delta=0.90))
end <- proc.time() - start


# ---- REAL DATA RSTAN ---- 

df <- final.data[[1]][[1]][[1]]

# Data import
data.stan <- list(N=nrow(df),
                  counts=df$count,
                  branches=df$length)

# Stan modeling 
start <- proc.time()
stan.fit <- stan("~/vindels/2_within-host/stan_modeling/rate-new.stan",
                 data= data.stan, 
                 chains=1,
                 iter=1000,
                 control=list(adapt_delta=0.90))
end <- proc.time() - start


# ---- REAL DATA RSTAN ITERATED ----

data <- final.data[[1]][[1]]


options(mc.cores = parallel::detectCores()-2)
# Stan modeling 
start <- proc.time()

vmean <- list()
vmed <- list()
vneff <- list()

for (v in 1:5){
  data <- final.data[[1]][[v]]
  mean <- c()
  median <- c()
  neff <- c()
  
  for (i in 1:200){
    # Data import
    df <- data[[i]]
    stan.df <- list(N=nrow(df),
               counts=df$count,
               branches=df$length)
    
    stan.fit <- stan("~/vindels/2_within-host/stan_modeling/rate-new.stan",
                     data=stan.df, 
                     chains=1,
                     iter=10000,
                     control=list(adapt_delta=0.90),
                     verbose=F)
    sum <- summary(stan.fit)$summary
    mean[i] <- exp(sum[1,1]) * (365 / median(df$vlen))
    median[i] <- exp(sum[1,6]) * (365 / median(df$vlen))
    neff[i] <- sum[1,'n_eff']
    print(paste0("Finished: ", i))
  }
  vmean[[v]] <- mean
  vmed[[v]] <- median
  vneff[[v]] <- neff
}


end <- proc.time() - start

lens <- lapply(1:5, function(v){
  sapply(1:200, function(x){
    median(final.data[[1]][[v]][[x]]$vlen)
  })
})

final <- lapply(1:5, function(v){
  sapply(1:200, function(x){
    exp(vmean[[v]][x]) * (365 / median(lens[[v]][x]))
  })
})



w <- extract(stan.fit)
plot(density(w$rate))
line(v=rate)
