
# LOAD 9_6_unfixed here 

# type 
# vloop
# data frame

final.data <- list(itip, iint, dtip, dint)

final.data <- lapply(final.data, function(x){
  
  # Filter for length 0 
  x <- x[x$length>0,]
  
  # Extract patient strings
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

# CHECK FOR LENGTH 0 ENTRIES
test <- lapply(1:4, function(x){
  lapply(1:5, function(v){
    lapply(1:200, function(rep){
      final.data[[x]][[v]][[rep]]$length
    })
  })
})
sum(unlist(x) == 0)


# ---- Simulate Data JAGS ---- 

rate <- 0.00001
time <- round(rlnorm(50000,4.5,0.3))

lambda <- rate * time

y <- rpois(50000,lambda=lambda)

data <- data.frame(time, y)

# ---- STANDARD GLM ----
fit <- glm(y ~ 1, offset=log(time), data=df, family='poisson')
coef(fit)[[1]]

fit2 <- glm(y ~ offset(log(time)), data=df, family='poisson')
coef(fit2)[[1]]




mod_string <- " 
model {
  for (i in 1:N){
    counts[i] ~ dpois(lam[i])
    log(lam[i]) = int + log(length[i])
  }

  int ~ dnorm(0,5)
  
}"


# ---- RJAGS REAL DATA (ITERATED) ----
set.seed(113)


start <- proc.time()
params = c("int")
vmean <- list()
vmed <- list()
vneff <- list()

for (v in 1:5){
  data <- final.data[[1]][[v]]
  mean <- c()
  median <- c()
  neff <- c()
  
  for (r in 1:200){
    df = data[[r]]
    
    forJags <- list(counts = df$count,
                    length = df$length,
                    N = nrow(df))
    
    model.fit <- jags.model(textConnection(mod_string), 
                            data = forJags,
                            n.chains=1)
    update(model.fit, 500)
    
    mod.sim = coda.samples(model=model.fit,
                           variable.names=params,
                           n.iter=5000)
    
    mod.csim <- as.mcmc(do.call(rbind, mod.sim))
  }
  vmean[[v]] <- mean
  vmed[[v]] <- median
  vneff[[v]] <- neff
}



