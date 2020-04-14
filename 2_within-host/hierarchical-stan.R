library(rstan)

rstan_options(auto_write = TRUE)
n.cores <- parallel::detectCores()
options(mc.cores = n.cores - 2 )

set.seed(1234)


# ---- Simulate "true" data ----

np = 20  # Number of patients
n = 100   # Number of data per patient


mm <- 0.25
lambda.true = rlnorm(n=np, 
                    meanlog = log(mm), 
                    sdlog =  0.3)  
mean(lambda.true)
hist(lambda.true, col='grey')

M <- matrix(nrow = n, ncol = np)

for(i in 1:np){
  M[,i] <- rexp(n=n, rate = lambda.true[i])
}

# ---- Stan Fit ----

data.stan <- list(np = np,
                  n = n,
                  M = M)

z <- stan('~/vindels/2_within-host/hierarchical-stan.stan', 
     data   = data.stan, 
     chains = 1, 
     iter   = 1e4)

rstan::summary(z)
#traceplot(z, alpha=0.7)

w <- extract(z)

par(mfrow=c(1,2))
plot(density(w[[1]]))
abline(v=log(mm), col="red")
plot(density(w[[2]]))
abline(v=median(mat), col="red")


plot(density(w$lambda))
