library(rstan)

rstan_options(auto_write = TRUE)
n.cores <- parallel::detectCores()
options(mc.cores = n.cores - 1 )

set.seed(1234)


# ---- Simulate "true" data ----

np = 2  # Number of patients
n = 10   # Number of data per patient


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

z <- stan('hierarchical-stan.stan', 
     data   = data.stan, 
     chains = 1, 
     iter   = 1e4)

rstan::summary(z)
traceplot(z, alpha=0.7)
w <- extract(z)

plot(density(w$mlog))
abline(v=mm)


plot(density(w$lambda))
