library("rjags")
library("data.table")

# iint <- as.data.frame(rbindlist(iint))
# itip <- as.data.frame(rbindlist(itip))
# dint <- as.data.frame(rbindlist(dint))
# dtip <- as.data.frame(rbindlist(dtip))

# LOAD 9_6_unfixed here 

# type 
  # vloop
    # data frame

final.data <- list(itip, iint, dtip, dint)

final.data <- lapply(final.data, function(x){
  
  # Filter for length 0 
  idx <- x$length == 0
  x[idx,'length'] <- 0.001
  
  # Extract patient strings
  res <- t(sapply(x$pat, function(str){
    
    parsed <- strsplit(str, "_")[[1]]
    # if(grepl("-b_", str)){
    #   parsed[2] = paste(as.numeric(parsed[2]) + 200)
    # }
    # parsed[1] <- strsplit(parsed[1], "-")[[1]][1]
    parsed[1] <- substr(parsed[1], 1, nchar(parsed[1])-2)
    parsed
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


# ---- Simulate Data---- 

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
                     control=list(adapt_delta=0.90))
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
    1000* exp(vmean[[v]][x]) * (365 / median(lens[[v]][x]))
  })
})

# ----- REAL DATA -----
# ------------------------

# ---- STANDARD GLM ----
df <- final.data[[1]][[1]]
fit <- glm(count ~ 1, offset=log(length), data=df, family='poisson')
coef(fit)[[1]]

fit2 <- glm(count ~ offset(log(length)), data=df, family='poisson')
coef(fit2)[[1]]



# ---- REAL DATA RSTAN ITERATED ----

data <- final.data[[1]][[1]]


options(mc.cores = parallel::detectCores()-2)
# Stan modeling 
start <- proc.time()

mean <- c()
median <- c()
neff <- c()


for (v in 1:5){
  # Load data 
  data <- final.data[[1]][[v]]
  data$rep <- as.numeric(data$rep)
  
  # Dimensions of matrices 
  npat = length(unique(data$pat))
  ntree = 50 # length(unique(data$rep))
  
  # Count number of branches per patient 
  sizes <- c(table(data[data$rep=="1",'pat']))
  
  # Split by patient 
  byPat <- split(data, data$pat)
  byPat <- lapply(byPat, function(x){
    x[order(x$rep),]
  })
  
  # Search for empty patients 
  pos = 1 
  sums <- c()
  for(i in 1:npat){
    start = pos
    end = pos + sizes[i]-1
    
    sums[i] <- sum(byPat[[i]]$count)
    
    pos = pos + sizes[i]
  }
  
  # Modify sizes to remove the empty patients
  sizes <- sizes[-which(sums==0)]
  npat <- npat - sum(sums==0)
  byPat <- byPat[-which(sums==0)]

  # Initialize matrices
  all.counts <- matrix(nrow=sum(sizes),
                       ncol=ntree)
  all.times <- matrix(nrow=sum(sizes),
                      ncol=ntree)

  # Load matrices
  pos = 1
  for(i in 1:npat){
    start = pos
    end = pos + sizes[i]-1
    print(paste0(start," ",end))
    all.counts[start:end,] = byPat[[i]]$count[sample(ntree)]
    all.times[start:end,] = log(byPat[[i]]$length[sample(ntree)])

    pos = pos + sizes[i]
  }

  stan.df <- list(npat=npat,
                  ntree=ntree,
                  sizes=sizes,
                  counts=all.counts,
                  l_branches=all.times)

  stan.fit <- stan("~/vindels/2_within-host/stan_modeling/rate-perpat-dc.stan",
                   data=stan.df,
                   chains=1,
                   iter=10000,
                   #control=list(adapt_delta=0.90)
                   )
  end <- proc.time() - start

  # Save summary
  sum <- summary(stan.fit)$summary

  # Load results
  mean[i] <- exp(sum[1,1]) * (365 / median(df$vlen))
  median[i] <- exp(sum[1,6]) * (365 / median(df$vlen))
  neff[i] <- sum[1,'n_eff']

  print(paste0("Finished: ", i))

}




lens <- lapply(1:5, function(v){
  sapply(1:200, function(x){
    median(final.data[[1]][[v]][[x]]$vlen)
  })
})

final <- lapply(1:5, function(v){
  sapply(1:200, function(x){
    1000* exp(vmean[[v]][x]) * (365 / median(lens[[v]][x]))
  })
})





w <- extract(stan.fit)
plot(density(w$rate))
line(v=rate)
