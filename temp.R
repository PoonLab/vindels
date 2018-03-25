> mle.result <- mle (obj.f, start=list(rate=1), method = "Brent", lower=1e-12, upper = 1)
> mle.result

Call:
  mle(minuslogl = obj.f, start = list(rate = 1), method = "Brent", 
      lower = 1e-12, upper = 1)

Coefficients:
  rate 
0.01421269 