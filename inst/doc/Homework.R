## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(comment = NA)

## -----------------------------------------------------------------------------
dataset1 = CO2
head(dataset1)

## ---- echo=FALSE--------------------------------------------------------------
plot(x = dataset1$conc[1:7], y = dataset1$uptake[1:7], type = "l",
     xlim = c(95, 1000), ylim = c(0, 50), xaxs = "i", yaxs = "i",
     main = "Carbon Dioxide Uptake Rates Versus Different Ambient Carbon Dioxide \nConcentrations in Each Enviroment",
     xlab = "ambient carbon dioxide concentrations (mL/L)", ylab = "carbon dioxide uptake rates (umol/m^2/sec)")
for(i in 2:12)
{
  lines(x = dataset1$conc[(7*i-6):(7*i)], y = dataset1$uptake[(7*i-6):(7*i)])
}

## -----------------------------------------------------------------------------
plot(lm(uptake ~ log(conc) + Type + Treatment, data = dataset1))

## -----------------------------------------------------------------------------
xtable::xtable(dataset1[c(1:7, 22:28), ])

## -----------------------------------------------------------------------------
set.seed(114)
U_3 = runif(1000)
sample_3 = 2 * (1-U_3) ^ (-1/2)

## -----------------------------------------------------------------------------
hist(sample_3, prob = TRUE,
     xlim = c(0, 20), breaks = 140,
     main = "Histogram of Samples From Pareto(2,2)", xlab = "x", ylab = "density")
lines(x = seq(from = 2, to = 25, by = 0.1), y = 8 / seq(from = 2, to = 25, by = 0.1) ^ 3)

## -----------------------------------------------------------------------------
beta_7 = function(n, a, b)
{
  x = c()
  count1 = 0
  count2 = 0
  while(count2 < n)
  {
    y = runif(1)
    u = runif(1)
    count1 = count1 + 1
    if(u < y ^ (a-1) * (1-y) ^ (b-1) * (max(a-1, b-1) / (a+b-2)) ^ (1 / (a+b-2)))
    {
      x = c(x, y)
      count2 = count2 + 1
    }
  }
  return(list(x, count1, count2))
}

## -----------------------------------------------------------------------------
set.seed(514)
sample_7 = beta_7(1000, 3, 2)
sample_7[[2]]

## -----------------------------------------------------------------------------
hist(sample_7[[1]], prob = TRUE,
     xlim = c(0, 1), breaks = 50,
     main = "Histogram of Samples From Beta(3,2)", xlab = "x", ylab = "density")
lines(x = seq(from = 0, to = 1, by = 0.01), y = 12 * seq(from = 0, to = 1, by = 0.01) ^ 2 * (1 - seq(from = 0, to = 1, by = 0.01)))

## -----------------------------------------------------------------------------
set.seed(1919)
lambda_12 = rgamma(1000, shape = 4, rate = 2)
sample_12 = rexp(1000, lambda_12)
summary(sample_12)

## -----------------------------------------------------------------------------
hist(sample_12+2, prob = TRUE,
     xlim = c(0, 9), breaks = 50,
     main = "Histogram of Samples From Pareto(4,2)", xlab = "x", ylab = "density")
lines(x = seq(from = 2, to = 10, by = 0.1), y = 64 / seq(from = 2, to = 10, by = 0.1) ^ 5)

## -----------------------------------------------------------------------------
quicksort = function(x, count = 0)
{
  if (length(x) <= 1) return(list(x, count))

  xi = sample(1:length(x), 1)
  flag = x[xi]
  temp = x[-xi]

  indexs = which(temp <= flag)
  count = count + length(indexs)
  x1 = quicksort(temp[indexs], count)
  count = x1[[2]]
  x1 = x1[[1]]

  indexb = which(temp > flag)
  count = count + length(indexb)
  x2 = quicksort(temp[indexb], count)
  count = x2[[2]]
  x2 = x2[[1]]

  return(list(c(x1, flag, x2), count))
}

## -----------------------------------------------------------------------------
sample1_1 = quicksort(sample(1:10^4))[[1]]
plot(x = 1: 10^4, y = sample1_1, type = "l",
     xaxs = "i", yaxs = "i", xlim = c(0, 10^4), ylim = c(0, 10^4),
     xlab = "1:10000", ylab = "fast-sorted-sample", main = "Fast Sort Randomly Permuted Data of 1 to 10000")

## ---- echo = FALSE------------------------------------------------------------
sample1_2 = quicksort(sample(1:(2 * 10^4)))[[1]]
plot(x = 1: (2 * 10^4), y = sample1_2, type = "l",
     xaxs = "i", yaxs = "i", xlim = c(0, 2 * 10^4), ylim = c(0, 2 * 10^4),
     xlab = "1:20000", ylab = "fast-sorted-sample", main = "Fast Sort Randomly Permuted Data of 1 to 20000")

## ---- echo = FALSE------------------------------------------------------------
sample1_4 = quicksort(sample(1:(4 * 10^4)))[[1]]
plot(x = 1: (4 * 10^4), y = sample1_4, type = "l",
     xaxs = "i", yaxs = "i", xlim = c(0, 4 * 10^4), ylim = c(0, 4 * 10^4),
     xlab = "1:40000", ylab = "fast-sorted-sample", main = "Fast Sort Randomly Permuted Data of 1 to 40000")

## ---- echo = FALSE------------------------------------------------------------
sample1_6 = quicksort(sample(1:(6 * 10^4)))[[1]]
plot(x = 1: (6 * 10^4), y = sample1_6, type = "l",
     xaxs = "i", yaxs = "i", xlim = c(0, 6 * 10^4), ylim = c(0, 6 * 10^4),
     xlab = "1:60000", ylab = "fast-sorted-sample", main = "Fast Sort Randomly Permuted Data of 1 to 60000")

## ---- echo = FALSE------------------------------------------------------------
sample1_8 = quicksort(sample(1:(8 * 10^4)))[[1]]
plot(x = 1: (8 * 10^4), y = sample1_8, type = "l",
     xaxs = "i", yaxs = "i", xlim = c(0, 8 * 10^4), ylim = c(0, 8 * 10^4),
     xlab = "1:80000", ylab = "fast-sorted-sample", main = "Fast Sort Randomly Permuted Data of 1 to 80000")

## -----------------------------------------------------------------------------
timing = function(n, times = 100)
{
  t = c()
  for(i in 1:times)
  {
    t = c(t, system.time(quicksort(sample(1:n)))[1])
  }
  return(mean(t))
}

## -----------------------------------------------------------------------------
set.seed(114514)
a_n = c(timing(10^4), timing(2 * 10^4),timing(4 * 10^4), timing(6 * 10^4), timing(8 * 10^4))
a_n

## -----------------------------------------------------------------------------
t_n = c(10^4 * log(10^4), 2 * 10^4 * log(2 * 10^4), 4 * 10^4 * log(4 * 10^4), 6 * 10^4 * log(6 * 10^4), 8 * 10^4 *log(8 * 10^4))
fit1_3 = lm(a_n ~ t_n)
summary(fit1_3)

## -----------------------------------------------------------------------------
plot(x = t_n, y = a_n, pch = 20,
     xlab = expression("t"["n"] * "=O(nlog(n))"), ylab = expression("a"["n"]), main = "Averaged Computation Time in the Fast Sort Algorithm \n on Different Size of Samples")
lines(x = t_n, y = predict(fit1_3))

## -----------------------------------------------------------------------------
(2 * exp(2) - 6 * exp(1) + 2) / (-exp(2) + 4 * exp(1) - 3)

## -----------------------------------------------------------------------------
set.seed(114514)
theta_MC = c()
theta_MCa = c()
for(i in 1:1000)
{
  sample5_7 = runif(2000)
  theta_MC = c(theta_MC, mean(exp(sample5_7)))
  theta_MCa = c(theta_MCa, mean(c(exp(sample5_7[1:1000]), exp(1 - sample5_7[1:1000]))))
}

## -----------------------------------------------------------------------------
plot(density(theta_MCa), col = "red",
     xlim = c(1.68, 1.76),
     xlab = "Estimate of θ", main = "Distributions of Estimates Using Two Monte Carlo Methods")
lines(density(theta_MC), col = "blue")
legend("topright", c("Simple MC", "Antithetic Variate"), lty = c(1, 1), col = c("blue", "red"))
abline(v = exp(1) - 1)

## -----------------------------------------------------------------------------
c(sd(theta_MC), sd(theta_MCa), sd(theta_MCa) / sd(theta_MC))

## -----------------------------------------------------------------------------
1 - (sd(theta_MCa) / sd(theta_MC)) ^ 2

## -----------------------------------------------------------------------------
(2 * exp(2) - 6 * exp(1) + 2) / (-exp(2) + 4 * exp(1) - 3)

## -----------------------------------------------------------------------------
x = seq(1, 5, by = 0.01)
plot(x = x, y = x^2 / sqrt(2*pi) * exp(-x^2 / 2), type = "l",
     ylim = c(0, 0.3), ylab = "value of function", main = "g(x)",
     xaxs = "i", yaxs = "i")

## -----------------------------------------------------------------------------
plot(x = x, y = 2.5 * x^2 / sqrt(2*pi) * exp(-x^2 / 2), type = "l",
     ylim = c(0, 1), ylab = "value of function", main = "(2.5)g(x) and It's 'Close'",
     xaxs = "i", yaxs = "i")
lines(x = x, y = 2 / sqrt(2*pi) * exp(-(1-x)^2 / 2), col = "red")
lines(x = x, y = 4 / gamma(2) * (x-1) * exp(-2 * (x-1)), col = "blue") #or dgamma(x-1, 2, 2)
legend("topright", c("2.5g(x)", expression("f"[1]*"(x)"), expression("f"[2]*"(x)")), col = c("black", "red", "blue"), lty = rep(1, 3))

## -----------------------------------------------------------------------------
g13 = function(x)
{
  return(x^2 / sqrt(2*pi) * exp(-x^2 / 2))
}

f113 = function(x)
{
  return(2 / sqrt(2*pi) * exp(-(1-x)^2 / 2))
}

f213 = function(x)
{
  return(4 / gamma(2) * (x-1) * exp(-2 * (x-1)))
}

## -----------------------------------------------------------------------------
set.seed(114514)
m13 = 100000
sample13_1 = abs(rnorm(m13)) + 1
grf13_1 = numeric(m13)
for(i in 1:m13)
{
  grf13_1[i] = g13(sample13_1[i]) / f113(sample13_1[i])
}

sample13_2 = rgamma(m13, 2, 2) + 1
grf13_2 = numeric(m13)
for(i in 1:m13)
{
  grf13_2[i] = g13(sample13_2[i]) / f213(sample13_2[i])
}

## -----------------------------------------------------------------------------
c(mean(grf13_1), mean(grf13_2))

## -----------------------------------------------------------------------------
integrate(g13, lower = 1, upper = "infinity")

## -----------------------------------------------------------------------------
c(var(grf13_1), var(grf13_2))

## -----------------------------------------------------------------------------
c(integrate(function(x) {return(x^4 / 4 / sqrt(2 * pi) * exp(-x^2 + (1-x) ^ 2 / 2))}, lower = 1, upper = "infinite")$value,
  integrate(function(x) {return(x^4 / (2 * pi) / 4 / (x-1) * exp(-x^2 + 2 * (x-1)))}, lower = 1.01, upper = "infinite")$value)

## -----------------------------------------------------------------------------
F_inverse15 = function(u, j)
{
  return(-log(exp(-j / 5) - u * exp(-j / 5) * (1 - exp(-1 / 5))))
}

g15 = function(x)
{
  return(exp(-x) / (1 + x ^ 2))
}

f15 = function(x, j)
{
  return(exp(-x) / exp(-j / 5) / (1 - exp(-1 / 5)))
}

## -----------------------------------------------------------------------------
set.seed(114514)
est15 = numeric(50)

for(i in 1:50)
{
  U15 = runif(10000)

  sample15_1 = numeric(2000)
  for(j in 1:2000)
  {
    X = F_inverse15(U15[j], 0)
    sample15_1[j] = g15(X) / f15(X, 0)
  }

  sample15_2 = numeric(2000)
  for(j in 1:2000)
  {
    X = F_inverse15(U15[j + 2000], 1)
    sample15_2[j] = g15(X) / f15(X, 1)
  }

  sample15_3 = numeric(2000)
  for(j in 1:2000)
  {
    X = F_inverse15(U15[j + 4000], 2)
    sample15_3[j] = g15(X) / f15(X, 2)
  }

  sample15_4 = numeric(2000)
  for(j in 1:2000)
  {
    X = F_inverse15(U15[j + 6000], 3)
    sample15_4[j] = g15(X) / f15(X, 3)
  }

  sample15_5 = numeric(2000)
  for(j in 1:2000)
  {
    X = F_inverse15(U15[j + 8000], 4)
    sample15_5[j] = g15(X) / f15(X, 4)
  }

  est15[i] = sum(c(mean(sample15_1), mean(sample15_2), mean(sample15_3), mean(sample15_4), mean(sample15_5)))
}

## -----------------------------------------------------------------------------
mean(est15)

## -----------------------------------------------------------------------------
var(est15)

## -----------------------------------------------------------------------------
0.0970314 ^ 2 # This data comes from Exercise 5.13.

## -----------------------------------------------------------------------------
sample4 = function(n, meanlog, sdlog) # data generation
{
  return(rlnorm(n, meanlog, sdlog))
}

## -----------------------------------------------------------------------------
#test sample4()
set.seed(114514)
hist(sample4(1000, 0, 0.5), breaks = 25, freq = FALSE)
lines(x = seq(0, 5, 0.01), y = dlnorm(seq(0, 5, 0.01), sdlog = 0.5))

## -----------------------------------------------------------------------------
CI4 = function(m = 1000, n = 20, meanlog = 0, sdlog = 1, alpha = 0.05) # data analysis
{
  mu_hat = numeric(m)
  sd_hat = numeric(m)
  UCL = numeric(m)
  DCL = numeric(m)
  y = numeric(m)
  for(i in 1:m)
  {
    data = sample4(n, meanlog, sdlog)
    mu_hat[i] = mean(log(data))
    sd_hat[i] = sd(log(data)) / sqrt(n)
    UCL[i] = mu_hat[i] + qt(1 - alpha / 2, n - 2) * sd_hat[i]
    DCL[i] = mu_hat[i] + qt(alpha / 2, n - 2) * sd_hat[i]
    y[i] = 1 - (meanlog < DCL[i]) - (meanlog > UCL[i])
  }
  return(list(mu = mu_hat, sd = sd_hat, UCL = UCL, DCL = DCL, CP = mean(y)))
}

## -----------------------------------------------------------------------------
#test CI4, we show a single estimate confidence interval. Following data is the estimate mean, standard deviation, DCL, UCL of the parameter.
set.seed(114514)
temp = CI4(m = 1)
c(temp$mu, temp$sd, temp$DCL, temp$UCL)
rm(temp)

## -----------------------------------------------------------------------------
resulting4 = function(m = 1000, n = 20, meanlog = 0, sdlog = 1, alpha = 0.05) # result reporting (these parameters seem a bit confusing)
{
  temp = CI4(m, n, meanlog, sdlog, alpha)

  #plot, where we show some simulative CIs as examples.
  m1 = min(m, 20)
  plot(x = 1:m1, y = temp$UCL[1:m1], type = "l",
       ylim = c(min(temp$DCL[1:m1]), max(temp$UCL[1:m1])),
       xlab = paste("μ = ", meanlog, ", σ = ", sdlog), ylab = "CI of μ",
       main = "Some Simulative Confidence Intervals of a Lognormal Distribution")
  lines(x =1:20, y = temp$DCL[1:m1])
  abline(h = meanlog, col = "grey")

  #estimate CP
  return(temp$CP)
}

## -----------------------------------------------------------------------------
set.seed(114514)
print(paste("The estimated CP is", resulting4(m = 10000)))

## -----------------------------------------------------------------------------
set.seed(114514)
print(paste("The estimated CP is", resulting4(n = 50, meanlog = 2, sdlog = 10, alpha = 0.01)))

## -----------------------------------------------------------------------------
count5test8 = function(x, y) # data analysis
{
  X = x - mean(x)
  Y = y - mean(y)
  outx = sum(X > max(Y)) + sum(X < min(Y))
  outy = sum(Y > max(X)) + sum(Y < min(X)) # return 1 (reject H0) or 0 (do not reject H0)
  reject = as.integer(max(c(outx, outy)) > 5)
  return(reject)
}

Ftest8 = function(x, y, alpha = 0.95)
{
  reject = as.integer(var.test(x, y, alternative = "two.sided", conf.level = 1 - alpha)[[3]] < alpha)
  return(reject)
}

## -----------------------------------------------------------------------------
#test functions above, through a single hypothesis test.
set.seed(114514)
x = rlnorm(1000, 0, 1)
y = rlnorm(1000, 0, 1.2)
c(count5test8(x, y), Ftest8(x, y))
rm(x)
rm(y)

## -----------------------------------------------------------------------------
resulting8 = function(m = 1000, n = 20, sigma1 = sigma1, sigma2 = sigma2, alpha = 0.055) # result reporting
{
  power_c5 = numeric(length(n))
  power_F = numeric(length(n))
  for(j in 1:length(n))
  {
    reject1 = numeric(m)
    reject2 = numeric(m)
    for(i in 1:m)
    {
      x = rnorm(n[j], 0, sigma1)
      y = rnorm(n[j], 0, sigma2)
      reject1[i] = count5test8(x, y)
      reject2[i] = Ftest8(x, y, alpha)
    }
    power_c5[j] = mean(reject1)
    power_F[j] = mean(reject2)
    rm(reject1)
    rm(reject2)
  }
  return(cbind(n, power_c5, power_F))
}

## -----------------------------------------------------------------------------
set.seed(114514)
resulting8(n = c(20, 100, 500, 2000), sigma1 = 1, sigma2 = 1.5)

## -----------------------------------------------------------------------------
set.seed(114514)
sample4 = c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
B4 = 10000
lambdastar4 = numeric(B4)
lambda4 = 1 / mean(sample4)
for(i in 1:B4)
{
  lambdastar4[i] = 1 / mean(sample(sample4, replace=TRUE))
}
hist(lambdastar4, breaks = 20, freq = FALSE)
abline(v = lambda4, col = "red", lwd = 2)
print(c(lambda.bootest = mean(lambdastar4), lambda.MLE = lambda4,
        bias = mean(lambdastar4) - lambda4,
        se.boot = sd(lambdastar4)))

## -----------------------------------------------------------------------------
rm(lambdastar4)

## -----------------------------------------------------------------------------
set.seed(114514)
boot.mean = function(dat, index) {return(mean(dat[index]))}
sampleb5 = boot::boot(data = sample4, statistic = boot.mean, R = 2000)
ci5 = boot::boot.ci(boot.out = sampleb5, type = c("norm", "basic", "perc", "bca"))
ci.norm5 = ci5$norm[2:3]
ci.basic5 = ci5$basic[4:5]
ci.perc5 = ci5$perc[4:5]
ci.BCa5 = ci5$bca[4:5]
ans5 = rbind(ci.norm5, ci.basic5, ci.perc5, ci.BCa5)
rownames(ans5) = c("norm", "basic", "percentile", "BCa")
colnames(ans5) = c("confidence.l", "confidence.h")
print(ans5)

## -----------------------------------------------------------------------------
rm(ci5)
rm(ci.norm5)
rm(ci.basic5)
rm(ci.perc5)
rm(ci.BCa5)
rm(ans5)

## -----------------------------------------------------------------------------
boot.mean = function(dat, index)
{
  return(mean(dat[index]))
}

## -----------------------------------------------------------------------------
bootciA = function(mu, sigma, n = 10, m = 1000)
{
  ci.norm = matrix(NA, m, 2)
  ci.basic = matrix(NA, m, 2)
  ci.perc = matrix(NA, m, 2)

  ci.norm.cover = numeric(m)
  ci.basic.cover = numeric(m)
  ci.perc.cover = numeric(m)

  ci.norm.ml = numeric(m)
  ci.basic.ml = numeric(m)
  ci.perc.ml = numeric(m)

  ci.norm.mr = numeric(m)
  ci.basic.mr = numeric(m)
  ci.perc.mr = numeric(m)

  for(i in 1:m)
  {
    sam = rnorm(n, mu, sigma)
    samb = boot::boot(data = sam, statistic = boot.mean, R = 2000)
    ci = boot::boot.ci(boot.out = samb, type = c("norm", "basic", "perc"))
    ci.norm[i, ] = ci$norm[2:3]
    ci.basic[i, ] = ci$basic[4:5]
    ci.perc[i, ] = ci$percent[4:5]

    mu.sam = mean(sam)
    ci.norm.ml[i] = ci.norm[i, 2] < mu.sam
    ci.norm.mr[i] = ci.norm[i, 1] > mu.sam
    ci.norm.cover[i] = 1 - ci.norm.ml[i] - ci.norm.mr[i]

    ci.basic.ml[i] = ci.basic[i, 2] < mu.sam
    ci.basic.mr[i] = ci.basic[i, 1] > mu.sam
    ci.basic.cover[i] = 1 - ci.basic.ml[i] - ci.basic.mr[i]

    ci.perc.ml[i] = ci.perc[i, 2] < mu.sam
    ci.perc.mr[i] = ci.perc[i, 1] > mu.sam
    ci.perc.cover[i] = 1 - ci.perc.ml[i] - ci.perc.mr[i]

  }
  return(rbind(c(mean(ci.norm.cover), mean(ci.norm.ml), mean(ci.norm.mr)),
               c(mean(ci.basic.cover), mean(ci.basic.ml), mean(ci.basic.mr)),
               c(mean(ci.perc.cover), mean(ci.perc.ml), mean(ci.perc.mr))))
}

## -----------------------------------------------------------------------------
set.seed(114514)
ansA = bootciA(mu = 2, sigma = 2)

## -----------------------------------------------------------------------------
colnames(ansA) = c("coveraged.p", "missed.l", "missed.r")
rownames(ansA) = c("norm", "basic", "perc")
ansA

## -----------------------------------------------------------------------------
rm(ansA)

## -----------------------------------------------------------------------------
principal.jack = function(dat8) # Calculate the estimate propotion of the first principal component of a matrix using jacknife
{
  n8 = dim(dat8)[1]
  p8 = dim(dat8)[2]
  lambda8 = eigen(cov(dat8))$values
  theta8 = lambda8[1] / sum(lambda8) # normal estimate value

  lambda.jack8 = matrix(nrow = n8, ncol = p8)
  theta.jack8 = numeric(n8)
  for(i in 1:n8)
  {
    lambda.jack8[i, ] = eigen(cov(dat8[-i, ]))$values
    theta.jack8[i] = lambda.jack8[i, 1] / sum(lambda.jack8[i, ])
  }

  return(c("theta.hat" = theta8,
           "theta.jack.mean" = mean(theta.jack8),
           "bias.jack" = (n8 - 1) * (mean(theta.jack8) - theta8),
           "sd.jack" = sqrt((n8 - 1) / n8 * sum((theta.jack8 - mean(theta.jack8)) ^ 2))))
}

## -----------------------------------------------------------------------------
dat8 = bootstrap::scor
principal.jack(dat8)

## -----------------------------------------------------------------------------
rm(dat8)

## -----------------------------------------------------------------------------
leaveto = function(dat11) # Using the leave-two-out CV to calculate the MSE of each model above.
{
  n11 = dim(dat11)[1]

  eps11_1 = matrix(nrow = n11 * (n11 - 1) / 2, ncol = 2)
  count = 0
  for(i in 1:(n11-1))
  {
    for(j in (i+1):n11)
    {
      count = count + 1
      fit = lm(y ~ x, data = data.frame(y = dat11$magnetic[-c(i, j)],
                                        x = dat11$chemical[-c(i, j)]))
      eps11_1[count, ] = dat11$magnetic[c(i, j)] - predict(fit,
                                                           newdata = data.frame(x = dat11$chemical[c(i, j)]))
    }
  }
  sigma2.hat11_1 = mean(eps11_1^2)

  eps11_2 = matrix(nrow = n11 * (n11 - 1) / 2, ncol = 2)
  count = 0
  for(i in 1:(n11-1))
  {
    for(j in (i+1):n11)
    {
      count = count + 1
      fit = lm(y ~ x + x2,
               data = data.frame(y = dat11$magnetic[-c(i, j)],
                                 x = dat11$chemical[-c(i, j)],
                                 x2 =  dat11$chemical[-c(i, j)] ^ 2))
      eps11_2[count, ] = dat11$magnetic[c(i, j)] - predict(fit,
                                                           newdata = data.frame(x = dat11$chemical[c(i, j)],
                                                                                x2 = dat11$chemical[c(i, j)] ^ 2))
    }
  }
  sigma2.hat11_2 = mean(eps11_2^2)

  eps11_3 = matrix(nrow = n11 * (n11 - 1) / 2, ncol = 2)
  count = 0
  for(i in 1:(n11-1))
  {
    for(j in (i+1):n11)
    {
      count = count + 1
      fit = lm(y ~ x, data = data.frame(y = log(dat11$magnetic[-c(i, j)]),
                                        x = dat11$chemical[-c(i, j)]))
      eps11_3[count, ] = dat11$magnetic[c(i, j)] - exp(predict(fit,
                                                               newdata = data.frame(x = dat11$chemical[c(i, j)])))
    }
  }
  sigma2.hat11_3 = mean(eps11_3^2)

  eps11_4 = matrix(nrow = n11 * (n11 - 1) / 2, ncol = 2)
  count = 0
  for(i in 1:(n11-1))
  {
    for(j in (i+1):n11)
    {
      count = count + 1
      fit = lm(y ~ x, data = data.frame(y = log(dat11$magnetic[-c(i, j)]),
                                        x = log(dat11$chemical[-c(i, j)])))
      eps11_4[count, ] = dat11$magnetic[c(i, j)] - exp(predict(fit,
                                                               newdata = data.frame(x = log(dat11$chemical[c(i, j)]))))
    }
  }
  sigma2.hat11_4 = mean(eps11_4^2)

  return(matrix(c(sigma2.hat11_1, sigma2.hat11_2, sigma2.hat11_3, sigma2.hat11_4), nrow = 1,
                dimnames = list(c("sigma2.hat"), c("Model I", "Model II", "Model III", "Model IV"))))
}

## -----------------------------------------------------------------------------
dat11 = DAAG::ironslag
leaveto(dat11)

## -----------------------------------------------------------------------------
rm(dat11)

## -----------------------------------------------------------------------------
spearman.per = function(x, y, m = 10000)
{
  rho = cor(x, y, method = "spearman")
  p.value = cor.test(x, y, method = "spearman")$p.value
  n1 = length(x)
  n2 = length(y)

  z = c(x, y)
  rho.per = numeric(m)
  for(i in 1:m)
  {
    z.per = sample(z)
    rho.per[i] = cor(z.per[1:n1], z.per[(n1+1):(n1+n2)], method = "spearman")
  }

  p.value.per = sum(abs(rho) < abs(rho.per)) / m
  return(c(rho = rho, p.value = p.value, p.value.per = p.value.per))
}

## -----------------------------------------------------------------------------
set.seed(114514)
x = 1:10
y = x + rnorm(10, 0, 3)
rbind(x, y)

## -----------------------------------------------------------------------------
set.seed(114514)
spearman.per(x, y)

## -----------------------------------------------------------------------------
rm(x)
rm(y)

## -----------------------------------------------------------------------------
RWMetroplois = function(ini, sigma, m = 10000)
{
  x = numeric(m)
  x[1] = ini
  count = 0
  for(t in 2:m)
  {
    y = rnorm(1, x[t - 1], sigma) #generate
    u = runif(1)
    if(u <= min(exp(abs(x[t - 1]) - abs(y)), 1))
    {
      x[t] = y
    }
    else
    {
      x[t] = x[t - 1]
      count = count + 1
    }
  }
  return(list(x = x, rejectrate = count / m))
}

## -----------------------------------------------------------------------------
set.seed(114514)
x0 = 5
x41 = RWMetroplois(x0, 0.05)
x42 = RWMetroplois(x0, 0.25)
x43 = RWMetroplois(x0, 1)
x44 = RWMetroplois(x0, 5)

## -----------------------------------------------------------------------------
rej4 = cbind(c(0.05, 0.25, 1, 5), c(x41$rejectrate, x42$rejectrate, x43$rejectrate, x44$rejectrate))
colnames(rej4) = c("sigma", "reject.rate")
rej4

## ---- eval=FALSE--------------------------------------------------------------
#  par(mfrow = c(2, 2))
#  plot(x41$x, type = "l", xlab = "sigma = 0.05", ylab = "Xt")
#  plot(x42$x, type = "l", xlab = "sigma = 0.25", ylab = "Xt")
#  plot(x43$x, type = "l", xlab = "sigma = 1", ylab = "Xt")
#  plot(x44$x, type = "l", xlab = "sigma = 5", ylab = "Xt")

## -----------------------------------------------------------------------------
GelmanRubin = function(phi)
{
  k = nrow(phi)
  n = ncol(phi)

  phi.chain.mean = rowMeans(phi)
  Bn = n * var(phi.chain.mean)

  phi.chain.variance = apply(phi, 1, "var")
  Wn = mean((n - 1) / n * phi.chain.variance)

  phi.variance.hat = (n - 1) / n * Wn + 1 / n * Bn

  GR = sqrt(phi.variance.hat / Wn)
  return(GR)
}

## -----------------------------------------------------------------------------
set.seed(114514)
k = 20
m = 10000
sigma = c(0.05, 0.25, 1, 5)
n.converge = numeric(3)
for(i in 1:4)
{
  ini = rnorm(k, 0, 5)
  phi = data.frame()
  for(j in 1:k)
  {
    phi = rbind(phi, RWMetroplois(ini[j], sigma[i], m)$x)
  }
  if(i == 1)
  {
    GR1 = GelmanRubin(as.matrix(phi))
  }
  else
  {
    for(n in 2:m)
    {
      GR = GelmanRubin(as.matrix(phi)[, 1:n])
      if(GR < 1.2)
      {
        n.converge[i-1] = n
        break
      }
    }
  }
}

## -----------------------------------------------------------------------------
GR1

## -----------------------------------------------------------------------------
as.matrix(cbind(sigma = c(0.25, 1, 5), step.converge = n.converge))

## -----------------------------------------------------------------------------
rm(phi)
rm(x41)
rm(x42)
rm(x43)
rm(x44)

## -----------------------------------------------------------------------------
Gibbs = function(ini, mu1, mu2, sigma1, sigma2, rho, n = 10000, burn = 1000)
{
  x = matrix(0, n, 2)
  s1 = sqrt(1 - rho ^ 2) * sigma1
  s2 = sqrt(1 - rho ^ 2) * sigma2

  x[1, ] = ini
  for (i in 2:n)
  {
    x2 = x[i - 1, 2]
    m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
    x[i, 1] = rnorm(1, m1, s1)
    x1 = x[i, 1]
    m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
    x[i, 2] = rnorm(1, m2, s2)
  }

  return(x[(burn + 1):n, ])
}

## -----------------------------------------------------------------------------
set.seed(114514)
x7 = Gibbs(ini = c(5, 5), 0, 0, 1, 1, 0.9)
plot(x7, cex = .5, xlab = "X", ylab = "Y", pch = 20)

## -----------------------------------------------------------------------------
fit7 = lm(y ~ x, data = data.frame(y = x7[, 2], x = x7[, 1]))
summary(fit7)

## -----------------------------------------------------------------------------
qqplot = function(fit)
{
  eps = fit$residuals
  n = length(eps)
  q = qnorm(seq(0, 1, by = 1 / n)[-1] - 0.5 / n)
  plot(x = q, y = sort(eps), pch = 20, xlab = "qnorm", ylab = "residuals")
}
qqplot(fit7)

## -----------------------------------------------------------------------------
set.seed(114514)
k = 20
m = 10000
ini = matrix(c(rnorm(k, -5, 5), rnorm(k, 5, 5)), ncol = 2)
phix = data.frame()
phiy = data.frame()
for(j in 1:k)
{
  temp = Gibbs(ini[j, ], 0, 0, 1, 1, 0.9, burn = 0) # Do not remove.
  phix = rbind(phix, temp[, 1])
  phiy = rbind(phiy, temp[, 2])
}

for(n in 2:m)
{
  GR1 = GelmanRubin(as.matrix(phix)[, 1:n])
  GR2 = GelmanRubin(as.matrix(phiy)[, 1:n])
  if(GR1 < 1.2 & GR2 < 1.2)
  {
    break
  }
}

## -----------------------------------------------------------------------------
n

## -----------------------------------------------------------------------------
rm(x7)
rm(phix)
rm(phiy)

## -----------------------------------------------------------------------------
Sample = function(N = 20, al, be, gam, intercept = c(0, 0)) # data generator
{
  eps = rnorm(2*N, 0, 1)
  X = runif(N)
  M = intercept[1] + al * X + eps[1:N]
  Y = intercept[2] + be * M + gam * X + eps[-(1:N)]
  return(data.frame(cbind(X, M, Y)))
}

## -----------------------------------------------------------------------------
T.stat = function(dat)
{
  colnames(dat) = c("X", "M", "Y")
  fit1 = lm(M ~ X, data = dat)
  fit2 = lm(Y ~ M + X, data = dat)
  T = fit1$coef[2] * fit2$coef[2] / sqrt(fit1$coef[2] ^ 2 *  (summary(fit1)$coef[2,2]) ^ 2 + fit2$coef[2] ^ 2 * (summary(fit2)$coef[2,2]) ^ 2)
  return(T)
  }

## -----------------------------------------------------------------------------
testalpha = function(dat, R, con)
{
  N = nrow(dat)

  T0 = T.stat(dat)

  Tper = numeric(R)
  for(i in 1:R)
  {
    ix = sample(1:N)
    Tper[i] = T.stat(data.frame(cbind(dat$X[ix], dat$M, dat$Y)))
  }

  p = mean(abs(Tper) >= abs(T0))
  return(p < con)
}

## -----------------------------------------------------------------------------
testbeta = function(dat, R, con)
{
  N = nrow(dat)

  T0 = T.stat(dat)

  Tper = numeric(R)
  for(i in 1:R)
  {
    ix = sample(1:N)
    Tper[i] = T.stat(data.frame(cbind(dat$X, dat$M, dat$Y[ix])))
  }

  p = mean(abs(Tper) >= abs(T0))
  return(p < con)
}

## -----------------------------------------------------------------------------
testalphabeta = function(dat, R, con)
{
  N = nrow(dat)

  T0 = T.stat(dat)

  Tper = numeric(R)
  for(i in 1:R)
  {
    ix = sample(1:N)
    Tper[i] = T.stat(data.frame(cbind(dat$X, dat$M[ix], dat$Y)))
  }

  p = mean(abs(Tper) >= abs(T0))
  return(p < con)
}

## -----------------------------------------------------------------------------
TypeIError = function(al, be, gam, R = 59, N = 100, con = 0.05)
{
  I.error1 = numeric(N)
  I.error2 = numeric(N)
  I.error3 = numeric(N)
  for(i in 1:N)
  {
    dat = Sample(al = al, be = be, gam = gam)
    I.error1[i] = testalpha(dat, R = R, con = con)
    I.error2[i] = testbeta(dat, R = R, con = con)
    I.error3[i] = testalphabeta(dat, R = R, con = con)
  }
  return(c(mean(I.error1), mean(I.error2), mean(I.error3)))
}

## -----------------------------------------------------------------------------
set.seed(114514)
TypeIError(0, 0, 1)

## -----------------------------------------------------------------------------
set.seed(114514)
TypeIError(0, 1, 1)

## -----------------------------------------------------------------------------
set.seed(114514)
TypeIError(1, 0, 1)

## -----------------------------------------------------------------------------
Expit = function(N, b1, b2, b3, f0)
{
  X1 = rpois(N, 1)
  X2 = rexp(N, 1)
  X3 = sample(c(0, 1), N, replace = TRUE)
  m = length(f0)
  alpha = numeric(m)
  for(i in 1:m)
  {
    g = function(alpha)
    {
      p = 1 / (1 + exp(- (alpha + b1 * X1 + b2 * X2 + b3 * X3)))
      return(mean(p) - f0[i])
    }
    alpha[i] = uniroot(g, c(-20, 0))$root
  }
  return(alpha)
}

## -----------------------------------------------------------------------------
set.seed(114514)
alpha = Expit(10^6, 0, 1, -1, c(10^-1, 10^-2, 10^-3, 10^-4))
print(alpha)

## -----------------------------------------------------------------------------
plot(x = alpha, y = c(10^-1, 10^-2, 10^-3, 10^-4), pch = 20, ylab = "f0")

## -----------------------------------------------------------------------------
rm(alpha)

## -----------------------------------------------------------------------------
MLE = function(u, v)
{
  loglh = function(lambda)
  {
    return(-sum(u) + sum((v - u) / (exp(lambda * (v - u)) - 1)))
  }
  return(uniroot(loglh, c(0, max(v))))
}

## -----------------------------------------------------------------------------
EM = function(u, v, lambda0 = 1)
{
  iter = 0
  eps = 10 ^ -6
  lambda.next = numeric(1)
  lambda = lambda0
  while(1)
  {
    iter = iter + 1
    lambda.next = 1 / (1 / lambda + mean(u) - mean((v - u) / (exp(lambda * (v - u)) - 1)))
    if(abs(lambda.next - lambda) < eps)
    {
      break
    }
    lambda =lambda.next
  }
  return(c(lambda.hat = lambda.next, iteration = iter))
}

## -----------------------------------------------------------------------------
u = c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
v = c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)

## -----------------------------------------------------------------------------
MLE(u, v)

## -----------------------------------------------------------------------------
EM(u, v)

## -----------------------------------------------------------------------------
rm(u)
rm(v)

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
dat11 = data.frame(matrix(c(1, 2, 3, 4, 5.5, 7), nrow = 3, ncol = 2))
colnames(dat11) = c("X1", "X2")
dat11
sapply(dat11, scale01)

## -----------------------------------------------------------------------------
dat12 = cbind(c("a", "b", "c"), dat11, c("x", "y", "z"))
colnames(dat12) = c("X1", "X2", "X3", "X4")
dat12
sapply(dplyr::select_if(dat12, is.numeric), scale01)

## -----------------------------------------------------------------------------
vapply(dat11, sd, numeric(1))

## -----------------------------------------------------------------------------
vapply(dat12[ ,which(vapply(dat12, is.numeric, logical(1)) == TRUE)], sd, numeric(1))

## -----------------------------------------------------------------------------
rm(dat11)
rm(dat12)

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('NumericMatrix gibbsC(double inix, double iniy, double mu1, double mu2, double sigma1, double sigma2, double rho, int n = 10000, int burn = 1000){

NumericMatrix a(n+burn, 2);
double s1 = sqrt(1 - rho * rho) * sigma1, s2 = sqrt(1 - rho * rho) * sigma2 ;

a(0,0) = inix, a(0,1) = iniy;

for (int i = 1; i <= n+burn-1; i++)
{
  double tempy = a(i-1,1);
  double m1 = mu1 + rho * (tempy - mu2) * sigma1 / sigma2;
  a(i,0) = rnorm(1, m1, s1)[0];

  double tempx = a(i,0);
  double m2 = mu2 + rho * (tempx - mu1) * sigma2/sigma1;
  a(i,1) = rnorm(1, m2, s2)[0];
}

return (a);
}')

## -----------------------------------------------------------------------------
gibbsR = function(ini, mu1, mu2, sigma1, sigma2, rho, n = 10000, burn = 1000)
{
  x = matrix(0, n + burn, 2)
  s1 = sqrt(1 - rho ^ 2) * sigma1
  s2 = sqrt(1 - rho ^ 2) * sigma2

  x[1, ] = ini
  for (i in 2:(n + burn))
  {
    x2 = x[i - 1, 2]
    m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
    x[i, 1] = rnorm(1, m1, s1)
    x1 = x[i, 1]
    m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
    x[i, 2] = rnorm(1, m2, s2)
  }

  return(x)
}

## -----------------------------------------------------------------------------
set.seed(114514)
C = gibbsC(100, 10, 0, 0, 1, 1, 0.9)[-(1:1000), ]
R = gibbsR(c(100, 10), 0, 0, 1, 1, 0.9)[-(1:1000), ]

## -----------------------------------------------------------------------------
plot(sort(C[,1]), sort(R[,1]), main="qqplot of x1", xlab="XC", ylab="XR")
plot(sort(C[,2]), sort(R[,2]), main="qqplot of x2", xlab="XC", ylab="XR")

## -----------------------------------------------------------------------------
ts = microbenchmark::microbenchmark(gibbC = gibbsC(100, 10, 0, 0, 1, 1, 0.9), gibbR = gibbsR(c(100, 10), 0, 0, 1, 1, 0.9))
summary(ts)[, c(1, 3, 5, 6)]

## -----------------------------------------------------------------------------
rm(C)
rm(R)

