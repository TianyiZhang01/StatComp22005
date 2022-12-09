#' @title The group LARS method to select factors in liner regression problem raised by Yuan and Lin(2004)
#' @description The group LARS method to select factors in liner regression problem raised by Yuan and Lin(2004)
#' @param X the design matrix. Some details about X are shown in group
#' @param Y the response variable
#' @param sigma the standard deviation of the noise term to calculate Cp statistic of selected models. If sigma = -1 the OLS estimate deviation will be used in the progress. Initial -1
#' @param group a numeric vector of the quantity of variables in each factor(group). Thus the variables in one factor are asked to be arranged contiguously in the design matrix X
#' @param scale logical. If TRUE, each variable is standardized to have unit L2 norm, otherwise it is left alone.
#' @return a list abput the solution path, including the selectived coefficient, the Cp statistic of each model and the time the method cost.
#' @examples
#' \dontrun{
#'   Y = coef %*% X + rnorm(n, 0, 0.5)
#'   spath = gLARS(X, Y, sigma = 0.5, group = rep(2, 15))
#'  }
#' @export

#Group LARS
gLARS = function(X, Y, sigma = -1, group, scale = FALSE)
{

  st = Sys.time()
  #计时开始

  X = scale(X, scale = scale)
  Y = scale(Y, scale = scale)
  #中心化

  n = dim(X)[1]
  M = dim(X)[2]
  J = length(group)
  beta = 0
  k = 0
  r = Y
  activeset = numeric(0)
  Xacts = numeric(0)
  rstorage = matrix(0, nrow = n, ncol = J+1)
  betastorage = matrix(0, nrow = M, ncol = J)

  betals = MASS::ginv(t(X) %*% X) %*% t(X) %*% Y
  if(sigma == -1)
  {
    res = Y - X %*% betals
    sigma = sum(res ^ 2) / (n - M)
  }

  I = numeric(J)
  I[1] = 0
  flag = 0
  for(j in 1:J)
  {
    flag = flag + group[j]
    I[j+1] = flag
  }
  #初始化

  temp = 0
  index = 0
  for(i in 1:J)
  {
    corr = numeric(J)
    corr[i] = sum((t(X[, (I[i] + 1):I[i + 1]]) %*% r) ^ 2) / group[i]
    if(corr[i] > temp)
    {
      temp = corr[i]
      index = i
    }
  }
  activeset = sort(c(activeset, index))
  Xacts = sort(c(Xacts, (I[index] + 1):I[index + 1]))
  #初始化

  while(1)
  {
    direc = MASS::ginv(t(X[, Xacts]) %*% X[, Xacts]) %*% t(X[, Xacts]) %*% r
    direc0 = numeric(M)
    l = 0
    for(i in Xacts)
    {
      l = l + 1
      direc0[i] = direc[l]
    }
    move = X %*% direc0
    #选择方向

    alpha = rep(1, J)
    for(j in (1:J)[-activeset])
    {
      f = function(x)
      {
        return(sum((t(X[, (I[j] + 1):I[j + 1]]) %*% (r - x * move)) ^ 2 / group[j]) - sum((t(X[, (I[activeset[1]] + 1):I[activeset[1] + 1]]) %*% (r - x * move)) ^ 2 / group[activeset[1]]))
      }
      if(f(-0.1) * f(1) > 0)
      {
        alpha[j] = 1
      } else
      {
        alpha[j] = stats::uniroot(f, c(-0.1, 1))$root
      }
    }
    alphamin = min(alpha)
    jstar = which(alpha == alphamin)
    activeset = sort(union(activeset, jstar))
    Xacts = sort(union(Xacts, (I[jstar[1]] + 1):I[jstar[1] + 1]))
    #选择距离

    beta = beta + alpha[jstar[1]] * direc0
    r = Y - X %*% beta
    k = k+1
    #更新

    betastorage[, k] = beta
    if(alphamin == 1) break
    #储存路径、检查是否完成计算
  }

  betastorage = cbind(rep(0, M), betastorage)
  #补上空模型

  rstorage = X %*% betastorage
  for(i in 1:(J+1))
  {
    rstorage[, i] = Y - rstorage[, i]
  }
  #计算残差

  df = 0
  betagr = matrix(0, nrow = J, ncol = J)
  for(i in 1:J)
  {
    for(j in 1:J)
    {
      betagr[j, i] = sqrt(sum((betastorage[(I[j] + 1):I[j + 1], i+1] - betastorage[(I[j] + 1):I[j + 1], i])^2))
    }
  }
  for(i in 1:J)
  {
    temp = betastorage[, i]
    betasl = numeric(J)
    betagrr = numeric(J)
    for(j in 1:J)
    {
      betasl[j] = 1 - isTRUE(all.equal(sum(temp[(I[j] + 1):I[j + 1]] ^ 2), 0))
      betagrr[j] = sum(betagr[j, 1:i]) / sum(betagr[j, ]) * (group[j] - 1)
    }
    df = c(df, sum(betasl) + sum(betagrr))
  }

  Cp = numeric(J+1)
  for(i in 1:(J+1))
  {
    Cp[i] = sum(rstorage[, i] ^ 2) / sigma ^ 2 - n + 2 * df[i]
  }
  #计算Cp值

  rate = numeric(J+1)
  b0 = sum(abs(betals))
  for(j in 1:(J+1))
  {
    rate[j] = sum(abs(betastorage[, j])) / b0
  }
  #计算长度


  #Rinverse = MASS::ginv(qr.R(qr))
  #betastorage = Rinverse %*% betastorage
  #正交化还原

  fi = Sys.time()
  t = fi - st
  #计时结束

  return(list(beta = betastorage, residuals = rstorage, Cp = Cp, t = t, rate = rate))
}
