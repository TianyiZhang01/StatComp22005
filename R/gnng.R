#' @title The group LARS method to select factors in liner regression problem raised by Yuan and Lin(2004)
#' @description The group non-negetive garrotte method to select factors in liner regression problem raised by Yuan and Lin(2004)
#' @param X the design matrix. Some details about X are shown in group
#' @param Y the response variable
#' @param sigma the standard deviation of the noise term to calculate Cp statistic of selected models. If sigma = -1 the OLS estimate deviation will be used in the progress. Initial -1
#' @param group a numeric vector of the quantity of variables in each factor(group). Thus the variables in one factor are asked to be arranged contiguously in the design matrix X
#' @param scale logical. If TRUE, each variable is standardized to have unit L2 norm, otherwise it is left alone.
#' @return a list abput the solution path, including the selectived coefficient, the Cp statistic of each model and the time the method cost.
#' @examples
#' \dontrun{
#'   Y = coef %*% X + rnorm(n, 0, 0.5)
#'   spath = gnng(X, Y, sigma = 0.5, group = rep(2, 15))
#'  }
#' @export

#Group non-negative Garrotte
gnng = function(X, Y, sigma = -1, group, scale = FALSE)
{

  st = Sys.time()
  #计时开始

  n = dim(X)[1]
  M = dim(X)[2]
  J = length(group)
  d = rep(0, J)
  k = 0
  r = Y
  activeset = numeric(0)
  rstorage = matrix(nrow = n, ncol = 0)
  dstorage = matrix(nrow = J, ncol = 0)
  betastorage = matrix(nrow = M, ncol = 0)

  I = numeric(J)
  I[1] = 0
  flag = 0
  for(j in 1:J)
  {
    flag = flag + group[j]
    I[j+1] = flag
  }
  #初始化

  X = scale(X, scale = scale)
  Y = scale(Y, scale = scale)
  #中心化、正交化

  betals = MASS::ginv(t(X) %*% X) %*% t(X) %*% Y
  if(sigma == -1)
  {
    res = Y - X %*% betals
    sigma = sum(res ^ 2) / (n - M)
  }
  Z = matrix(0, nrow = n, ncol = J)

  for(j in 1:J)
  {
    Z[, j] = X[, (I[j]+1):I[j+1]] %*% betals[(I[j]+1):I[j+1]]
  }
  #OLS

  temp = 0
  index = 0
  for(i in 1:J)
  {
    corr = numeric(J)
    corr[i] = t(Z[, i]) %*% r / group[i]
    if(corr[i] > temp)
    {
      temp = corr[i]
      index = i
    }
  }
  activeset = sort(c(activeset, index))
  #初始化

  while(1)
  {
    direc = MASS::ginv(t(Z[, activeset]) %*% Z[, activeset]) %*% t(Z[, activeset]) %*% r
    direc0 = numeric(J)
    l = 0
    for(i in activeset)
    {
      l = l + 1
      direc0[i] = direc[l]
    }
    move = Z %*% direc0
    #选择方向

    alpha = rep(1, J)
    for(j in (1:J)[-activeset])
    {
      f = function(x)
      {
        return(t(Z[, j]) %*% (r - x * move) / group[j] - t(Z[, activeset[1]]) %*% (r - x * move) / group[activeset[1]])
      }
      if(f(-0.1) * f(1.1) > 0)
      {
        alpha[j] = 1
      } else
      {
        alpha[j] = stats::uniroot(f, c(-0.1, 1.1))$root
      }
    }
    for(j in activeset)
    {
      betaj = -d[j] / MASS::ginv(t(Z[, j]) %*% Z[, j]) %*% t(Z[, j]) %*% r
      alpha[j] = min(betaj, 1)
    }

    if(max(alpha) <= 0) alpha = rep(1, J)
    alphamin = min(alpha[which(alpha > 0)])
    jstar = which(alpha == alphamin)
    d = d + alphamin * direc0

    if(length(activeset) > 0 & length(intersect(jstar, activeset)) == length(jstar))
    {
      activeset = sort(setdiff(activeset, jstar))
    } else
    {
      activeset = sort(union(activeset, jstar))
    }
    #选择距离

    r = Y - Z %*% d
    k = k+1
    #更新

    dstorage = cbind(dstorage, d)
    if(alphamin == 1) break
    #储存路径、检查是否完成计算
  }

  dstorage = cbind(rep(0, J), dstorage)
  #补上空模型

  rstorage = Z %*% dstorage
  for(i in 1:dim(dstorage)[2])
  {
    rstorage[, i] = Y - rstorage[, i]
    dbeta = c()
    for(j in 1:J)
    {
      dbeta = c(dbeta, dstorage[j, i] * betals[(I[j]+1):I[j+1]])
    }
    betastorage = cbind(betastorage, dbeta)
  }
  #计算残差、β

  df = 0
  K = dim(dstorage)[2]
  for(i in 1:K)
  {
    temp = dstorage[, i]
    dsl = numeric(J)
    d = numeric(J)
    for(j in 1:J)
    {
      dsl[j] = 1 - isTRUE(all.equal(temp[j], 0))
      d[j] = temp[j] * (group[j] - 2)
    }
    df = c(df, sum(dsl) + sum(d))
  }

  Cp = numeric(J+1)
  for(i in 1:(J+1))
  {
    Cp[i] = sum(rstorage[, i] ^ 2) / sigma ^ 2 - n + 2 * df[i]
  }
  #计算Cp值

  rate = numeric(K)
  b0 = sum(abs(betals))
  for(j in 1:K)
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
