## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  #生成正态随机变量矩阵
#  l <- 0.2; n <- 50; p <- 30
#  theta <- c(rep(log(0.5*p),5),rep(0,p-5))
#  sigma <- matrix(0,nrow = p,ncol = p)
#  for (i in 1:p) {
#    for (j in 1:p) {
#      sigma[i,j] <- l^abs(i-j)
#    }
#  }
#  w <- MASS::mvrnorm(n,mu=theta,Sigma = sigma)
#  
#  #将正态样本转化为对数成分数据z
#  x <- w
#  for (i in 1:n) {
#    x[i,] <- exp(w[i,])/sum(exp(w[i,]))
#  }
#  z <- log(x)
#  
#  #生成响应变量y
#  sig <- 0.5; beta_star <- c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2,rep(0,p-8))
#  y <- z %*% beta_star+rnorm(n,mean = 0, sd = sig)

## ----1------------------------------------------------------------------------
library(SA24204132)
data(data)
y <- data[1]
z <- data[2:31]

## ----2------------------------------------------------------------------------
mu <- 1
lambda1 <- r_clambda(y,z,mu)
lambda1

## ----3------------------------------------------------------------------------
mu <- 1
lambda2 <- r_slambda(y,z,mu,low = 0, up = 0.15, by = 0.01)
lambda2

## ----4------------------------------------------------------------------------
mu <- 1; lambda <- lambda1
r_cdmlasso(y,z,lambda,mu)

## ----6------------------------------------------------------------------------
set.seed(132); options(warn = -1); simulations <- 10
#设定要生成的n个p维正态分布随机变量的信息
l <- 0.2; n <- 50; p <- 30
theta <- c(rep(log(0.5*p),5),rep(0,p-5))
sigma <- matrix(0,nrow = p,ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    sigma[i,j] <- l^abs(i-j)
  }
}
#m次mc随机实验
for (m in 1:simulations) {
  mc_beta <- matrix(0,nrow = m, ncol = p)
  mc_PE <- mc_l1 <- mc_l2 <- mc_linf <- mc_FP <- mc_FN <- numeric(m)
  w <- MASS::mvrnorm(n,mu=theta,Sigma = sigma)
  #将w样本转化为成分数据x
  x <- w
  for (i in 1:n) {
    x[i,] <- exp(w[i,])/sum(exp(w[i,]))
  }
  #按照模型1、2方式转换成分数据x为对数数据z
  z <- log(x)
  #生成响应变量y
  mu <- 1; sig <- 0.5; beta_star <- c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2,rep(0,p-8))
  y <- z %*% beta_star+rnorm(n,mean = 0, sd = sig)
  lambda <- c_clambda(y,z,mu)
  mc_beta[m,] <- c_cdmlasso(y,z,lambda,mu)
  
  #计算measuring指标
  bias <- mc_beta[m,]-beta_star
  mc_PE[m] <- (norm((y-(z %*% mc_beta[m,])),type = '2')^2)/n
  mc_l1[m] <- sum(abs(bias))
  mc_l2[m] <- norm(bias,type = '2')
  mc_linf[m] <- max(abs(bias))
  mc_FP[m] <- sum(beta_star[which(mc_beta[m,] != 0)] == 0)
  mc_FN[m] <- sum(mc_beta[m,][which(beta_star != 0)] == 0)
}
print(round(c('rho'=l,
              'samplesize'=n,
              'dimension'=p,
              'PE'=mean(mc_PE),
              'l1loss'=mean(mc_l1),
              'l2loss'=mean(mc_l2),
              'linfloss'=mean(mc_linf),
              'FP'=mean(mc_FP),
              'FN'=mean(mc_FN)),2))

