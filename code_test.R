# Asymptotic test of separability for matrix elliptical data

# For generating matrix-t data
library(MixMatrix)

# For approximating the CDF of sums of scaled chi-squares
library(CompQuadForm)

# For general matrix data operations
library(tensorBSS)


# Generates a matrix-t sample
# scaled such that E(X_11^2) = 1
gen_matrix_t <- function(n, p1, p2, nu){
  x <- rmatrixt(n = n, df = nu, mean = matrix(0, p1, p2))*sqrt(nu - 2)
}



# Estimates the two constants (m4 + m2) and (3m4 - m2) required in Theorem 12 in the paper
estim_m_constants <- function(x){
  p1 <- dim(x)[1]
  p2 <- dim(x)[2]
  n <- dim(x)[3]
  est_2 <- MLmatrixnorm(x, tol = 1e-8)
  S1 <- est_2$U
  S2 <- est_2$V*est_2$var
  
  eig_S1 <- eigen(S1)
  S1_invsqrt <- eig_S1$vectors%*%diag(eig_S1$values^(-1/2))%*%t(eig_S1$vectors)
  
  eig_S2 <- eigen(S2)
  S2_invsqrt <- eig_S2$vectors%*%diag(eig_S2$values^(-1/2))%*%t(eig_S2$vectors)
  
  x_stand <- tensorCentering(x)
  x_stand <- tensorTransform2(x_stand, list(S1_invsqrt, S2_invsqrt), c(1, 2))
  
  d1 <- mean(x_stand^2)
  
  d2 <- 0
  for(i in 1:n){
    d2 <- d2 + (sum(x_stand[, , i]^2))^2
  }
  d2 <- d2/(n*p1*p2)
  
  d3 <- mean(x_stand^4)
  
  t1 <- (d2 + (p1*p2 - 2*p1 - 2*p2)*d3/3)/((p1 - 1)*(p2 - 1)*d1^2)
  t2 <- (3*d2 - (p1 + 2)*(p2 + 2)*d3/3)/((p1 - 1)*(p2 - 1)*d1^2)
  
  c(t1, t2)
}




# Performs the asymptotic test of covariance separability for elliptical data
# and returns the p-value for
# H0 = Sigma is separable
#
# x        = the data as an array of size p1 x p2 x n
# gaussian = whether the data is Gaussian,
#           if TRUE, a simpler version of the limiting distribution is used
ellip_sep_test <- function(x, gaussian = FALSE){
  
  
  p1 <- dim(x)[1]
  p2 <- dim(x)[2]
  n <- dim(x)[3]
  
  # Estimate the constants required in the limiting distribution
  t_const <- estim_m_constants(x)

  # Full vec-covariance and its scaled version
  S <- cov(t(apply(x, 3, c)))
  R <- S/(det(S)^(1/(p1*p2)))
  
  # Matrix normal MLE and its scaled versions
  est_2 <- MLmatrixnorm(x, tol = 1e-8)
  S1 <- est_2$U
  S2 <- est_2$V*est_2$var
  R1 <- S1/(det(S1)^(1/p1))
  R2 <- S2/(det(S2)^(1/p2))
  
  # Test statistic computation
  eig_R <- eigen(R)
  R_invhalf <- eig_R$vectors%*%diag(eig_R$values^(-1/2))%*%t(eig_R$vectors)
  
  test_stat <- n*norm(R_invhalf%*%(R2%x%R1)%*%R_invhalf - diag(p1*p2), type = "F")^2
  
  if(gaussian){
    pval <- pchisq(test_stat/2, df = 0.25*(p1 + 2)*(p1 - 1)*(p2 + 2)*(p2 - 1) + 0.25*p1*(p1 - 1)*p2*(p2 - 1), lower.tail = FALSE)
  } else {
    pval <- suppressWarnings(imhof(test_stat, lambda = t_const, h = c(0.25*(p1 + 2)*(p1 - 1)*(p2 + 2)*(p2 - 1), 0.25*p1*(p1 - 1)*p2*(p2 - 1)))$Qq)
  }
  
  pval
}




# Performs the Gaussian LRT of covariance separability for elliptical data
# and returns the p-value for
# H0 = Sigma is separable
#
# x = the data as an array of size p1 x p2 x n
lrt_sep_test <- function(x){
  p1 <- dim(x)[1]
  p2 <- dim(x)[2]
  n <- dim(x)[3]
  
  # Full vec-covariance
  S <- cov(t(apply(x, 3, c)))
  
  # Matrix normal MLE
  est_2 <- MLmatrixnorm(x, tol = 1e-8)
  S1 <- est_2$U
  S2 <- est_2$V*est_2$var
  
  # Test statistic computation
  test_stat <- n*(log(det(S2%x%S1)) - log(det(S)) + sum(S*(solve(S2)%x%solve(S1))) - p1*p2)
  
  pval <- pchisq(test_stat, df = 0.5*p1*p2*(p1*p2 + 1) - 0.5*p1*(p1 + 1) - 0.5*p2*(p2 + 1) + 1, lower.tail = FALSE)
  
  pval
}



# Examples

n <- 1000
p1 <- 3
p2 <- 3
nu <- 5


# Separable non-Gaussian data
# LRT rejects too often

x <- rmatrixt(n = n, df = nu, mean = matrix(0, p1, p2))*sqrt(nu - 2)
ellip_sep_test(x)
lrt_sep_test(x)


# Non-separable non-Gaussian data
# LRT has better power (at cost of wrong level)

x <- rmatrixt(n = n, df = nu, mean = matrix(0, p1, p2))*sqrt(nu - 2)
x[1, 1, ] <- 1.2*x[1, 1, ]
ellip_sep_test(x)
lrt_sep_test(x)


# Separable Gaussian data
# Similar performance

x <- rmatrixnorm(n, mean = matrix(0, p1, p2))
ellip_sep_test(x)
lrt_sep_test(x)


# Non-separable Gaussian data
# Similar performance

x <- rmatrixnorm(n, mean = matrix(0, p1, p2))
x[1, 1, ] <- 1.2*x[1, 1, ]
ellip_sep_test(x)
lrt_sep_test(x)



