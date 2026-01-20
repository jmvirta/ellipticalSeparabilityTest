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



# Auxiliary function that estimates the two constants (m4 + m2)/betaF^2 and
# (3m4 - m2)/betaF^2 required in Theorems 12, 13 in the paper
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




# Auxiliary function that computes the projection matrix parts of
# the pseudo-inverse of \Omega^{1/2} B \Omega^{1/2}
compute_pseudoinverse_parts <- function(p1, p2){
  
  K1 <- matrix(0, p1^2*p2^2, p1^2*p2^2)
  for(i in 1:p1){
    for(j in 1:p1){
      Eij <- matrix(0, p1, p1)
      Eij[i, j] <- 1
      K1 <- K1 + diag(p2)%x%Eij%x%diag(p2)%x%t(Eij)
    }
  }
  
  K2 <- matrix(0, p1^2*p2^2, p1^2*p2^2)
  for(i in 1:p2){
    for(j in 1:p2){
      Eij <- matrix(0, p2, p2)
      Eij[i, j] <- 1
      K2 <- K2 + Eij%x%diag(p1)%x%t(Eij)%x%diag(p1)
    }
  }
  
  V1 <- matrix(0, p1^2*p2^2, p1^2*p2^2)
  for(i in 1:p1){
    for(j in 1:p1){
      Eij <- matrix(0, p1, p1)
      Eij[i, j] <- 1
      V1 <- V1 + diag(p2)%x%Eij%x%diag(p2)%x%Eij
    }
  }
  
  V2 <- matrix(0, p1^2*p2^2, p1^2*p2^2)
  for(i in 1:p2){
    for(j in 1:p2){
      Eij <- matrix(0, p2, p2)
      Eij[i, j] <- 1
      V2 <- V2 + Eij%x%diag(p1)%x%Eij%x%diag(p1)
    }
  }
  
  G1 <- 0.25*(diag(p1^2*p2^2) + K1 + K2 + K1%*%K2)
  G2 <- 0.25*(diag(p1^2*p2^2) - K1 - K2 + K1%*%K2)
  
  
  p <- p1*p2
  Q1 <- diag(p1^2) - (1/p1)*c(diag(p1))%*%t(c(diag(p1)))
  R1 <- (1/p2)*Q1%*%(t(c(diag(p2)))%x%diag(p1^2))%*%(diag(p2)%x%commutation.matrix(p2, p1)%x%diag(p1))
  
  Q2 <- diag(p2^2) - (1/p2)*c(diag(p2))%*%t(c(diag(p2)))
  R2 <- (1/p1)*Q2%*%(t(c(diag(p1)))%x%diag(p2^2))%*%(diag(p1)%x%commutation.matrix(p1, p2)%x%diag(p2))%*%(commutation.matrix(p1, p2)%x%commutation.matrix(p1, p2))
  
  Qp <- diag(p^2) - (1/p)*c(diag(p))%*%t(c(diag(p)))
  
  B0 <- ((diag(p2) %x% commutation.matrix(p1, p2) %x% diag(p1))%*%((R2 %x% c(diag(p1))) + (c(diag(p2)) %x% R1)) - Qp)
  
  return(list(E1 = B0%*%G1%*%t(B0), E2 = B0%*%G2%*%t(B0)))
}







# Performs the asymptotic test of covariance separability for elliptical data
# using the squared norm test statistic and returns the p-value for
# H0 = Sigma is separable
#
# x        = the data as an array of size p1 x p2 x n
# gaussian = whether the data is Gaussian,
#           if TRUE, a simpler version of the limiting distribution is used
ellip_norm_sep_test <- function(x, gaussian = FALSE){
  
  
  p1 <- dim(x)[1]
  p2 <- dim(x)[2]
  n <- dim(x)[3]
  
  # Estimate the constants required in the limiting distribution
  t_const <- estim_m_constants(x)
  
  if(t_const[2] < 0){
    t_const[2] <- 0
  }

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







# Performs the asymptotic test of covariance separability for elliptical data
# using the Wald test statistic and returns the p-value for
# H0 = Sigma is separable
#
# x        = the data as an array of size p1 x p2 x n
# gaussian = whether the data is Gaussian,
#           if TRUE, the fourth moments are not estimated
ellip_wald_sep_test <- function(x, gaussian = FALSE){
  
  p1 <- dim(x)[1]
  p2 <- dim(x)[2]
  n <- dim(x)[3]
  
  # Estimate the constants required in the limiting distribution
  if(!gaussian){
    t_const <- estim_m_constants(x)
    if(t_const[2] < 0){
      t_const[2] <- 0
    }
  } else {
    t_const <- c(2, 2)
  }

  # Pseudoinverse formation
  orth_parts <- compute_pseudoinverse_parts(p1, p2)
  if(t_const[2] > 1e-6){
    pseudo_i <- (1/t_const[1])*orth_parts$E1 + (1/t_const[2])*orth_parts$E2
  } else {
    pseudo_i <- (1/t_const[1])*orth_parts$E1
  }
  
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
  V_diff_mat <- R_invhalf%*%(R2%x%R1)%*%R_invhalf - diag(p1*p2)

  test_stat <- c(n*t(c(V_diff_mat))%*%pseudo_i%*%c(V_diff_mat))
  
  pval <- pchisq(test_stat, df = 0.25*(p1 + 2)*(p1 - 1)*(p2 + 2)*(p2 - 1) + 0.25*p1*(p1 - 1)*p2*(p2 - 1), lower.tail = FALSE)
  
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
ellip_norm_sep_test(x)
ellip_wald_sep_test(x)
lrt_sep_test(x)


# Non-separable non-Gaussian data
# LRT has better power (at cost of wrong level)

x <- rmatrixt(n = n, df = nu, mean = matrix(0, p1, p2))*sqrt(nu - 2)
x[1, 1, ] <- 1.2*x[1, 1, ]
ellip_norm_sep_test(x)
ellip_wald_sep_test(x)
lrt_sep_test(x)


# Separable Gaussian data
# Similar performance

x <- rmatrixnorm(n, mean = matrix(0, p1, p2))
ellip_norm_sep_test(x)
ellip_wald_sep_test(x)
lrt_sep_test(x)


# Non-separable Gaussian data
# Similar performance

x <- rmatrixnorm(n, mean = matrix(0, p1, p2))
x[1, 1, ] <- 1.2*x[1, 1, ]
ellip_norm_sep_test(x)
ellip_wald_sep_test(x)
lrt_sep_test(x)
