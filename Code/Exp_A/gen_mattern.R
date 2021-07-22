rho <- 3
sigma <- 1
nu <- 1/3
K <- 4

mattern.cov <- function(d, nu=5/2, rho=1, sigma=1) {
  tmp <- sqrt(2*nu)*d/rho
  a <- besselK(x=tmp, nu=nu) #, expon.scaled = FALSE)
  b <- tmp^nu
  return( sigma^2 * 2^(1-nu) * b * a / gamma(nu) )
}


compute_mattern_eigen <- function(nu, rho, sigma, K){
  # "compute" eigenfunctions numerically
  x <- seq(0, 1, length=500)
  mm <- matrix(NA, length(x), length(x))
  for(i in 1:length(x)) mm[i,] <- mattern.cov(abs(x[i]-x), nu=nu, 
                                              rho= rho, sigma=sigma)
  diag(mm) <- sigma^2
  mm.svd <- svd(mm)$u
  delta <- diff(x)[1]
  phis <- vector('list', K)
  for(j in 1:K){ 
    phis[[j]] <- approxfun(x, mm.svd[,j]/sqrt(delta))
  }
  
  return(phis)
}


generate_mattern <- function(mu, phis, K){
  
  # eigenvalues of the covariance function
  lambdas <- c(0.8, 0.3, 0.2, 0.1)
  # eigenfunctions
  
  fun <- function(xi_matrix, tt){
    xx <- matrix(1,nrow = nrow(xi_matrix), ncol = 1)%*%mu(tt)
    for(j in 1:K){
      xx <- xx + sqrt(lambdas[j]) *as.matrix(xi_matrix[,j])%*%(phis[[j]](tt))
    }
    return(xx)
  }
  return(fun)
}



mu <- function(tt){
  2*sin(tt*pi)*exp(1-tt)
}

phis <- compute_mattern_eigen(nu, rho, sigma, K)

library(rlist)
saveRDS(list(mu = mu, phis = phis), file = "Code/Exp_A/mattern.rds")

