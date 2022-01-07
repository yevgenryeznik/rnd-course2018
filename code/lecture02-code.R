library(tidyverse)
library(nleqslv)

# D-optimal allocation
doptimal <- function(sigma_vec) {
  # K -- number of treatment groups
  K <- length(sigma_vec)
  
  rep(1, K)/K
}

# A-optimal allocation
aoptimal <- function(sigma_vec) {
  sigma_vec/sum(sigma_vec)
}


# D-A-optimal allocation
daoptimal <- function(sigma_vec) {
  obj_fcn <- function(rho, sigma) {
    rho_ <- rho^2/sum(rho^2)
    one <- cbind(rep(1, length(sigma)-1))
    det(diag(sigma[-1]^2/rho_[-1])+sigma[1]^2/rho_[1]*one%*%t(one))
  }
  rho <- optim(rep(1, length(sigma_vec)), obj_fcn, sigma = sigma_vec)$par
  rho^2/sum(rho^2)
}


# A-A-optimal allocation
aaoptimal <- function(sigma_vec) {
  # K -- number of treatment groups
  K <- length(sigma_vec)

  fctr <- c(sqrt(K-1), rep(1, K-1))
  fctr*sigma_vec/sum(fctr*sigma_vec)
}

# D-efficiency
defficiency <- function(rho1, rho2, sigma_vec) {
  # K -- number of treatment groups
  K <- length(sigma_vec)
  
  obj_fcn <- map_dbl(list(rho1, rho2), ~ prod(sigma_vec^2/.))

  (obj_fcn[1]/obj_fcn[2])^(1/K)
}


# A-efficiency
aefficiency <- function(rho1, rho2, sigma_vec) {
  # K -- number of treatment groups
  K <- length(sigma_vec)
  
  obj_fcn <- map_dbl(list(rho1, rho2), ~ sum(sigma_vec^2/.))

  (obj_fcn[1]/obj_fcn[2])
}


# D-A-efficiency
daefficiency <- function(rho1, rho2, sigma_vec) {
  # K -- number of treatment groups
  K <- length(sigma_vec)
  one <- cbind(rep(1, K-1))
  
  obj_fcn <- map_dbl(list(rho1, rho2), ~ det(diag(sigma_vec[-1]^2/.[-1]) + 
                                               sigma_vec[1]^2/.[1]*one%*%t(one)))

  (obj_fcn[1]/obj_fcn[2])^(1/K)
}


# A-A-efficiency
aaefficiency <- function(rho1, rho2, sigma_vec) {
  # K -- number of treatment groups
  K <- length(sigma_vec)
  one <- cbind(rep(1, K-1))
  
  obj_fcn <- map_dbl(list(rho1, rho2), ~ sum(sigma_vec[-1]^2/.[-1]) + (K-1)*sigma_vec[1]^2/.[1])
  
  (obj_fcn[1]/obj_fcn[2])
}


# Fixed trace optimal allocation (ftoptimal)
ftoptimal <- function(prob_vec) {
  # number of treatments
  K <- length(prob_vec)
  
  p <- prob_vec*c((K-1), rep(1, K-1))
  
  sqrt(p)/sum(sqrt(p))
}


# simulate data to fit logistic regression
# d <- c(-1.5, 0, 1.5)
# alpha  <- 1
# beta <- 1
simulate_logistic <- function(nsbj, d, alpha, beta) {
  # doses assigns to subjetcs
  dose <- sample(d, nsbj, TRUE)
  
  # number of subjects per dose
  n <- map_dbl(d, ~ sum(as.numeric(. == dose)))
  
  # responses
  Y <- rbinom(length(dose), 1, prob = 1/(1+exp(-alpha-beta*dose))) 
  
  # number of events per dose
  x <- map_dbl(d, ~ sum(Y[. == dose]))
  
  data_frame(d = d, n = n, x = x)
}

# fit 2-parameter logistic regression
fit_logistic <- function(d, n, x) {
  # define score equations
  score_equations <- function(theta, d, n, x) {
    equation <- numeric(2)
    
    alpha <- theta[1]
    beta <- theta[2]
    
    prob <- 1/(1+exp(-alpha-beta*d))
    
    equation[1] <- sum(x - n*prob)
    equation[2] <- sum((x - n*prob)*d)
    
    return(equation)
  }
  
  theta_mle <- nleqslv(c(0, 0), function(theta) {score_equations(theta, d, n, x)}, control=list(btol=.001))$x
  list(alpha = theta_mle[1], beta = theta_mle[2])
}

# simulation
# nsim <- 1000
# map(seq_len(nsim), ~ {
#   obs <- simulate_logistic(100, d, alpha, beta)
#   coeffs <- fit_logistic(obs$d, obs$n, obs$x)
#   
#   data_frame(alpha = coeffs$alpha, beta = coeffs$beta)
# }) %>% 
#   bind_rows() %>% 
#   gather(coefficient, value) %>% 
#   group_by(coefficient) %>% 
#   summarise(mean = mean(value), sd = sd(value))
  
  
