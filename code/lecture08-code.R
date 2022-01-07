library(tidyverse)
library(nleqslv)


# simulate data to fit logistic regression
# example:
#  d <- c(5, 10, 30, 50, 75)
#  alpha <- 0.05
#  beta <- 0.005
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
  
  ## MLE's are the solutions of the score equations
  theta_mle <- nleqslv(c(0, 0), function(theta) {score_equations(theta, d, n, x)}, control=list(btol=.001))$x
  alpha_hat = theta_mle[1] 
  beta_hat = theta_mle[2]
  
  ## Fisher Information matrix
  psi <- function(d, alpha, beta) {
    exp(-(alpha+beta*d))/(1+exp(-(alpha+beta*d)))^2
  }
  
  FIM <- rbind(
    c(sum(n*psi(d, alpha_hat, beta_hat)), sum(n*d*psi(d, alpha_hat, beta_hat))),
    c(sum(n*d*psi(d, alpha_hat, beta_hat)), sum(n*d^2*psi(d, alpha_hat, beta_hat)))
  )
  
  ## Variance-Covariance Matrix
  Cov = solve(FIM)
  
  ## P(d): function to estimate a probability of toxicity given MLE's (alpha_hat, beta_hat)
  pr_hat = function(dose) {
    1/(1+exp(-(alpha_hat+beta_hat*dose)))
  }
  
  ## CI's of P(d): function to calcultae CI's, given confidence level
  ci_pr <- function(dose, level = 0.95) {
    # variance of (alpha_hat + beta_hat x d)
    Var <- Cov[1,1]+dose^2*Cov[2,2]+2*dose*Cov[1,2]

    # CI of (alpha + beta x d)
    R <- (alpha_hat+beta_hat*dose)-qnorm(c(1-(1-level)/2, (1-level)/2))*sqrt(Var)
    
    # CI of P(d)
    c(1/(1+exp(-R[1])), 1/(1+exp(-R[2])))
  }
  
  # CI of D50 = -alpha/beta
  ci_D50 <- function(level = 0.95) {
    D50_hat <- -alpha_hat/beta_hat
    
    # variance of D50_hat = -alpha_hat/beta_hat
    Var <- 1/beta_hat^2*(Cov[1,1]+D50_hat^2*Cov[2,2]+2*D50_hat*Cov[1,2])
                         
    # CI of D50
    D50_hat-qnorm(c(1-(1-level)/2, (1-level)/2))*sqrt(Var)
  }
    
  # CI of D100G = (g-alpha)/beta, where g = log(G/(1-G)), 0 < G < 1:
  # function to calculate CI of D100G, given confidence level
  ci_D100G <- function(G, level = 0.95) {
    g <- log(G/(1-G))
    D100G_hat <- (g-alpha_hat)/beta_hat
    
    # variance of D100G_hat = (g-alpha_hat)/beta_hat
    Var <- 1/beta_hat^2*(Cov[1,1]+D100G_hat^2*Cov[2,2]+2*D100G_hat*Cov[1,2])
    
    # CI of D100G
    list(
      D100G_hat = D100G_hat,
      ci = D100G_hat-qnorm(c(1-(1-level)/2, (1-level)/2))*sqrt(Var)
    )
    
  }

  list(
    alpha_hat = alpha_hat, 
    beta_hat = beta_hat,
    FIM = FIM,
    Cov = Cov,
    pr_hat = pr_hat,
    ci_pr = ci_pr,
    D50_hat = -alpha_hat/beta_hat,
    ci_D50 = ci_D50,
    ci_D100G = ci_D100G
    )
}



## D-optimal design for logistic model
doptimal_logistic <- function(alpha, beta){
  cval <- 1.5434
  #cval <- log(0.82/(1-0.82))
  list(
    d = (c(-cval, cval)-alpha)/beta,
    rho = rep(0.5, 2)
  )  
}

## A-optimal design for logistic model
aoptimal_logistic <- function(alpha, beta, type){
  if (type == "symmetric") {
    ofv <- function(x) {
     (alpha^2+beta^2)/(x^2*exp(x)/(1+exp(x))^2)+1/(exp(x)/(1+exp(x))^2)
    }
    opt <- nlminb(runif(1), ofv, control = list(step.min = 0.1, iter.max=1500))
    cval <- opt$par
    d <- (c(-cval, cval)-alpha)/beta
    rho <- c(0.5, 0.5)
  }
  else if (type == "non-symmetric") {
    ofv <- function(x){
      c1 <- x[1]
      C1 <- exp(-c1)/(1+exp(-c1))^2
      
      c2 <- x[2]
      C2 <- exp(-c2)/(1+exp(-c2))^2
      
      rho <- x[3]
      
      (rho*C1*(beta^2+(c1-alpha)^2) + (1-rho)*C2*(beta^2+(c2-alpha)^2))/
        ((rho*C1+(1-rho)*C2)*(rho*c1^2*C1+(1-rho)*c2^2*C2)-(rho*c1*C1+(1-rho)*c2*C2)^2)
    }
    opt <- nlminb(runif(3), ofv, control = list(step.min = 0.1, iter.max=1500))
    x <- opt$par
    d <- (c(x[1], x[2])-alpha)/beta
    rho <- c(x[3], 1-x[3])
  }
  else {
    stop("The input variable type may accept one of the values c('symmetric', 'non-symmetric')")
  }
  list(
    d = sort(d),
    rho = c(rho[d == min(d)], 1-rho[d == min(d)])
  )  
}

## A-optimal design for logistic model to estimate D100G and beta
aoptimal_logistic_G <- function(alpha, beta, G){
  
  g <- log(G/(1-G))
  ofv <- function(x) {
    1/beta^2*((beta^4+g^2)/(x^2*exp(x)/(1+exp(x))^2)+1/(exp(x)/(1+exp(x))^2))
  }
  opt <- nlminb(runif(1), ofv, control = list(step.min = 0.05, iter.max=1500))
  cval <- opt$par
  d <- (c(-cval, cval)-alpha)/beta
  rho <- c(0.5, 0.5)
    
  list(
    d = sort(d),
    rho = rho
  )  
}

# Directinal derivative function for D-optimality (sensitivity function)
phi <- function(dose, xi, theta) {
  # we do map over given doses
  map_dbl(dose, ~ {
    # Full Fisher Information Matrix
    M <- seq_along(xi$d) %>% map(~ {
      d <- xi$d[.]
      rho <- xi$rho[.]
      #w <- w^2/sum(xi$w)^2
      psi <- exp(-(theta[1]+theta[2]*d))/(1+exp(-(theta[1]+theta[2]*d)))^2
      rho*rbind(
        c(psi, psi*d),
        c(psi*d, psi*d^2)
      )
    }) %>% 
      reduce(`+`)
    
    # Fisher Information at dose given
    psi <- exp(-(theta[1]+theta[2]*.))/(1+exp(-(theta[1]+theta[2]*.)))^2
    Md <- rbind(
      c(psi, psi*.),
      c(psi*., psi*.^2)
    )
    
    # function value
    sum(diag(crossprod(solve(M), Md)))-2
  })
}


# Random Walk Rule
# - d -- dose range
# - N -- sample size
# - G -- target toxicity level (0 < G <= 0.5)
# - (alpha, beta) -- logistic model parameters
# - shiow_plot -- an option which tells the programm to plot RWR (if TRUE) or not (if FALSE)
rwr_design <- function(d, N, G, alpha, beta, show_plot = FALSE) {
  # bias coin probability
  b <- G/(1-G)
  
  # responses
  Y <- vector("integer", N)
  
  # dose assignments
  X <- vector("numeric", N)
  
  X[1] <- d[1]
  u <- runif(1)
  Y[1] <- if_else(u <= 1/(1+exp(-(alpha+beta*X[1]))), 1, 0)
  for (j in seq_len(N-1)) {
    u <- runif(1)
    if (X[j] == d[1]) { # lowest dose
      X[j+1] <- if_else(Y[j] == 1, d[1], if_else(u <= b, d[2], d[1]))
    }
    else if (X[j] == d[length(d)]){
      X[j+1] <- if_else(Y[j] == 1, d[length(d)-1], d[length(d)])
    }
    else {
      k <- seq_along(d)[d == X[j]]
      X[j+1] <- if_else(Y[j] == 1, d[k-1], if_else(u <= b, d[k+1], d[k]))
    }
    Y[j+1] <- if_else(u <= 1/(1+exp(-(alpha+beta*X[j+1]))), 1, 0)
  }
  return(data_frame(X, Y))
}



