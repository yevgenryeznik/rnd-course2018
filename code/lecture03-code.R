library(tidyverse)

# Complete Randomization
CRD <- function(nsbj) {
  trt <- vector("integer", nsbj)  # vector of treatment assignments vs. allocation step
  fi <- vector("numeric", nsbj)   # vector of partial FI vs. allocation step
  imb <- vector("numeric", nsbj)  # vector of imbalance vs. allocation step
  loss <- vector("numeric", nsbj) # vector of loss vs. allocation step
  
  prob <- 0.5
  for (j in seq_len(nsbj)) {
    trt[j] <- as.numeric(runif(1) <= prob)
    fi[j] <- abs(prob-0.5)/0.5
    imb[j] <- 2*sum(trt[seq_len(j)] == 1) - j
    loss[j] <- imb[j]^2/j
  }
  
  # output
  tibble::tibble(
    procedure = "CRD",
    subject = seq_len(nsbj),
    FI = map_dbl(subject, ~ mean(fi[seq_len(.)])),
    imbalance = imb, 
    Loss = loss
  )
}


# Efron's BCD
Efron_BCD <- function(nsbj, p) {
  trt <- vector("integer", nsbj)  # vector of treatment assignments vs. allocation step
  fi <- vector("numeric", nsbj)   # vector of partial FI vs. allocation step
  imb <- vector("numeric", nsbj)  # vector of imbalance vs. allocation step
  loss <- vector("numeric", nsbj) # vector of loss vs. allocation step
  
  for (j in seq_len(nsbj)) {
    if (j == 1) {prob <- 0.5}
    else {
      prob <- 1/2+sign(imb[j-1])*(1/2-p)
    }
    
    trt[j] <- as.numeric(runif(1) <= prob)
    fi[j] <- abs(prob-0.5)/0.5
    imb[j] <- 2*sum(trt[seq_len(j)] == 1) - j
    loss[j] <- imb[j]^2/j
  }
  
  # output
  tibble::tibble(
    procedure = paste0("Efron's BCD (p = ", round(p, 3), ")"),
    subject = seq_len(nsbj),
    FI = map_dbl(subject, ~ mean(fi[seq_len(.)])),
    imbalance = imb, 
    Loss = loss
  )
}

# Adjustable BCD
ABCD <- function(nsbj, a) {
  trt <- vector("integer", nsbj)  # vector of treatment assignments vs. allocation step
  fi <- vector("numeric", nsbj)   # vector of partial FI vs. allocation step
  imb <- vector("numeric", nsbj)  # vector of imbalance vs. allocation step
  loss <- vector("numeric", nsbj) # vector of loss vs. allocation step
  
  for (j in seq_len(nsbj)) {
    if (j == 1) {prob <- 0.5}
    else {
      if (imb[j-1] <= -1) {prob <- abs(imb[j-1])^a/(abs(imb[j-1])^a+1)}
      else if (imb[j-1] >= 1) {prob <- 1/(abs(imb[j-1])^a+1)}
      else {prob <- 0.5}
    }
    
    trt[j] <- as.numeric(runif(1) <= prob)
    fi[j] <- abs(prob-0.5)/0.5
    imb[j] <- 2*sum(trt[seq_len(j)] == 1) - j
    loss[j] <- imb[j]^2/j
  }
  
  # output
  tibble::tibble(
    procedure = paste0("ABCD (a = ", a, ")"),
    subject = seq_len(nsbj),
    FI = map_dbl(subject, ~ mean(fi[seq_len(.)])),
    imbalance = imb, 
    Loss = loss
  )
}


# Ehrenfest Urn Design
Ehrenfest_UD <- function(nsbj, w) {
  trt <- vector("integer", nsbj)  # vector of treatment assignments vs. allocation step
  fi <- vector("numeric", nsbj)   # vector of partial FI vs. allocation step
  imb <- vector("numeric", nsbj)  # vector of imbalance vs. allocation step
  loss <- vector("numeric", nsbj) # vector of loss vs. allocation step
  
  for (j in seq_len(nsbj)) {
    if (j == 1) {prob <- 0.5}
    else {
      if (imb[j-1] <= -w/2) {prob <- 1}
      else if (imb[j-1] >= w/2) {prob <- 0}
      else {prob <- 0.5-imb[j-1]/w}
    }
    
    trt[j] <- as.numeric(runif(1) <= prob)
    fi[j] <- abs(prob-0.5)/0.5
    imb[j] <- 2*sum(trt[seq_len(j)] == 1) - j
    loss[j] <- imb[j]^2/j
  }
  
  # output
  tibble::tibble(
    procedure = paste0("ABCD (a = ", a, ")"),
    subject = seq_len(nsbj),
    FI = map_dbl(subject, ~ mean(fi[seq_len(.)])),
    imbalance = imb, 
    Loss = loss
  )
}

# Big Stick
Big_Stick <- function(nsbj, b) {
  trt <- vector("integer", nsbj)  # vector of treatment assignments vs. allocation step
  fi <- vector("numeric", nsbj)   # vector of partial FI vs. allocation step
  imb <- vector("numeric", nsbj)  # vector of imbalance vs. allocation step
  loss <- vector("numeric", nsbj) # vector of loss vs. allocation step
  
  for (j in seq_len(nsbj)) {
    if (j == 1) {prob <- 0.5}
    else {
      if (abs(imb[j-1]) < b) {prob <- 0.5}
      else if (imb[j-1] == b) {prob <- 0}
      else {prob <- 1}
    }
    
    trt[j] <- as.numeric(runif(1) <= prob)
    fi[j] <- abs(prob-0.5)/0.5
    imb[j] <- 2*sum(trt[seq_len(j)] == 1) - j
    loss[j] <- imb[j]^2/j
  }
  
  # output
  tibble::tibble(
    procedure = paste0("ABCD (a = ", a, ")"),
    subject = seq_len(nsbj),
    FI = map_dbl(subject, ~ mean(fi[seq_len(.)])),
    imbalance = imb, 
    Loss = loss
  )
}


# Generalized BCD
GBCD <- function(nsbj, rho) {
  trt <- vector("integer", nsbj)  # vector of treatment assignments vs. allocation step
  fi <- vector("numeric", nsbj)   # vector of partial FI vs. allocation step
  imb <- vector("numeric", nsbj)  # vector of imbalance vs. allocation step
  loss <- vector("numeric", nsbj) # vector of loss vs. allocation step
  
  for (j in seq_len(nsbj)) {
    if (j == 1) {prob <- 0.5}
    else {
      prob <- (1-imb[j-1]/(j-1))^rho/((1-imb[j-1]/(j-1))^rho + (1+imb[j-1]/(j-1))^rho)
    }
    
    trt[j] <- as.numeric(runif(1) <= prob)
    fi[j] <- abs(prob-0.5)/0.5
    imb[j] <- 2*sum(trt[seq_len(j)] == 1) - j
    loss[j] <- imb[j]^2/j
  }
  
  # output
  tibble::tibble(
    procedure = paste0("GBCD (rho = ", rho, ")"),
    subject = seq_len(nsbj),
    FI = map_dbl(subject, ~ mean(fi[seq_len(.)])),
    imbalance = imb, 
    Loss = loss
  )
}


# Wei's UD
Wei_UD <- function(nsbj, alpha, beta){
  trt <- vector("integer", nsbj)  # vector of treatment assignments vs. allocation step
  fi <- vector("numeric", nsbj)   # vector of partial FI vs. allocation step
  imb <- vector("numeric", nsbj)  # vector of imbalance vs. allocation step
  loss <- vector("numeric", nsbj) # vector of loss vs. allocation step
  
  for (j in seq_len(nsbj)) {
    if (j == 1) {prob <- 0.5}
    else {
      #prob <- 1/2*(1-beta*imb[j-1]/(2*alpha+beta*(j-1)))
      prob <- (1-beta*sum(trt[seq_len(j)] == 1)/(2*alpha+beta*(j-1)))
    }

    trt[j] <- as.numeric(runif(1) <= prob)
    fi[j] <- abs(prob-0.5)/0.5
    imb[j] <- 2*sum(trt[seq_len(j)] == 1) - j
    loss[j] <- imb[j]^2/j
  }
  
  # output
  tibble::tibble(
    procedure = paste0("Wei's UD (alpha = ", alpha, ", beta = ", beta, ")"),
    subject = seq_len(nsbj),
    FI = map_dbl(subject, ~ mean(fi[seq_len(.)])),
    imbalance = imb, 
    Loss = loss
  )
} 


## ===== Simulation study =====
nsim <- 1000
nsbj <- 200

# CRD
crd <- map(seq_len(nsim), ~ CRD(nsbj)) %>% 
  bind_rows() %>% 
  mutate(procedure = paste("1.", procedure))

# Efron's BCD (p = 2/3)
efron_bcd <- map(seq_len(nsim), ~ Efron_BCD(nsbj, 2/3)) %>% 
  bind_rows() %>% 
  mutate(procedure = "2. Efron's BCD (p = 2/3)")

# Adjustable BCD with a = 1
abcd <- map(seq_len(nsim), ~ ABCD(nsbj, 1)) %>% 
  bind_rows() %>% 
  mutate(procedure = paste("3.", procedure))

# Wei's UD with alpha = 0, beta = 1
wei_ud <- map(seq_len(nsim), ~ Wei_UD(nsbj, 0, 1)) %>% 
  bind_rows() %>% 
  mutate(procedure = paste("4.", procedure))

# Generalized BCD 
gbcd <- bind_rows(
  # with rho = 2
  map(seq_len(nsim), ~ GBCD(nsbj, 2)) %>% 
    bind_rows() %>% 
    mutate(procedure = paste("5.", procedure)),
  
  # with rho = 5
  map(seq_len(nsim), ~ GBCD(nsbj, 5)) %>% 
    bind_rows() %>% 
    mutate(procedure = paste("6.", procedure)),
  
  # with rho = 20
  map(seq_len(nsim), ~ GBCD(nsbj, 20)) %>% 
    bind_rows() %>% 
    mutate(procedure = paste("7.", procedure))
)

simulation_study <- bind_rows(
  crd, 
  efron_bcd,
  abcd,
  wei_ud,
  gbcd
) %>% 
  group_by(procedure, subject) %>% 
  summarise(
    `2. Average FI` = mean(FI),
    `1. Average Loss` = mean(Loss)
  )

save(simulation_study, file = "./data/lecture03-sim-study.Rda")
load(file = "./data/lecture03-sim-study.Rda")
