library(tidyverse)
library(nleqslv)

## ===== auxiliary functions =====
# treatment assignment
assign_trt <- function(prob, trt_group){

  # cumulative probabilities
  cumulative_prob <- cumsum(c(0, prob))
  
  # draw a random number from the Uniform(0, 1)
  u <- runif(1)
  
  trt_group[seq_along(prob) %>% 
    map_lgl(~ {
      cumulative_prob[.] < u && u <= cumulative_prob[.+1]
    })]
}

# function to calculate allocation, given j (subject id) and vector of treatment assignments
get_allocation <- function(trt, ntrt) {
  map_dbl(seq_len(ntrt), ~ {
    sum(. == trt)
  })  
}

## ==========


# Complete randomization (CRD)
CRD <- function(w, nsbj) {
  # setup input
  rho <- w/sum(w)                                  # target allocation proportions
  trt <- vector("integer", nsbj)                   # vector of treatment assignments vs. allocation step
  imb <- vector("numeric", nsbj)                   # vector of imbalance vs. allocation step
  prob <- matrix(0, nrow = nsbj, ncol = length(w)) # matrix of treatment assignments
  
  for (j in seq_along(trt)) {
    trt[j] <- assign_trt(rho, seq_along(rho))
    N <- get_allocation(trt[seq_len(j)], length(w))
    imb[j] <- sqrt(sum((N-j*rho)^2))
    prob[j,] <- rho
  }
  
  # output
  list(
    op = tibble::tibble(
      procedure = "CRD",
      subject = seq_len(nsbj),
      imbalance = imb,
      FI = rep(0, nsbj)
    ),
    prob = as_tibble::tibble(prob) %>% 
      set_names(map_chr(seq_along(w), ~ paste0("pi[", ., "]")))%>% 
      add_column(subject = seq_len(nsbj), .before = 1) %>% 
      add_column(procedure = "CRD", .before = 1)
  )
}


# Random Allocation Rule (RAR)
RAR <- function(w, nsbj) {
  # setup input
  rho <- w/sum(w)                                  # target allocation proportions
  trt <- vector("integer", nsbj)                   # vector of treatment assignments vs. allocation step
  imb <- vector("numeric", nsbj)                   # vector of imbalance vs. allocation step
  fi <- vector("numeric", nsbj)                    # vector of FI vs. allocation step
  prob <- matrix(0, nrow = nsbj, ncol = length(w)) # matrix of treatment assignments
  
  for (j in seq_along(trt)) {
    if (j == 1) {prob[j,] <- w/sum(w)}
    else {prob[j,] <- (w-N)/(sum(w)-(j-1))}
    
    trt[j] <- assign_trt(prob[j,], seq_along(prob[j,]))
    N <- get_allocation(trt[seq_len(j)], length(w))
    imb[j] <- sqrt(sum((N-j*rho)^2))
    fi[j] <- sqrt(sum((prob[j,]-rho)^2))
  }
  
  # output
  list(
    op = tibble::tibble(
      procedure = "RAR",
      subject = seq_len(nsbj),
      imbalance = imb,
      FI = map_dbl(subject, ~ mean(fi[seq_len(.)]))
    ),
    prob = as_tibble::tibble(prob) %>% 
      set_names(map_chr(seq_along(w), ~ paste0("pi[", ., "]"))) %>% 
      add_column(subject = seq_len(nsbj), .before = 1) %>% 
      add_column(procedure = "RAR", .before = 1)
  )
}

# Truncated Multinomial Design
TMD <- function(w, nsbj) {
  # setup input
  rho <- w/sum(w)                                  # target allocation proportions
  trt <- vector("integer", nsbj)                   # vector of treatment assignments vs. allocation step
  imb <- vector("numeric", nsbj)                   # vector of imbalance vs. allocation step
  fi <- vector("numeric", nsbj)                    # vector of FI vs. allocation step
  prob <- matrix(0, nrow = nsbj, ncol = length(w)) # matrix of treatment assignments
  
  trt_group <- seq_along(rho)
  id <- rep(TRUE, length(trt_group))
  N <- rep(0, length(trt_group))
  
  for (j in seq_along(trt)) {
    prob[j,id] <- w[id]/sum(w[id])
    trt[j] <- assign_trt(prob[j,id], trt_group[id])
    N <- get_allocation(trt[seq_len(j)], length(w))
    # N[id] <- map_dbl(trt_group[id], ~ {
    #   sum(. == trt[seq_len(j)])
    # })
    imb[j] <- sqrt(sum((N-j*rho)^2))
    fi[j] <- sqrt(sum((prob[j,]-rho)^2))
    id <- !(N == w)
  }
  
  # output
  list(
    op = tibble::tibble(
      procedure = "TMD",
      subject = seq_len(nsbj),
      imbalance = imb,
      FI = map_dbl(subject, ~ mean(fi[seq_len(.)]))
    ),
    prob = as_tibble::tibble(prob) %>% 
      set_names(map_chr(seq_along(w), ~ paste0("pi[", ., "]"))) %>% 
      add_column(subject = seq_len(nsbj), .before = 1) %>% 
      add_column(procedure = "TMD", .before = 1)
  )
}

# Block Urn Design (BUD)
BUD <- function(w, nsbj, lambda) {
  # setup input
  rho <- w/sum(w)                                  # target allocation proportions
  trt <- vector("integer", nsbj)                   # vector of treatment assignments vs. allocation step
  imb <- vector("numeric", nsbj)                   # vector of imbalance vs. allocation step
  fi <- vector("numeric", nsbj)                    # vector of FI vs. allocation step
  prob <- matrix(0, nrow = nsbj, ncol = length(w)) # matrix of treatment assignments
  
  N <- rep(0, length(w))
  for (j in seq_along(trt)) {
    k <- floor(N/w)
    prob[j,] <- (w*(lambda+k)-N)/(sum(w)*(lambda+k)-(j-1))
    
    trt[j] <- assign_trt(prob[j,], seq_along(prob[j,]))
    N <- get_allocation(trt[seq_len(j)], length(w))
    imb[j] <- sqrt(sum((N-j*rho)^2))
    fi[j] <- sqrt(sum((prob[j,]-rho)^2))
    id <- !(N == w)
  }
  
  # output
  list(
    op = tibble::tibble(
      procedure = paste0("BUD (",lambda, ")"),
      subject = seq_len(nsbj),
      imbalance = imb,
      FI = map_dbl(subject, ~ mean(fi[seq_len(.)]))
    ),
    prob = as_tibble::tibble(prob) %>% 
      set_names(map_chr(seq_along(w), ~ paste0("pi[", ., "]")))%>% 
      add_column(subject = seq_len(nsbj), .before = 1) %>% 
      add_column(procedure = paste0("BUD (",lambda, ")"), .before = 1)
  )
    
}


# Mass Weighted Urn Design (MWUD)
MWUD <- function(w, nsbj, alpha) {
  # setup input
  rho <- w/sum(w)                                  # target allocation proportions
  trt <- vector("integer", nsbj)                   # vector of treatment assignments vs. allocation step
  imb <- vector("numeric", nsbj)                   # vector of imbalance vs. allocation step
  fi <- vector("numeric", nsbj)                    # vector of FI vs. allocation step
  prob <- matrix(0, nrow = nsbj, ncol = length(w)) # matrix of treatment assignments
  
  N <- rep(0, length(w))
  for (j in seq_along(trt)) {
    prob[j,] <- map_dbl(seq_along(rho), ~ max(alpha*rho[.] - N[.] + (j-1)*rho[.], 0))  
    prob[j,] <- prob[j,]/sum(prob[j,])
    
    trt[j] <- assign_trt(prob[j,], seq_along(prob[j,]))
    N <- get_allocation(trt[seq_len(j)], length(w))
    imb[j] <- sqrt(sum((N-j*rho)^2))
    fi[j] <- sqrt(sum((prob[j,]-rho)^2))
  }
  
  # output
  list(
    op = tibble::tibble(
      procedure = paste0("MWUD (",alpha, ")"),
      subject = seq_len(nsbj),
      imbalance = imb,
      FI = map_dbl(subject, ~ mean(fi[seq_len(.)]))
    ),
    prob = as_tibble::tibble(prob) %>% 
      set_names(map_chr(seq_along(w), ~ paste0("pi[", ., "]")))%>% 
      add_column(subject = seq_len(nsbj), .before = 1) %>% 
      add_column(procedure = paste0("MWUD (",alpha, ")"), .before = 1)
  )
   
  
}


# Drop-the-Loser Rule
DL <- function(w, nsbj, a){
  # setup input
  rho <- w/sum(w)                                  # target allocation proportions
  trt <- vector("integer", nsbj)                   # vector of treatment assignments vs. allocation step
  imb <- vector("numeric", nsbj)                   # vector of imbalance vs. allocation step
  fi <- vector("numeric", nsbj)                    # vector of FI vs. allocation step
  prob <- matrix(0, nrow = nsbj, ncol = length(w)) # matrix of treatment assignments
  
  urn <- c(1, w)
  for ( j in seq_along(trt)) {
    flag <- TRUE
    while(flag) {
      ball <- sample(c(0, seq_along(w)), 1, prob = urn/sum(urn)) # sample a ball from the urn
      if (ball == 0) { # "immigration ball is drawn"
        urn[-1] <- urn[-1]+a*w # additional balls are put into the urn
      }
      else { # the corresponding treatment (ball number) is assigned
        trt[j] <- ball
        prob[j,] <- urn[-1]/sum(urn[-1])
        urn[ball+1] <- urn[ball+1]-1
        flag <- FALSE
      }
    }
    
    N <- get_allocation(trt[seq_len(j)], length(w))
    imb[j] <- sqrt(sum((N-j*rho)^2))
    fi[j] <- sqrt(sum((prob[j,]-rho)^2))
  }
  
  # output
  list(
    op = tibble::tibble(
      procedure = paste0("DL (",a, ")"),
      subject = seq_len(nsbj),
      imbalance = imb,
      FI = map_dbl(subject, ~ mean(fi[seq_len(.)]))
    ),
    prob = as_tibble::tibble(prob) %>% 
      set_names(map_chr(seq_along(w), ~ paste0("pi[", ., "]")))%>% 
      add_column(subject = seq_len(nsbj), .before = 1) %>% 
      add_column(procedure = paste0("DL (",a, ")"), .before = 1)
  )
    
}


# Doubly-Adaptive Biased Coin Design
DBCD <- function(w, nsbj, gm){
  # setup input
  rho <- w/sum(w)                                  # target allocation proportions
  trt <- vector("integer", nsbj)                   # vector of treatment assignments vs. allocation step
  imb <- vector("numeric", nsbj)                   # vector of imbalance vs. allocation step
  fi <- vector("numeric", nsbj)                    # vector of FI vs. allocation step
  prob <- matrix(0, nrow = nsbj, ncol = length(w)) # matrix of treatment assignments
  
  N <- rep(0, length(w))
  
  for (j in seq_along(trt)) {
    if (all(N > 0)) {
      prob[j,] <- rho*(rho/(N/j))^gm/sum(rho*(rho/(N/j))^gm)
    }
    else {
      prob[j,] <- rho
    }
    
    trt[j] <- assign_trt(prob[j,], seq_along(prob[j,]))
    N <- get_allocation(trt[seq_len(j)], length(w))
    imb[j] <- sqrt(sum((N-j*rho)^2))
    fi[j] <- sqrt(sum((prob[j,]-rho)^2))
  } 
  # output
  list(
    op = tibble::tibble(
      procedure = paste0("DBCD (",gm, ")"),
      subject = seq_len(nsbj),
      imbalance = imb,
      FI = map_dbl(subject, ~ mean(fi[seq_len(.)]))
    ),
    prob = as_tibble::tibble(prob) %>% 
      set_names(map_chr(seq_along(w), ~ paste0("pi[", ., "]")))%>% 
      add_column(subject = seq_len(nsbj), .before = 1) %>% 
      add_column(procedure = paste0("DBCD (",gm, ")"), .before = 1)
  )
  
  
}


# Maximum Entropy Constrained Balance Randomization
MaxEnt <- function(w, nsbj, eta) {
  if (eta < 0 || eta > 1) {
    stop(paste("The value of randomizaton procedure parameter must belong to [0, 1]: current value equals to", eta)) 
  }
  
  # setup input
  rho <- w/sum(w)                                  # target allocation proportions
  trt <- vector("integer", nsbj)                   # vector of treatment assignments vs. allocation step
  imb <- vector("numeric", nsbj)                   # vector of imbalance vs. allocation step
  fi <- vector("numeric", nsbj)                    # vector of FI vs. allocation step
  prob <- matrix(0, nrow = nsbj, ncol = length(w)) # matrix of treatment assignments
  
  N <- rep(0, length(w))
  
  for(j in seq_along(trt)) {
    # calculate "hypothetical" imbalance
    B <- map_dbl(seq_along(w), ~ {
      N1 <- N
      N1[.] <- N1[.]+1
      max(abs(N1/j-rho))
      #N1[.] <- N1[.]-1
    })
    if (eta != 1) {
      if (var(B) <= 1e-16) {prob[j,] <- rho}
      else {
        # we have to find a zero of a one-variable (mu) function.

        # function to find zero of
        obj <- function(mu, eta, B, rho) {
         min(B)*eta + (1-eta)*sum(rho*B) - sum(B*rho*exp(-mu*B))/sum(rho*exp(-mu*B))
        }
        
        # function's derivative
        dobj <- function(mu, eta, B, rho) {
          sum(B*B*rho*exp(-mu*(B-min(B))))/sum(rho*exp(-mu*(B-min(B)))) -
            (sum(B*rho*exp(-mu*(B-min(B))))/sum(rho*exp(-mu*(B-min(B)))))^2
        }
        
        mu <- nleqslv(5, function(x) obj(x, eta, B, rho))$x
        prob[j,] <- rho*exp(-mu*B)/sum(rho*exp(-mu*B))
      }
    }
    else {
      prob[j,B == min(B)] <- rho[B == min(B)]/sum(rho[B == min(B)])
    }
    trt[j] <- assign_trt(prob[j,], seq_along(prob[j,]))
    N <- get_allocation(trt[seq_len(j)], length(w))
    imb[j] <- sqrt(sum((N-j*rho)^2))
    fi[j] <- sqrt(sum((prob[j,]-rho)^2))
  }
  
  # output
  list(
    op = tibble::tibble(
      procedure = paste0("MaxEnt (",eta, ")"),
      subject = seq_len(nsbj),
      imbalance = imb,
      FI = map_dbl(subject, ~ mean(fi[seq_len(.)]))
    ),
    prob = as_tibble::tibble(prob) %>% 
      set_names(map_chr(seq_along(w), ~ paste0("pi[", ., "]"))) %>% 
      add_column(subject = seq_len(nsbj), .before = 1) %>% 
      add_column(procedure = paste0("MaxEnt (",eta, ")"), .before = 1)
  )
   
}

## ===== simulations =====
w <- c(182, 55, 55)
nsbj <- 292
nsim <- 1000

crd_ls <- seq_len(nsim) %>% 
  map(~ CRD(w, nsbj)) 

rar_ls <- seq_len(nsim) %>% 
  map(~ RAR(w, nsbj)) 

tmd_ls <- seq_len(nsim) %>% 
  map(~ TMD(w, nsbj)) 

bud_ls2 <- seq_len(nsim) %>% 
  map(~ BUD(w, nsbj, lambda = 2)) 
  
bud_ls3 <- seq_len(nsim) %>% 
  map(~ BUD(w, nsbj, lambda = 3)) 

mwud_ls2 <- seq_len(nsim) %>% 
  map(~ MWUD(w, nsbj, alpha = 2)) 

mwud_ls4 <- seq_len(nsim) %>% 
  map(~ MWUD(w, nsbj, alpha = 4)) 

mwud_ls6 <- seq_len(nsim) %>% 
  map(~ MWUD(w, nsbj, alpha = 6)) 

mwud_ls8 <- seq_len(nsim) %>% 
  map(~ MWUD(w, nsbj, alpha = 8)) 

dl_ls2 <- seq_len(nsim) %>% 
  map(~ DL(w, nsbj, a = 2)) 

dl_ls4 <- seq_len(nsim) %>% 
  map(~ DL(w, nsbj, a = 4)) 

dl_ls6 <- seq_len(nsim) %>% 
  map(~ DL(w, nsbj, a = 6)) 

dl_ls8 <- seq_len(nsim) %>% 
  map(~ DL(w, nsbj, a = 8)) 

dbcd_ls1 <- seq_len(nsim) %>% 
  map(~ DBCD(w, nsbj, gm = 1)) 

dbcd_ls2 <- seq_len(nsim) %>% 
  map(~ DBCD(w, nsbj, gm = 2)) 

dbcd_ls5 <- seq_len(nsim) %>% 
  map(~ DBCD(w, nsbj, gm = 5)) 

dbcd_ls10 <- seq_len(nsim) %>% 
  map(~ DBCD(w, nsbj, gm = 10))

max_ent_ls005 <- seq_len(nsim) %>% 
  map(~ MaxEnt(w, nsbj, eta = 0.05)) 

max_ent_ls025 <- seq_len(nsim) %>% 
  map(~ MaxEnt(w, nsbj, eta = 0.25)) 

max_ent_ls05 <- seq_len(nsim) %>% 
  map(~ MaxEnt(w, nsbj, eta = 0.5)) 

max_ent_ls1 <- seq_len(nsim) %>% 
  map(~ MaxEnt(w, nsbj, eta = 1)) 

op <- list(crd_ls, rar_ls, tmd_ls, bud_ls2, bud_ls3, mwud_ls2, mwud_ls4, mwud_ls6, mwud_ls6, mwud_ls8,
           dl_ls2, dl_ls4, dl_ls6, dl_ls8, dbcd_ls1, dbcd_ls2, dbcd_ls5, dbcd_ls10, max_ent_ls005, 
           max_ent_ls025, max_ent_ls05, max_ent_ls1) %>% map(~ {
  ls <- .
  ls %>% map(~ .$op) %>% 
    bind_rows() %>% 
    group_by(procedure, subject) %>% 
    summarise(
      MPM = mean(imbalance),
      AFI = mean(FI)
    ) %>% 
    mutate(CMPM = map_dbl(seq_along(MPM), ~ mean(MPM[seq_len(.)])))
  }) %>% 
  bind_rows() %>% 
  ungroup() %>% 
  mutate(
    design = map_chr(procedure, ~ str_split(., " ")[[1]][1])
  ) 

save(op, file = "./data/lecture04-op.Rda")

load(file = "./data/lecture04-op.Rda")
# plot imbalance vs. allocation step
op %>% 
  ggplot(aes(x = subject, y = CMPM, color = procedure, shape = procedure))+
  geom_point(size = 1.25)+
  geom_line(size = 0.5)+
  xlab("allocation step")+
  ylab("cumulative imbalance")+
  scale_shape_manual(values = c(
    "BUD (2)" = 0, "BUD (3)" = 1,
    "CRD" = 3, "RAR" = 4, "TMD" = 8,
    "DBCD (1)" = 0, "DBCD (2)" = 1, "DBCD (5)" = 2, "DBCD (10)" = 3,
    "MWUD (2)" = 0, "MWUD (4)" = 1, "MWUD (6)" = 2, "MWUD (8)" = 3, 
    "DL (2)" = 0, "DL (4)" = 1, "DL (6)" = 2, "DL (8)" = 3, 
    "MaxEnt (0.05)" = 0, "MaxEnt (0.25)" = 1, "MaxEnt (0.5)" = 2, "MaxEnt (1)" = 3
  )) +
  scale_color_manual(values = c(
    "BUD (2)" = "red", "BUD (3)" = "darkblue",
    "CRD" = "red", "RAR" = "purple", "TMD" = "violet",
    "DBCD (1)" = "red", "DBCD (2)" = "darkblue", "DBCD (5)" = "darkgreen", "DBCD (10)" = "darkorange",
    "MWUD (2)" = "red", "MWUD (4)" = "darkblue", "MWUD (6)" = "darkgreen", "MWUD (8)" = "darkorange", 
    "DL (2)" = "red", "DL (4)" = "darkblue", "DL (6)" = "darkgreen", "DL (8)" = "darkorange", 
    "MaxEnt (0.05)" = "red", "MaxEnt (0.25)" = "darkblue", "MaxEnt (0.5)" = "darkgreen", "MaxEnt (1)" = "darkorange"
  )) +
  facet_wrap(~ design, ncol = 1, scales = "free_y")+
  theme(
    axis.text = element_text(family = "Helvetica", face = "bold", size = 14),
    axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
    strip.text = element_text(family = "Helvetica", face = "bold", size = 18),
    legend.title = element_blank(),
    legend.text = element_text(family = "Helvetica", face = "bold", size = 12),
    legend.position = "right"
  )
ggsave("./figures/lecture04-fig04-imbalance.png", width = 16, height = 9, units = "in")

# plot average FI vs. allocation step
op %>% 
  ggplot(aes(x = subject, y = AFI, color = procedure, shape = procedure))+
  geom_point(size = 1.25)+
  geom_line(size = 0.5)+
  xlab("allocation step")+
  ylab("average forcing index")+
  scale_shape_manual(values = c(
    "BUD (2)" = 0, "BUD (3)" = 1,
    "CRD" = 3, "RAR" = 4, "TMD" = 8,
    "DBCD (1)" = 0, "DBCD (2)" = 1, "DBCD (5)" = 2, "DBCD (10)" = 3,
    "MWUD (2)" = 0, "MWUD (4)" = 1, "MWUD (6)" = 2, "MWUD (8)" = 3, 
    "DL (2)" = 0, "DL (4)" = 1, "DL (6)" = 2, "DL (8)" = 3, 
    "MaxEnt (0.05)" = 0, "MaxEnt (0.25)" = 1, "MaxEnt (0.5)" = 2, "MaxEnt (1)" = 3
  )) +
  scale_color_manual(values = c(
    "BUD (2)" = "red", "BUD (3)" = "darkblue",
    "CRD" = "red", "RAR" = "purple", "TMD" = "violet",
    "DBCD (1)" = "red", "DBCD (2)" = "darkblue", "DBCD (5)" = "darkgreen", "DBCD (10)" = "darkorange",
    "MWUD (2)" = "red", "MWUD (4)" = "darkblue", "MWUD (6)" = "darkgreen", "MWUD (8)" = "darkorange", 
    "DL (2)" = "red", "DL (4)" = "darkblue", "DL (6)" = "darkgreen", "DL (8)" = "darkorange", 
    "MaxEnt (0.05)" = "red", "MaxEnt (0.25)" = "darkblue", "MaxEnt (0.5)" = "darkgreen", "MaxEnt (1)" = "darkorange"
  )) +
  facet_wrap(~ design, ncol = 1, scales = "free_y")+
  theme(
    axis.text = element_text(family = "Helvetica", face = "bold", size = 14),
    axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
    strip.text = element_text(family = "Helvetica", face = "bold", size = 18),
    legend.title = element_blank(),
    legend.text = element_text(family = "Helvetica", face = "bold", size = 12),
    legend.position = "right"
  )
ggsave("./figures/lecture04-fig05-afi.png", width = 16, height = 9, units = "in")


# average FI vs. average cumulative imbalance
op %>% 
  filter(subject %in% c(73, 146, 219, 292)) %>% 
  ggplot()+
  geom_point(aes(x = CMPM, y = AFI, fill = procedure, color = procedure, shape = design), size = 3.25)+
  #geom_line(size = 0.5)+
  xlab("allocation step")+
  ylab("average forcing index")+
  scale_shape_manual(values = c(
    "BUD" = 21,
    "CRD" = 22, "RAR" = 23, "TMD" = 24,
    "DBCD" = 25,
    "MWUD" = 3, 
    "DL" = 4, 
    "MaxEnt" = 8
  )) +
  scale_color_manual(values = c(
    "BUD (2)" = "red", "BUD (3)" = "darkblue",
    "CRD" = "red", "RAR" = "purple", "TMD" = "violet",
    "DBCD (1)" = "red", "DBCD (2)" = "darkblue", "DBCD (5)" = "darkgreen", "DBCD (10)" = "darkorange",
    "MWUD (2)" = "red", "MWUD (4)" = "darkblue", "MWUD (6)" = "darkgreen", "MWUD (8)" = "darkorange", 
    "DL (2)" = "red", "DL (4)" = "darkblue", "DL (6)" = "darkgreen", "DL (8)" = "darkorange", 
    "MaxEnt (0.05)" = "red", "MaxEnt (0.25)" = "darkblue", "MaxEnt (0.5)" = "darkgreen", "MaxEnt (1)" = "darkorange"
  )) +
  scale_fill_manual(values = c(
    "BUD (2)" = "red", "BUD (3)" = "darkblue",
    "CRD" = "red", "RAR" = "purple", "TMD" = "violet",
    "DBCD (1)" = "red", "DBCD (2)" = "darkblue", "DBCD (5)" = "darkgreen", "DBCD (10)" = "darkorange",
    "MWUD (2)" = "red", "MWUD (4)" = "darkblue", "MWUD (6)" = "darkgreen", "MWUD (8)" = "darkorange", 
    "DL (2)" = "red", "DL (4)" = "darkblue", "DL (6)" = "darkgreen", "DL (8)" = "darkorange", 
    "MaxEnt (0.05)" = "red", "MaxEnt (0.25)" = "darkblue", "MaxEnt (0.5)" = "darkgreen", "MaxEnt (1)" = "darkorange"
  )) +
  facet_wrap(~ subject, ncol = 2, scales = "free_y")+
  theme(
    axis.text = element_text(family = "Helvetica", face = "bold", size = 14),
    axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
    strip.text = element_text(family = "Helvetica", face = "bold", size = 18),
    legend.title = element_blank(),
    legend.text = element_text(family = "Helvetica", face = "bold", size = 12),
    legend.position = "right"
  )
ggsave("./figures/lecture04-fig06-afi-vs-imb.png", width = 16, height = 9, units = "in")


# rank - distance to the origin

rank_ls <- c(73, 146, 219, 292) %>% map(~ {
  nsbj <- .
  op %>% 
    filter(subject == nsbj) %>% 
    mutate(
      dist = sqrt(AFI^2+CMPM^2),
      out = paste(procedure, "/", round(dist, 3))
    ) %>% 
    arrange(dist) %>% 
    select(out) %>% 
    add_column(rank = seq_len(nrow(.)), .before = 1) %>%
    set_names(c("rank", paste0("procedure / distance to the origin (n = ", nsbj, ")")))
})

rank_df <- inner_join(
  inner_join(
    rank_ls[[1]],
    rank_ls[[2]]
  ),
  inner_join(
    rank_ls[[3]],
    rank_ls[[4]]
  )
)
save(rank_df, file="./data/rank_df.Rda")


# ARP property

prob <- list(crd_ls, rar_ls, tmd_ls, bud_ls2, bud_ls3, mwud_ls2, mwud_ls4, mwud_ls6, mwud_ls6, mwud_ls8,
           dl_ls2, dl_ls4, dl_ls6, dl_ls8, dbcd_ls1, dbcd_ls2, dbcd_ls5, dbcd_ls10, max_ent_ls005, 
           max_ent_ls025, max_ent_ls05, max_ent_ls1) %>% map(~ {
             ls <- .
             ls %>% map(~ .$prob) %>% 
               bind_rows()
           }) %>% 
  bind_rows() %>% 
  gather(key, value, -procedure, -subject) %>% 
  group_by(procedure, key, subject) %>% 
  summarise(
    value = mean(value)
  )

save(prob, file = "./data/lecture04-prob.Rda")

prob %>% 
  ggplot(aes(x = subject, y = value, color = key))+
  geom_point(size = 1)+
  geom_line(size = 0.5)+
  xlab("allocation step")+
  ylab("probability of treatment assignmet")+
  facet_wrap(~ procedure, ncol = 7, scales = "free_y")+
  theme(
    axis.text = element_text(family = "Helvetica", face = "bold", size = 14),
    axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
    strip.text = element_text(family = "Helvetica", face = "bold", size = 18),
    legend.title = element_blank(),
    legend.text = element_text(family = "Helvetica", face = "bold", size = 12),
    legend.position = "right"
  )
ggsave("./figures/lecture04-fig07-arp.png", width = 16, height = 9, units = "in")
