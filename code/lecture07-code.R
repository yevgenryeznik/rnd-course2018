library(tidyverse)
library(animation)


## ===== an R code to create a gif =====
# 1) ImageMagick has to be installed first: https://www.imagemagick.org/script/index.php
# 2) animation package has to be installed from GitHib: run the command 
#       devtools::install_github('yihui/animation')

Pr <- 0.7
Y <- rbinom(50, 1, Pr)
# fucntion to calculate a posterior
posterior <- function(p, X, n, a, b) {
  dbeta(p, X+a, n-X+b)
}

a <- 1
b <- 1

p <- seq(0, 1, by = 0.01)
f <- posterior(p, 0, 0, a, b)

low <- qbeta(0.025, 0, 0)
upp <- qbeta(0.975, 1, 0)
  
saveGIF({
  ani.options(nmax = 10)
  plot(p, f, 
       main = substitute(paste("n = ", n, ": ", p %~% Beta (a, b), ", 95% CI for p: (", low, ",", upp, ")"), list(n = 0, a = a, b = b, low = low, upp = upp)), 
       xlab = "Probability of success ( p )",
       ylab = "p. d. f. of p ( f(p) )",
       #ylim = c(0, 30),
       type="l")
  lines(c(Pr, Pr), c(0, 30), col = "red", lwd = 1.2, lty = "dashed")
  grid()
  
  for(j in seq_along(Y)) {
    ani.pause()
    
    n <- j
    X <- sum(Y[seq_len(n)])
    f <- posterior(p, X, n, a, b)
    
    a <- X + a
    b <- n-X+b
    
    low <- round(qbeta(0.025, X, n-X), 2)
    upp <- round(qbeta(0.975, X+1, n-X), 2)
    
    plot(p, f, 
         main = substitute(paste("n = ", n, ": ", p %~% Beta (a, b), ", 95% CI for p: (", low, ",", upp, ")"), list(n = n, a = a, b = b, low = low, upp = upp)), 
         xlab = "Probability of success ( p )",
         ylab = "p. d. f. of p ( f(p) )",
         #ylim = c(0, 30),
         type="l")
    lines(c(Pr, Pr), c(0, 30), col = "red", lwd = 1.2, lty = "dashed")
    grid()
  } 
}, interval = 0.3, movie.name = 'bayesian_demo.gif', ani.width = 600, ani.height = 600)


## Cumulative distribution function (CDF) and probability density function (PDF) of 
## Delta = pE-pC

# pE ~ Beta(a1, b1)
# pC ~ Beta(a2, b2)
# Pr(Delta = pE-pC <= x) = delta_cdf(x, a1, b1, a2, b2)

# CDF
delta_cdf <- function(x, a1, b1, a2, b2) {
  seq_along(x) %>% map_dbl(~ {
    if (x[.] <= -1) 0
    else if (-1 < x[.] & x[.] <= 0) integrate(function(t) pbeta(x[.]+t, a2, b2)*dbeta(t, a1, b1), -x[.], 1)$value
    else if (0 < x[.] & x[.] < 1) integrate(function(t) pbeta(x[.]+t, a2, b2)*dbeta(t, a1, b1), 0, 1-x[.])$value + 
      integrate(function(t) dbeta(t, a1, b1), 1-x[.], 1)$value
    else 1
  })
}

# PDF
delta_pdf <- function(x, a1, b1, a2, b2) {
  seq_along(x) %>% map_dbl(~ {
    if (x[.] <= -1) 0
    else if (-1 < x[.] & x[.] <= 0) integrate(function(t) dbeta(x[.] + t, a2, b2)*dbeta(t, a1, b1), -x[.], 1)$value + pbeta(0, a2, b2)*dbeta(-x[.], a1, b1)
    else if (0 < x[.] & x[.] <= 1) integrate(function(t) dbeta(x[.] + t, a2, b2)*dbeta(t, a1, b1), 0, 1 - x[.])$value - pbeta(1, a2, b2)*dbeta(1-x[.], a1, b1)+
      dbeta(1-x[.], a1, b1)
    else 0
  })
} 


nsbj <- 10
trt <- sample(c(0, 1), nsbj, TRUE, c(0.5, 0.5))
sum(trt==1)

a1 <- 1/3
b1 <- 1/3

a2 <- 1/3 
b2 <- 1/3

Y <- vector("integer", nsbj)
Y[trt == 1] <- rbinom(sum(trt==1), 1, 0.3)
Y[trt == 0] <- rbinom(sum(trt==0), 1, 0.7) 
  
Delta <- seq(-0.99, 0.99, by = 0.05)
saveGIF({
  ani.options(nmax = 50)
  
plot(Delta, delta_pdf(Delta, a1, b1, a2, b2), 
     main = substitute(paste(n[1], "=", n1, ", ", n[2], "=", n2, ": ", p[1] %~% Beta (a1, b1), ", ", p[1] %~% Beta (a2, b2)), list(n1 = 0, n2 = 0, a1 = round(a1, 3), b1 = round(b1, 3), a2 = round(a2, 3), b2 = round(b2, 3))), 
     xlab = "Probability of success ( p )",
     ylab = "p. d. f. of p ( f(p) )",
     ylim = c(0, 6.5),
     type="l")
grid()

for(j in seq_along(Y)) {
  ani.pause()
  
  n <- j
  n1 <- sum(trt[seq_len(n)])
  n2 <- sum(1-trt[seq_len(n)])
  
  X1 <- sum(Y[seq_len(n)]*trt[seq_len(n)])
  X2 <- sum(Y[seq_len(n)]*(1-trt[seq_len(n)]))

  a1 <- X1 + a1
  b1 <- n-X1+b1
  
  a2 <- X2 + a1
  b2 <- n-X2+b2

  plot(Delta, delta_pdf(Delta, a1, b1, a2, b2), 
       main = substitute(paste(n[1], "=", n1, ", ", n[2], "=", n2, ": ", p[1] %~% Beta (a1, b1), ", ", p[1] %~% Beta (a2, b2)), list(n1 = n1, n2 = n2, a1 = round(a1, 3), b1 = round(b1, 3), a2 = round(a2, 3), b2 = round(b2, 3))), 
       xlab = expression(Delta == p[E]-p[C]),
       ylab = "p. d. f. of p ( f(p) )",
       ylim = c(0, 6.5),
       type="l")
  grid()
} 
}, interval = 1, movie.name = 'bayesian_demo.gif', ani.width = 600, ani.height = 600)


## ===== Simulations =====
# CRD vs. BAR 

CRD <- function(nsbj, Y) {
  delta <- vector("integer", nsbj)      # vector of treatment assignments
  resp <- vector("integer", nsbj)       # vector of responses
  
  for (j in seq_along(delta)) {
    delta[j] <- as.integer(runif(1) <= 0.5)
    resp[j] <- delta[j]*Y$E[j] + (1-delta[j])*Y$C[j]
  }
  pE_hat <- sum(delta*resp)/sum(delta)
  pC_hat <- sum((1-delta)*resp)/sum(1-delta)
  
  # output
  list(
    pE = pE_hat,            # estimated probability success on E
    pC = pC_hat,            # estimated probability success on C
    nE = sum(delta),        # NE(n) -- number of subjects assigned to E
    prop = sum(delta)/nsbj, # NE(n)/n -- proportion of subjects assigned to E
    TNF = sum(1-resp),      # Total Number of Failures, TNF(n)
    reject = as.integer(    # rejection of H0: pE = pC
      abs(pE_hat-pC_hat)/sqrt(pE_hat*(1-pE_hat)/sum(delta) + pC_hat*(1-pC_hat)/sum(1-delta)) > qnorm(0.975)
    )
  )
}


BAR <- function(nsbj, Y, lambda = 1, m0 = 10) {
  delta <- vector("integer", nsbj)      # vector of treatment assignments
  resp <- vector("integer", nsbj)       # vector of responses
  
  a1 <- 1
  b1 <- 1
  
  a2 <- 1
  b2 <- 1
  
  for (j in seq_along(delta)) {
    # PBD
    if (j <= m0) prob <- (m0/2-sum(delta[seq_len(j-1)]))/(m0-(j-1))
    # BAR
    else {
      Pi <- delta_cdf(0, a1, b1, a2, b2)
      # lambda = 1 => Thompson
      if (class(lambda) == "function") {
        lmb <- lambda(j, nsbj)
      }
      else {
        lmb <- lambda
      }
      prob <- Pi^lmb/(Pi^lmb + abs(1-Pi)^lmb)
    }
    
      
    delta[j] <- as.integer(runif(1) <= prob)
    resp[j] <- delta[j]*Y$E[j] + (1-delta[j])*Y$C[j]
    
    X1 <- sum(resp[seq_len(j)]*delta[seq_len(j)])
    X2 <- sum(resp[seq_len(j)]*(1-delta[seq_len(j)]))
    
    a1 <- X1 + a1
    b1 <- j - X1 + b1
    
    a2 <- X2 + a1
    b2 <- j - X2 + b2
  }
  
  pE_hat <- sum(delta*resp)/sum(delta)
  pC_hat <- sum((1-delta)*resp)/sum(1-delta)
  
  # output
  list(
    pE = pE_hat,            # estimated probability success on E
    pC = pC_hat,            # estimated probability success on C
    nE = sum(delta),        # NE(n) -- number of subjects assigned to E
    prop = sum(delta)/nsbj, # NE(n)/n -- proportion of subjects assigned to E
    TNF = sum(1-resp),      # Total Number of Failures, TNF(n)
    reject = as.integer(    # rejection of H0: pE = pC
      abs(pE_hat-pC_hat)/sqrt(pE_hat*(1-pE_hat)/sum(delta) + pC_hat*(1-pC_hat)/sum(1-delta)) > qnorm(0.975)
    )
  )
}

DBCD_BAR <- function(nsbj, Y, xi0, lambda = 0.5, gm = 2, m0 = 10) {
  delta <- vector("integer", nsbj)      # vector of treatment assignments
  resp <- vector("integer", nsbj)       # vector of responses
  
  a1 <- 1
  b1 <- 1
  
  a2 <- 1
  b2 <- 1
  
  for (j in seq_along(delta)) {
    # PBD
    if (j <= m0) prob <- (m0/2-sum(delta[seq_len(j-1)]))/(m0-(j-1))
    # BAR
    else {
      Pi <- delta_cdf(0, a1, b1, a2, b2)
      # lambda = 1 => Thompson
      if (class(lambda) == "function") {
        lmb <- lambda(j, nsbj)
      }
      else {
        lmb <- lambda
      }
      rho <- Pi^lmb/(Pi^lmb + abs(1-Pi)^lmb)
      if (rho <= 1-xi0) {rho_ <- 1-xi0}
      else if (rho >= xi0) {rho_ <- xi0}
      else {rho_ <- rho}
      
      prob <- rho_*(rho_/(sum(delta[seq_len(j-1)])/(j-1)))^gm/
        (rho_*(rho_/(sum(delta[seq_len(j-1)])/(j-1)))^gm + (1-rho_)*((1-rho_)/(sum(1-delta[seq_len(j-1)])/(j-1)))^gm)
    }
    
    
    delta[j] <- as.integer(runif(1) <= prob)
    resp[j] <- delta[j]*Y$E[j] + (1-delta[j])*Y$C[j]
    
    X1 <- sum(resp[seq_len(j)]*delta[seq_len(j)])
    X2 <- sum(resp[seq_len(j)]*(1-delta[seq_len(j)]))
    
    a1 <- X1 + a1
    b1 <- j - X1 + b1
    
    a2 <- X2 + a1
    b2 <- j - X2 + b2
    # print(a1)
    # print(b1)
    # print(a2)
    # print(b2)
    # print("=======")
  }
  
  pE_hat <- sum(delta*resp)/sum(delta)
  pC_hat <- sum((1-delta)*resp)/sum(1-delta)
  
  # output
  list(
    pE = pE_hat,            # estimated probability success on E
    pC = pC_hat,            # estimated probability success on C
    nE = sum(delta),        # NE(n) -- number of subjects assigned to E
    prop = sum(delta)/nsbj, # NE(n)/n -- proportion of subjects assigned to E
    TNF = sum(1-resp),      # Total Number of Failures, TNF(n)
    reject = as.integer(    # rejection of H0: pE = pC
      abs(pE_hat-pC_hat)/sqrt(pE_hat*(1-pE_hat)/sum(delta) + pC_hat*(1-pC_hat)/sum(1-delta)) > qnorm(0.975)
    )
  )
}
## Uncomment the code below, if there are no corresponding outputs in ./figures and ./data

# nsim <- 10000              # number of simulations
# nsbj <- 200                # number of subjects
# Y_df <- map(seq_len(nsim), ~
#               data_frame(
#                 E = rbinom(nsbj, 1, 0.4),
#                 C = rbinom(nsbj, 1, 0.2)
#               )
# )                          # subjects' responses
# 
# crd <- map(Y_df, ~ CRD(nsbj, .))
# bar1 <- map(Y_df, ~ BAR(nsbj, .))
# bar05 <- map(Y_df, ~ BAR(nsbj, ., lambda = 0.5))
# bar_lmb <- map(Y_df, ~ BAR(nsbj, ., lambda = function(m, N) m/(2*N)))
# 
# # DBCD
# dbcd_bar_080 <- map(Y_df, ~ DBCD_BAR(nsbj, ., xi0 = 0.80))
# dbcd_bar_075 <- map(Y_df, ~ DBCD_BAR(nsbj, ., xi0 = 0.75))
# dbcd_bar_066 <- map(Y_df, ~ DBCD_BAR(nsbj, ., xi0 = 2/3))
# 
# ## ===== Analysis =====
# # CRD vs. Thompson summary
# bar_summary <- map2(
#   list(crd, bar1),
#   list("1. CRD", "2. BAR (Thompson)"),
#   ~ {
#     df <- .x
#     design <- .y
#     df %>% map(~ {
#       data_frame(
#         design = design,
#         nE = .$nE,
#         prop = .$prop,
#         TNF = .$TNF,
#         reject = .$reject
#       )
#     }) %>%
#       bind_rows()
#   }
# ) %>%
#   bind_rows()
# 
# # boxplots
# bar_summary %>% 
#   select(design, `Allocation proportion` = prop, `Total number of Failures` = TNF) %>% 
#   gather(key, value, -design) %>%
#   ggplot(aes(x = factor(design), y = value))+
#     geom_boxplot()+
#     xlab("")+
#     ylab("")+
#     facet_wrap( ~ key, scales = "free_y")+
#     theme(
#       axis.text = element_text(family = "Helvetica", face = "bold", size = 14),
#       axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
#       strip.text = element_text(family = "Helvetica", face = "bold", size = 18)
#     ) 
# 
# ggsave("./figures/lecture07-fig08-crd-vs-thompson.png", width = 16, height = 9, units = "in")
# 
# bar_summary %>% 
#   group_by(design) %>% 
#   summarise(
#     mean_nE = round(mean(nE)),
#     min_nE = round(min(nE)),
#     max_nE = round(max(nE)),
#     mean_prop = mean(prop),
#     var_prop = var(prop),
#     sd_prop = sd(prop),
#     Power = mean(reject),
#     TNF_mean = mean(TNF),
#     TNF_sd = sd(TNF)
#   ) %>%
#   mutate(
#     mean_prop = round(mean_prop, 2),
#     var_prop = round(var_prop, 4),
#     Power = round(Power, 2),
#     TNF_mean = round(TNF_mean),
#     TNF_sd = round(TNF_sd)
#   )
# #####  
# 
# # Exploring the Impact of lambda
# bar_summary2 <- map2(
#   list(crd, bar1, bar05, bar_lmb),
#   list("1. CRD", "2. lmb = 1 (Thompson)", "3. lmb = 1/2", "4. lmb = m/(2n)"),
#   ~ {
#     df <- .x
#     design <- .y
#     df %>% map(~ {
#       data_frame(
#         design = design,
#         nE = .$nE,
#         prop = .$prop,
#         TNF = .$TNF,
#         reject = .$reject
#       )
#     }) %>%
#       bind_rows()
#   }
# ) %>%
#   bind_rows()
# 
# # boxplots
# bar_summary2 %>% 
#   select(design, `Allocation proportion` = prop, `Total number of Failures` = TNF) %>% 
#   gather(key, value, -design) %>%
#   ggplot(aes(x = factor(design), y = value))+
#   geom_boxplot()+
#   xlab("")+
#   ylab("")+
#   facet_wrap( ~ key, scales = "free_y")+
#   theme(
#     axis.text = element_text(family = "Helvetica", face = "bold", size = 14),
#     axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
#     strip.text = element_text(family = "Helvetica", face = "bold", size = 18)
#   ) 
# 
# ggsave("./figures/lecture07-fig09-bar-vs-lambda.png", width = 16, height = 9, units = "in")
# 
# bar_summary2 %>% 
#   group_by(design) %>% 
#   summarise(
#     mean_nE = round(mean(nE)),
#     min_nE = round(min(nE)),
#     max_nE = round(max(nE)),
#     mean_prop = mean(prop),
#     var_prop = var(prop),
#     sd_prop = sd(prop),
#     Power = mean(reject),
#     TNF_mean = mean(TNF),
#     TNF_sd = sd(TNF)
#   ) %>%
#   mutate(
#     mean_prop = round(mean_prop, 2),
#     var_prop = round(var_prop, 4),
#     Power = round(Power, 2),
#     TNF_mean = round(TNF_mean),
#     TNF_sd = round(TNF_sd)
#   )
# 
# save(bar_summary, bar_summary2, file = "./data/lecture07-bar-summary.Rda")

#####

