library(tidyverse)
library(nleqslv)

## ===== Different Optimal Allocation function ===== ##

Neyman <- function(pE, pC) {
  sqrt(pE*(1-pE))/(sqrt(pE*(1-pE)) + sqrt(pC*(1-pC)))
}

RSIHR <- function(pE, pC) {
  sqrt(pE)/(sqrt(pE) + sqrt(pC))
}

urn <- function(pE, pC) {
  (1-pC)/(2-pE-pC)
}

compound <- function(pE, pC, omega = 0.5) {
  eq <- function(x, p1, p2, q1, q2, omega_) {
    omega_/(1-omega_)*(p1-p2)/min(q1, q2)*(sqrt(p2*q2/(p1*q1))+1)^2-
      ((p2*q2/(p1*q1)-1)*x^2+2*x-1)/(x*(1-x))^2
  }
  nleqslv(0.5, eq, p1 = pE, p2 = pC, q1 = 1-pE, q2 = 1-pC, omega_ = omega)$x
} 

# ===== functions to calculate power and ENF (expected number of tretament failures) ===== #
## ===== ===== ##

Power <- function(pE, pC, rho, n, alpha = 0.05) {
  Delta <- (pE-pC)/sqrt(pE*(1-pE)/(n*rho)+pC*(1-pC)/(n*(1-rho)))
  1-pnorm(qnorm(1-alpha/2) - abs(Delta))
}

ENF <- function(pE, pC, rho, n) {
  n*((1-pE)*rho + (1-pC)*(1-rho))
}


## ===== Compare different allocations ===== ##
## Uncomment the code below, if there are no outputs in ./figures and ./data


# cmp_allocations <- function(pE = seq(0.2, 0.9, by = 0.05), pC = 0.4, omega = 0.5) {
#   tibble::tibble(
#     `pE` = pE,
#     `1) 1:1` = 0.5,
#     `2) 2:1` = 2/3,
#     `3) Neyman` = map_dbl(pE, ~ Neyman(., pC)),
#     `4) RSIHR` = map_dbl(pE, ~ RSIHR(., pC)),
#     `5) Urn` = map_dbl(pE, ~ urn(., pC)),
#     `6) Compound` = map_dbl(pE, ~ compound(., pC, omega))
#   ) %>% 
#     gather(allocation, value, -pE) %>% 
#     ggplot(aes(x = pE, y = value, color = allocation, linetype = allocation))+
#       geom_line(size = 1.25)+
#       scale_x_continuous(limits = c(0.2, 0.9), breaks = seq(0.2, 0.9, by = 0.1))+
#       scale_y_continuous(limits = c(0.2, 0.9), breaks = seq(0.1, 0.9, by = 0.1))+
#       scale_color_manual(values = c("black", "darkgrey", "green", "darkblue", "red", "violet"))+
#       scale_linetype_manual(values = c("solid", "solid", "solid", "longdash", "dashed", "dotdash"))+
#       labs(
#         title = "Plot of Allocation Proportion to E",
#         subtitle = "(Probability of success on C = 0.40)",
#         x = "Probability of success on E",
#         y = "Allocation proportion"
#       )+
#       theme(
#         title = element_text(family = "Helvetica", face = "bold", size = 16),
#         axis.title = element_text(family = "Helvetica", face = "bold", size = 14),
#         axis.text  = element_text(family = "Helvetica", face = "bold", size = 14),
#         legend.title = element_blank(),
#         legend.text = element_text(family = "Helvetica", face = "bold", size = 14),
#         legend.position = c(0.3, 0.8),
#         legend.background = element_rect(fill = NA)
#       )
#   ggsave("./figures/lecture06-fig02-allocation.png", width = 8, height = 8, units = "in")    
# }
# 
# ## ===== Power [ENF] table ===== ##
# 
# pE <- c(0.2, 0.5, 0.6, 0.7, 0.8)
# pC <- rep(0.4, length(pE))
# 
# tibble::tibble(
#   pE,
#   pC,
#   `1:1` = paste0(
#     map2_dbl(pE, pC, ~ round(Power(.x, .y, 0.5, 100), 2)),
#     " [",
#     map2_dbl(pE, pC, ~ round(ENF(.x, .y, 0.5, 100))),
#     "]"
#   ),
#   `2:1` = paste0(
#     map2_dbl(pE, pC, ~ round(Power(.x, .y, 2/3, 100), 2)),
#     " [",
#     map2_dbl(pE, pC, ~ round(ENF(.x, .y, 2/3, 100))),
#     "]"
#   ),
#   Neyman = paste0(
#     map2_dbl(pE, pC, ~ round(Power(.x, .y, Neyman(.x, .y), 100), 2)),
#     " [",
#     map2_dbl(pE, pC, ~ round(ENF(.x, .y, Neyman(.x, .y), 100))),
#     "]"
#   ),
#   RSIHR = paste0(
#     map2_dbl(pE, pC, ~ round(Power(.x, .y, RSIHR(.x, .y), 100), 2)),
#     " [",
#     map2_dbl(pE, pC, ~ round(ENF(.x, .y, RSIHR(.x, .y), 100))),
#     "]"
#   ),
#   Urn = paste0(
#     map2_dbl(pE, pC, ~ round(Power(.x, .y, urn(.x, .y), 100), 2)),
#     " [",
#     map2_dbl(pE, pC, ~ round(ENF(.x, .y, urn(.x, .y), 100))),
#     "]"
#   ),
#   `CO ($\\omega$ = 0.5)` = paste0(
#     map2_dbl(pE, pC, ~ round(Power(.x, .y, compound(.x, .y, 0.5), 100), 2)),
#     " [",
#     map2_dbl(pE, pC, ~ round(ENF(.x, .y, compound(.x, .y, 0.5), 100))),
#     "]"
#   )
# )
# 
# 
# ## ===== function to calculate asymptotic variances of allocation proportion =====
# 
# # Target allocation is RSIHR
# # Procedures are: SMLE, DBCD, ERADE
# sigma2_fcn <- function(pE, pC, procedure, parameter=NA) {
#   # lower bound for the variance of RAR procedure targeting RSIHR
#   B <- 1/(4*(sqrt(pE)+sqrt(pC))^3)*((pC*(1-pE))/sqrt(pE)+(pE*(1-pC))/sqrt(pC))
#   
#   if (procedure == "SMLE") {
#     sqrt(pE*pC)/(sqrt(pE)+sqrt(pC))^2 + 2*B
#   }
#   else if (procedure == "DBCD") {
#     gm <- parameter
#     sqrt(pE*pC)/((1+2*gm)*(sqrt(pE)+sqrt(pC))^2) + 2*(1+gm)/(1+2*gm)*B
#   }
#   else if (procedure == "ERADE") {
#     B
#   }
#   else {
#     stop("the value of input parameter 'procedure' must be one of c('SMLE', 'DBCD', 'ERADE')")
#   }
# }
# 
# ## ===== Compare variances ===== ##
# cmp_variances <- function(pE = seq(0.2, 0.9, by = 0.05), pC = 0.4) {
#   tibble::tibble(
#     `pE` = pE,
#     `1) 1:1` = 0.25,
#     `2) SMLE` = map_dbl(pE, ~ sigma2_fcn(., pC, "SMLE")),
#     `3) DBCD (2)` = map_dbl(pE, ~ sigma2_fcn(., pC, "DBCD", 2)),
#     `4) ERADE` = map_dbl(pE, ~ sigma2_fcn(., pC, "ERADE"))
#   ) %>% 
#     gather(allocation, value, -pE) %>% 
#     ggplot(aes(x = pE, y = value, color = allocation, linetype = allocation))+
#     geom_line(size = 1.25)+
#     scale_x_continuous(limits = c(0.2, 0.9), breaks = seq(0.2, 0.9, by = 0.1))+
#     scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))+
#     scale_color_manual(values = c("black", "darkblue", "red", "violet"))+
#     scale_linetype_manual(values = c("solid", "longdash", "dashed", "dotdash"))+
#     labs(
#       title = "Asymptotic Variance of the Sample",
#       subtitle = "(Probability of success on C = 0.40)",
#       x = "Probability of success on E",
#       y = "Asymptotic variance"
#     )+
#     theme(
#       title = element_text(family = "Helvetica", face = "bold", size = 16),
#       axis.title = element_text(family = "Helvetica", face = "bold", size = 14),
#       axis.text  = element_text(family = "Helvetica", face = "bold", size = 14),
#       legend.title = element_blank(),
#       legend.text = element_text(family = "Helvetica", face = "bold", size = 14),
#       legend.position = c(0.3, 0.8),
#       legend.background = element_rect(fill = NA)
#     )
#   ggsave("./figures/lecture06-fig03-variance.png", width = 8, height = 8, units = "in")    
# }
# 
# 
# ## ====== simulations ===== ##
# 
# # Complete randomization (CRD)
# CRD <- function(nsbj, Y) {
#   # setup input
#   ntrt <- ncol(Y)                       # number of treatments  
#   rho <- rep(1, ntrt)/ntrt              # target allocation proportions
#   delta <- vector("integer", nsbj)      # vector of treatment assignments
#   resp <- vector("integer", nsbj)       # vector of responses
#   
#   for (j in seq_along(delta)) {
#     delta[j] <- as.integer(runif(1) <= rho[1])
#     resp[j] <- delta[j]*Y$E[j] + (1-delta[j])*Y$C[j]
#   }
#   pE_hat <- sum(delta*resp)/sum(delta)
#   pC_hat <- sum((1-delta)*resp)/sum(1-delta)
#   
#   # output
#   list(
#     pE = pE_hat,            # estimated probability success on E
#     pC = pC_hat,            # estimated probability success on C
#     nE = sum(delta),        # NE(n) -- number of subjects assigned to E
#     prop = sum(delta)/nsbj, # NE(n)/n -- proportion of subjects assigned to E
#     TNF = sum(1-resp),      # Total Number of Failures, TNF(n)
#     reject = as.integer(    # rejection of H0: pE = pC
#       abs(pE_hat-pC_hat)/sqrt(pE_hat*(1-pE_hat)/sum(delta) + pC_hat*(1-pC_hat)/sum(1-delta)) > qnorm(0.975)
#     )
#   )
# }
# 
# # Permuted Block (CRD)
# # b -- block size
# PBD <- function(b, nsbj, Y) {
#   # setup input
#   ntrt <- ncol(Y)                       # number of treatments  
#   rho <- rep(1, ntrt)/ntrt              # target allocation proportions
#   delta <- vector("integer", nsbj)      # vector of treatment assignments
#   resp <- vector("integer", nsbj)       # vector of responses
#   
#   delta_ <- vector("integer", b)
#   
#   for (j in seq_along(delta)) {
#     m <- j %% b + b*as.integer(j %% b == 0) # sbject's ID in a block of size b
#     
#     if (m == 1) {
#       prob <- 0.5
#     }
#     else {
#       prob <- (b/2-sum(delta_[seq_len(m-1)]))/(b-(m-1))
#     }
#     
#     delta[j] <- as.integer(runif(1) <= prob)
#     delta_[m] <- delta[j]
#     
#     resp[j] <- delta[j]*Y$E[j] + (1-delta[j])*Y$C[j]
#     
#     if (m == b) {
#       delta_ <- vector("integer", b)
#     }
#   }
#   pE_hat <- sum(delta*resp)/sum(delta)
#   pC_hat <- sum((1-delta)*resp)/sum(1-delta)
#   
#   # output
#   list(
#     pE = pE_hat,            # estimated probability success on E
#     pC = pC_hat,            # estimated probability success on C
#     nE = sum(delta),        # NE(n) -- number of subjects assigned to E
#     prop = sum(delta)/nsbj, # NE(n)/n -- proportion of subjects assigned to E
#     TNF = sum(1-resp),      # Total Number of Failures, TNF(n)
#     reject = as.integer(    # rejection of H0: pE = pC
#       abs(pE_hat-pC_hat)/sqrt(pE_hat*(1-pE_hat)/sum(delta) + pC_hat*(1-pC_hat)/sum(1-delta)) > qnorm(0.975)
#     )
#   )
# }
# 
# 
# # RAR targeting optimal allocation with ERADE (nu)
# RAR <- function(allocation_fcn, m0, nsbj, Y, nu = 0.5) {
#   # setup input
#   ntrt <- ncol(Y)                       # number of treatments  
#   rho <- rep(1, ntrt)/ntrt              # target allocation proportions
#   delta <- vector("integer", nsbj)      # vector of treatment assignments
#   resp <- vector("integer", nsbj)       # vector of responses
#   
#   delta_ <- vector("integer", 2*m0)
#   
#   # first 2m0 subjects are randomized with PBD (block size = 2m0)
#   for (j in seq_len(2*m0)) {
#     m <- j %% (2*m0) + 2*m0*as.integer(j %% (2*m0) == 0) # sbject's ID in a block of size 2*m0
#     
#     if (m == 1) {
#       prob <- 0.5
#     }
#     else {
#       prob <- (m0-sum(delta_[seq_len(m-1)]))/(2*m0-(m-1))
#     }
#     delta_[m] <- as.integer(runif(1) <= prob)
#     
#     delta[j] <- delta_[m]
#     resp[j] <- delta[j]*Y$E[j] + (1-delta[j])*Y$C[j]
#   }
#   
#   # the rest subijects are rndomized with RAR
#   for (j in (seq_len(nsbj-2*m0)+2*m0)) {
#     # Bayesian estimates of success rates
#     pE_hat <- (sum(delta[seq_len(j-1)]*resp[seq_len(j-1)])+0.5)/(sum(delta[seq_len(j-1)])+1)
#     pC_hat <- (sum((1-delta[seq_len(j-1)])*resp[seq_len(j-1)])+0.5)/(sum(1-delta[seq_len(j-1)])+1)
# 
#     # updated target allocation
#     rho_hat <- allocation_fcn(pE_hat, pC_hat)
#     
#     # ERADE (nu)
#     prop <- sum(delta[seq_len(j-1)])/(j-1)
#     if (prop > rho_hat) {
#       prob <- nu*rho_hat
#     }
#     else if (prop < rho_hat) {
#       prob <- 1-nu*(1-rho_hat)
#     }
#     else {
#       prob <- rho_hat
#     }
#     
#     delta[j] <- as.integer(runif(1) <= prob)
#     resp[j] <- delta[j]*Y$E[j] + (1-delta[j])*Y$C[j]
#     
#   }
#   pE_hat <- sum(delta*resp)/sum(delta)
#   pC_hat <- sum((1-delta)*resp)/sum(1-delta)
#   
#   # output
#   list(
#     pE = pE_hat,            # estimated probability success on E
#     pC = pC_hat,            # estimated probability success on C
#     nE = sum(delta),        # NE(n) -- number of subjects assigned to E
#     prop = sum(delta)/nsbj, # NE(n)/n -- proportion of subjects assigned to E
#     TNF = sum(1-resp),      # Total Number of Failures, TNF(n)
#     reject = as.integer(    # rejection of H0: pE = pC
#       abs(pE_hat-pC_hat)/sqrt(pE_hat*(1-pE_hat)/sum(delta) + pC_hat*(1-pC_hat)/sum(1-delta)) > qnorm(0.975)
#     )
#   )
# }
# 
# 
# nsim <- 10000              # number of simulations
# nsbj <- 120                # number of subjects
# Y_df <- map(seq_len(nsim), ~
#   tibble::tibble(
#     E = rbinom(nsbj, 1, 0.7),
#     C = rbinom(nsbj, 1, 0.4)
#   )
# )                          # subjects' responses
# 
# crd <- map(Y_df, ~ CRD(nsbj, .))
# pbd <- map(Y_df, ~ PBD(4, nsbj, .))
# rar_neyman <- map(Y_df, ~ RAR(Neyman, 5, nsbj, .))
# rar_rsihr <- map(Y_df, ~ RAR(RSIHR, 5, nsbj, .))
# rar_urn <- map(Y_df, ~ RAR(urn, 5, nsbj, .))
# rar_compound <- map(Y_df, ~ RAR(compound, 5, nsbj, .))
# 
# # # summary of simulation data
# rar_summary <- map2(
#   list(crd, pbd, rar_neyman, rar_rsihr, rar_urn, rar_compound),
#   list("1. CRD", "2. PBD", "3. Neyman", "4. RSIHR", "5. Urn", "6. CO(omega=0.5)"),
#   ~ {
#     df <- .x
#     design <- .y
#     df %>% map(~ {
#       tibble::tibble(
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
#   bind_rows() %>%
#   group_by(design) %>%
#   summarise(
#     mean_nE = round(mean(nE)),
#     min_nE = round(min(nE)),
#     max_nE = round(max(nE)),
#     mean_prop = mean(prop),
#     var_prop = var(prop),
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
# # # boxplots of allocation proportions
# rar_boxplot <- map2(
#   list(crd, pbd, rar_neyman, rar_rsihr, rar_urn, rar_compound),
#   list("1. CRD", "2. PBD", "3. Neyman", "4. RSIHR", "5. Urn", "6. CO(omega=0.5)"),
#   ~ {
#     df <- .x
#     design <- .y
#     df %>% map(~ {
#       tibble::tibble(
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
#   bind_rows() %>%
#   select(design, prop) %>%
#   ggplot(aes(x = factor(design), y = prop))+
#     geom_boxplot()
# 
# rar_boxplot +
#   labs(
#     title = "Distribution of Allocation Proportion on E",
#     x = "",
#     y = "Allocation proportion"
#   )+
#   theme(
#     title = element_text(family = "Helvetica", face = "bold", size = 16),
#     axis.title = element_text(family = "Helvetica", face = "bold", size = 14),
#     axis.text  = element_text(family = "Helvetica", face = "bold", size = 14),
#     legend.title = element_blank(),
#     legend.text = element_text(family = "Helvetica", face = "bold", size = 14),
#     legend.position = c(0.3, 0.8),
#     legend.background = element_rect(fill = NA)
#   )
# ggsave("./figures/lecture06-fig04-sim-allocation.png", width = 8, height = 8, units = "in")


## ===== Redesign of the AZT trial ===== ##
# nsim <- 10000              # number of simulations
# nsbj <- 477                # number of subjects
# Y_df <- map(seq_len(nsim), ~
#               tibble::tibble(
#                 E = rbinom(nsbj, 1, 0.917), # AZT has 0.917 probaility of success
#                 C = rbinom(nsbj, 1, 0.745)  # PBL has 0.745 probaility of success
#               )
# )                          # subjects' responses
# 
# AZT_pbd <- map(Y_df, ~ PBD(4, nsbj, .))
# AZT_rar_compound <- map(Y_df, ~ RAR(compound, 5, nsbj, .))
# 
# 
# AZT_summary <- map2(
#   list(AZT_pbd, AZT_rar_compound),
#   list("1. PBD", "2. CO (omega=0.5) RAR"),
#   ~ {
#     df <- .x
#     design <- .y
#     df %>% map(~ {
#       tibble::tibble(
#         design = design,
#         nE = .$nE,
#         pE = .$pE,
#         pC = .$pC,
#         prop = .$prop,
#         TNF = .$TNF,
#         reject = .$reject
#       )
#     }) %>%
#       bind_rows()
#   }
# ) %>%
#   bind_rows() %>%
#   group_by(design) %>%
#   summarise(
#     Mean = round(mean(prop), 3),
#     SD = round(sd(prop), 3),
#     ave_nE = round(mean(nE)),
#     min_nE = round(min(nE)),
#     max_nE = round(max(nE)),
#     ave_TNF = round(mean(TNF)),
#     min_TNF = min(TNF),
#     max_TNF = max(TNF),
#     ave_pE = round(mean(pE), 3),
#     sd_pE = round(sd(pE), 3),
#     ave_pC = round(mean(pC), 3),
#     sd_pC = round(sd(pC), 3)
#   ) %>%
#   mutate(
#     Prop = paste0(Mean, " (", SD, ")"),
#     nE = paste0(ave_nE, " [", min_nE, " - ", max_nE, "]"),
#     TNF = paste0(ave_TNF, " [", min_TNF, " - ", max_TNF, "]"),
#     pE = paste0(ave_pE, " (", sd_pE, ")"),
#     pC = paste0(ave_pC, " (", sd_pC, ")")
#   ) %>%
#   select(Design = design, Prop, nE, TNF, pE, pC)
# 
# save(AZT_pbd, AZT_rar_compound, AZT_summary, file = "./data/lecture06-AZT.Rda")
