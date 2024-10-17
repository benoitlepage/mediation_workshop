### test programs for estimations based on IPTW

rm(list=ls())

df2 <- read.csv(file = "./data/df2.csv")
df2_int <- read.csv(file = "./data/df2_int.csv")


# ---------------------------------------------------------------------------- #
# 1) Estimation of the Average Total Effect (ATE) ----
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
## 1.1) IPTW for the ATE ----
# ---------------------------------------------------------------------------- #
## 1. Estimate g
g.L <- glm(A0_PM2.5 ~ L0_male + L0_soc_env,
           family = "binomial", data = df2_int)

## 2. Predict each individual's probability of being exposed to her own exposure
# predict the probabilities P(A0_PM2.5=1|L(0)) & P(A0_PM2.5=0|L(0))
pred.g1.L <- predict(g.L, type="response")
pred.g0.L <- 1 - pred.g1.L
# the predicted probability of the observed treatment A=a_i is :
gA.L <- rep(NA, nrow(df2_int))
gA.L[df2_int$A0_PM2.5 == 1] <- pred.g1.L[df2_int$A0_PM2.5 == 1]
gA.L[df2_int$A0_PM2.5 == 0] <- pred.g0.L[df2_int$A0_PM2.5 == 0]

## 3. Apply weights corresponding to the inverse of the predicted probability
wt <- 1 / gA.L

## 4. Use the empirical mean of the weighted outcome
# point estimates:
IPTW.death <- mean(wt * as.numeric(df2_int$A0_PM2.5 == 1) * df2_int$Y_death) -
  mean(wt * as.numeric(df2_int$A0_PM2.5 == 0) * df2_int$Y_death)
IPTW.death
# [1] 0.08224947

IPTW.qol <- mean(wt * as.numeric(df2_int$A0_PM2.5 == 1) * df2_int$Y_qol) -
  mean(wt * as.numeric(df2_int$A0_PM2.5 == 0) * df2_int$Y_qol)
IPTW.qol
# [1] -8.436797


# ---------------------------------------------------------------------------- #
## 1.2) Stabilized IPTW for the ATE ----
# ---------------------------------------------------------------------------- #

## 3. For example, applying g^*(A) = 1
## 4. Applying the stabilized estimator
# point estimates:
s.IPTW.death <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1) * df2_int$Y_death) /
                   mean(wt * as.numeric(df2_int$A0_PM2.5 == 1))) -
  (mean(wt * as.numeric(df2_int$A0_PM2.5 == 0) * df2_int$Y_death) /
     mean(wt * as.numeric(df2_int$A0_PM2.5 == 0)))
s.IPTW.death
# [1] 0.08294185

s.IPTW.qol <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1) * df2_int$Y_qol) /
                 mean(wt * as.numeric(df2_int$A0_PM2.5 == 1))) -
  (mean(wt * as.numeric(df2_int$A0_PM2.5 == 0) * df2_int$Y_qol) /
     mean(wt * as.numeric(df2_int$A0_PM2.5 == 0)))
s.IPTW.qol
# [1] -8.291992


# ---------------------------------------------------------------------------- #
# 2) Estimation of the Controlled direct effect (CDE) ----
# ---------------------------------------------------------------------------- #
## 2.1) IPTW for the CDE ----
# ---------------------------------------------------------------------------- #
## 1. Estimate gA and gM
gA.L <- glm(A0_PM2.5 ~ L0_male + L0_soc_env,
            family = "binomial", data = df2_int)
gM.L <- glm(M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
            family = "binomial", data = df2_int)

## 2. Predict each individual's probability of being exposed to her own exposure
# predict the probabilities P(A0_PM2.5=1|L(0)) & P(A0_PM2.5=0|L(0))
pred.gA1.L <- predict(gA.L, type = "response")
pred.gA0.L <- 1 - pred.gA1.L
# the predicted probability of the observed treatment A_i=a is :
gAobs.L <- rep(NA, nrow(df2_int))
gAobs.L[df2_int$A0_PM2.5 == 1] <- pred.gA1.L[df2_int$A0_PM2.5 == 1]
gAobs.L[df2_int$A0_PM2.5 == 0] <- pred.gA0.L[df2_int$A0_PM2.5 == 0]

# predict the probabilities P(M=1|L(0),A,L(1)) & P(M=0|L(0),A,L(1))
pred.gM1.L <- predict(gM.L, type = "response")
pred.gM0.L <- 1 - pred.gM1.L
# the predicted probability of the observed treatment M_i=m is :
gMobs.L <- rep(NA, nrow(df2_int))
gMobs.L[df2_int$M_diabetes == 1] <- pred.gM1.L[df2_int$M_diabetes == 1]
gMobs.L[df2_int$M_diabetes == 0] <- pred.gM0.L[df2_int$M_diabetes == 0]

## 3. Apply weights corresponding to the inverse of the predicted probability
wt_A <- 1 / gAobs.L
wt_M <- 1 / gMobs.L
wt <- wt_A * wt_M

## 4. Use the empirical mean of the weighted outcome
# point estimates of CDE, setting M=0
CDE_IPTW_m0_death <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                             df2_int$M_diabetes == 0) *
                             df2_int$Y_death) -
                        mean(wt * as.numeric(df2_int$A0_PM2.5==0 &
                                               df2_int$M_diabetes == 0) *
                               df2_int$Y_death))
CDE_IPTW_m0_death
# [1] 0.05874684

CDE_IPTW_m0_qol <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                           df2_int$M_diabetes == 0) *
                           df2_int$Y_qol) -
                      mean(wt * as.numeric(df2_int$A0_PM2.5==0 &
                                             df2_int$M_diabetes == 0) *
                             df2_int$Y_qol))
CDE_IPTW_m0_qol
# [1] -5.341138

# point estimates of CDE, setting M=1
CDE_IPTW_m1_death <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                             df2_int$M_diabetes == 1) *
                             df2_int$Y_death) -
                        mean(wt * as.numeric(df2_int$A0_PM2.5==0 &
                                               df2_int$M_diabetes == 1) *
                               df2_int$Y_death))
CDE_IPTW_m1_death
# [1] 0.101733

CDE_IPTW_m1_qol <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                           df2_int$M_diabetes == 1) *
                           df2_int$Y_qol) -
                      mean(wt * as.numeric(df2_int$A0_PM2.5==0 &
                                             df2_int$M_diabetes == 1) *
                             df2_int$Y_qol))
CDE_IPTW_m1_qol
# [1] -8.185866

# ---------------------------------------------------------------------------- #
## 2.2) Stabilized IPTW for the CDE ----
# ---------------------------------------------------------------------------- #
## 4. Applying the stabilized estimator
# point estimates of CDE, setting M=0:
CDE_sIPTW_m0_death <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                              df2_int$M_diabetes == 0) *
                              df2_int$Y_death) /
                         mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                                df2_int$M_diabetes == 0))) -
  (mean(wt * as.numeric(df2_int$A0_PM2.5 == 0 &
                          df2_int$M_diabetes == 0) *
          df2_int$Y_death) /
     mean(wt * as.numeric(df2_int$A0_PM2.5 == 0 &
                            df2_int$M_diabetes == 0)))
CDE_sIPTW_m0_death
# [1] 0.0601292

CDE_sIPTW_m0_qol <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                            df2_int$M_diabetes == 0) *
                            df2_int$Y_qol) /
                       mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                              df2_int$M_diabetes == 0))) -
  (mean(wt * as.numeric(df2_int$A0_PM2.5 == 0 &
                          df2_int$M_diabetes == 0) *
          df2_int$Y_qol) /
     mean(wt * as.numeric(df2_int$A0_PM2.5 == 0 &
                            df2_int$M_diabetes == 0)))
CDE_sIPTW_m0_qol
# [1] -4.966328

# point estimates of CDE, setting M=1:
CDE_sIPTW_m1_death <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                              df2_int$M_diabetes == 1) *
                              df2_int$Y_death) /
                         mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                                df2_int$M_diabetes == 1))) -
  (mean(wt * as.numeric(df2_int$A0_PM2.5 == 0 &
                          df2_int$M_diabetes == 1) *
          df2_int$Y_death) /
     mean(wt * as.numeric(df2_int$A0_PM2.5 == 0 &
                            df2_int$M_diabetes == 1)))
CDE_sIPTW_m1_death
# [1] 0.09030186

CDE_sIPTW_m1_qol <- (mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                            df2_int$M_diabetes == 1) *
                            df2_int$Y_qol) /
                       mean(wt * as.numeric(df2_int$A0_PM2.5 == 1 &
                                              df2_int$M_diabetes == 1))) -
  (mean(wt * as.numeric(df2_int$A0_PM2.5 == 0 &
                          df2_int$M_diabetes == 1) *
          df2_int$Y_qol) /
     mean(wt * as.numeric(df2_int$A0_PM2.5 == 0 &
                            df2_int$M_diabetes == 1)))
CDE_sIPTW_m1_qol
# [1] -10.03045



# ---------------------------------------------------------------------------- #
# 2) Estimation of the Conditional Randomized Direct & Indirect Effects ----
# ---------------------------------------------------------------------------- #
# from Zheng & van der Laan 2017

rm(list=ls())
df2_int <- read.csv(file = "./data/df2_int.csv")

## In order to estimate Psi^{a,a'}

## D^{a,a'} = Y * I(A=a) / p_A(a | l(0)) * [p_M(M | a',L(1),L(0)) / p_A(a | l(0)) * [p_M(M | a',L(1),L(0)) ]

## 1. Estimate gA and gM
gA <- glm(A0_PM2.5 ~ L0_male + L0_soc_env,
          family = "binomial", data = df2_int)
gM <- glm(M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
          family = "binomial", data = df2_int)

## 2. Predict each individual's probability of being exposed to her own exposure
# predict the probabilities P(A0_PM2.5=1|L(0)) & P(A0_PM2.5=0|L(0))
pred.gA1 <- predict(gA, type = "response")
pred.gA0 <- 1 - pred.gA1
# the predicted probability of the observed treatment A_i=a is :
gAobs <- rep(NA, nrow(df2_int))
gAobs[df2_int$A0_PM2.5 == 1] <- pred.gA1[df2_int$A0_PM2.5 == 1]
gAobs[df2_int$A0_PM2.5 == 0] <- pred.gA0[df2_int$A0_PM2.5 == 0]

# predict the probabilities
# P(M=1|L(0),A=1,L(1)) & P(M=0|L(0),A=1,L(1)) setting A = 1
# and P(M=1|L(0),A=0,L(1)) & P(M=0|L(0),A=0,L(1)) setting A = 0
data.Ais0 <- data.Ais1 <- df2_int
data.Ais0$A0_PM2.5 <- 0
data.Ais1$A0_PM2.5 <- 1

pred.gM1.Ais1 <- predict(gM, newdata = data.Ais1, type = "response")
pred.gM0.Ais1 <- 1 - pred.gM1.Ais1
# the predicted probability of the observed treatment M_i=m is :
gMobs.Ais1 <- rep(NA, nrow(df2_int))
gMobs.Ais1[df2_int$M_diabetes == 1] <- pred.gM1.Ais1[df2_int$M_diabetes == 1]
gMobs.Ais1[df2_int$M_diabetes == 0] <- pred.gM0.Ais1[df2_int$M_diabetes == 0]

pred.gM1.Ais0 <- predict(gM, newdata = data.Ais0, type = "response")
pred.gM0.Ais0 <- 1 - pred.gM1.Ais0
# the predicted probability of the observed treatment M_i=m is :
gMobs.Ais0 <- rep(NA, nrow(df2_int))
gMobs.Ais0[df2_int$M_diabetes == 1] <- pred.gM1.Ais0[df2_int$M_diabetes == 1]
gMobs.Ais0[df2_int$M_diabetes == 0] <- pred.gM0.Ais0[df2_int$M_diabetes == 0]

## 3. Calculate D^{a,a'}
D.death.11 <- (df2_int$Y_death * (I(df2_int$A0_PM2.5 == 1) / gAobs) *
                 (gMobs.Ais1 / gMobs.Ais1))
D.death.10 <- (df2_int$Y_death * (I(df2_int$A0_PM2.5 == 1) / gAobs) *
                 (gMobs.Ais0 / gMobs.Ais1))
D.death.00 <- (df2_int$Y_death * (I(df2_int$A0_PM2.5 == 0) / gAobs) *
                 (gMobs.Ais0 / gMobs.Ais0))

D.qol.11 <- (df2_int$Y_qol * (I(df2_int$A0_PM2.5 == 1) / gAobs) *
               (gMobs.Ais1 / gMobs.Ais1))
D.qol.10 <- (df2_int$Y_qol * (I(df2_int$A0_PM2.5 == 1) / gAobs) *
               (gMobs.Ais0 / gMobs.Ais1))
D.qol.00 <- (df2_int$Y_qol * (I(df2_int$A0_PM2.5 == 0) / gAobs) *
               (gMobs.Ais0 / gMobs.Ais0))

## 4. Calculate CRDE and CRIE
## For death
CRDE.death <- mean(D.death.10) - mean(D.death.00)
CRDE.death
# [1] 0.07405068
CRIE.death <- mean(D.death.11) - mean(D.death.10)
CRIE.death
# [1] 0.00819879

## For quality of life
CRDE.qol <- mean(D.qol.10) - mean(D.qol.00)
CRDE.qol
# [1] -7.563116
CRIE.qol <- mean(D.qol.11) - mean(D.qol.10)
CRIE.qol
# [1] -0.8736813

## 5. Calculate 95% CI based on the influence curve D^{a,a'}
# the variance of the estimator Psi^{a,a'} = mean(D^{a,a'}) is var(D^{a,a'}) / n)
N <- nrow(df2_int)
se.CRDE.death <- sqrt(var(D.death.10 - D.death.00) / N)
c(CRDE.death - qnorm(0.975) * se.CRDE.death,
  CRDE.death + qnorm(0.975) * se.CRDE.death)
# [1] 0.04141535 0.10668602
se.CRIE.death <- sqrt(var(D.death.11 - D.death.10) / N)
c(CRIE.death - qnorm(0.975) * se.CRIE.death,
  CRIE.death + qnorm(0.975) * se.CRIE.death)
# [1] 0.003380238 0.013017342

se.CRDE.qol <- sqrt(var(D.qol.10 - D.qol.00) / N)
c(CRDE.qol - qnorm(0.975) * se.CRDE.qol,
  CRDE.qol + qnorm(0.975) * se.CRDE.qol)
# [1] -11.748095  -3.378137
se.CRIE.qol <- sqrt(var(D.qol.11 - D.qol.10) / N)
c(CRIE.qol - qnorm(0.975) * se.CRIE.qol,
  CRIE.qol + qnorm(0.975) * se.CRIE.qol)
# [1] -1.4111404 -0.3362223


## check on simulations
# ---------------------------------------------------------------------------- #
### test on simulations for df2 ----
# ---------------------------------------------------------------------------- #
rm(list=ls())

# functions to generate data
param.causal.model.2 <- function(A.M.interaction = NULL) {
  # L0
  p_L0_male <- 0.5
  p_L0_soc_env <- 0.65

  # A: A0_PM2.5 <- rbinom( 0.05 + 0.04 * L0_male + 0.06 * L0_soc_env )
  b_A <- 0.05   # reference prevalence is 5%
  b_male_A <- 0.04  # + 0.04 for the effect of L0_male -> A0_PM2.5
  b_soc_env_A <- 0.06  # +0.06 for the effect of L0_soc_env -> A0_PM2.5

  # L1: L1 <- rbinom( 0.30 - 0.05 * L0_male + 0.08 * L0_soc_env +
  #                   0.2 * A0_PM2.5 )
  b_L1 <- 0.30   # reference prevalence is 30%
  b_male_L1 <- -0.05  # - 0.05 for the effect of L0_male -> L1
  b_soc_env_L1 <- +0.08 # + 0.08 for the effect of L0_soc_env -> L1
  b_A_L1 <- +0.2 # +0.2 for the effect of A0_PM2.5 -> L1

  # M: M_diabetes <- rbinom( 0.2 + 0.05 * L0_male + 0.06 * L0_soc_env +
  #                         0.2 * L1 + 0.1 * A0_PM2.5 )
  b_M <- 0.2 # reference prevalence is 20%
  b_male_M <- 0.05 # +0.05 for the effect of L0_male -> M_diabetes
  b_soc_env_M <- 0.06 # +0.06 for the effect of L0_soc_env -> M_diabetes
  b_A_M <- 0.1 # +0.10 for the effect of A0_PM2.5 -> M_diabetes
  b_L1_M <- 0.2 # +0.2 for the effect of L1 -> M_diabetes

  # Y binary: rbinom( 0.10 + 0.06 * L0_male + 0.04 * L0_soc_env +
  #                   0.05 * A0_PM2.5 + 0.07 * L1 + 0.08 * M_diabetes +
  #                   0.03 * A0_PM2.5 * M_diabetes * A.M.inter )
  b_Y <- 0.1 # reference prevalence is 10%
  b_male_Y <- 0.06 # +0.06 for the effect of L0_male -> Y
  b_soc_env_Y <- 0.04 # +0.04 for the effect of L0_soc_env -> Y
  b_A_Y <- 0.05 # 0.05 for the effect of A0_PM2.5 -> Y
  b_L1_Y <- 0.07 # +0.07 for the effect of L1 -> Y
  b_M_Y <- 0.08 # 0.08 for the effect of M_diabetes -> Y
  b_AM_Y <- 0.03 # 0.03 for the interaction effect A0_PM2.5 * M_diabetes -> Y

  # Y continuous: (75 - 1 * L0_male - 3 * L0_soc_env - 4 * A0_PM2.5 +
  #                -3.5 * L1 - 9 * M_diabetes +
  #                -5 * A0_PM2.5 * M_diabetes * A.M.inter ) + rnorm(N, mean = 0, sd = 10)
  mu_Y <- 75 # reference mean for QoL
  c_male_Y <- -1 # -1 for the effect of L0_male -> Y
  c_soc_env_Y <- -3 # -3 for the effect of L0_soc_env -> Y
  c_A_Y <- -4 # -4 for the effect of A0_PM2.5 -> Y
  c_L1_Y <- -5 # -5 for the effect of L1 -> Y
  c_M_Y <- -9 # -9 for the effect of M_diabetes -> Y
  c_AM_Y <- -5  # - 5 for the interaction effect A0_PM2.5 * M_diabetes  -> Y
  sd_Y <- 10 # standard deviation of the residuals

  # A*M interaction ?
  A.M.inter <- A.M.interaction

  coef <- c( p_L0_male = p_L0_male, p_L0_soc_env = p_L0_soc_env,
             b_A = b_A, b_male_A = b_male_A, b_soc_env_A = b_soc_env_A,
             b_L1 = b_L1, b_male_L1 = b_male_L1, b_soc_env_L1 = b_soc_env_L1,
             b_A_L1 = b_A_L1,
             b_M = b_M, b_male_M = b_male_M, b_soc_env_M = b_soc_env_M,
             b_L1_M = b_L1_M, b_A_M = b_A_M,
             b_Y = b_Y, b_male_Y = b_male_Y, b_soc_env_Y = b_soc_env_Y,
             b_A_Y = b_A_Y, b_L1_Y = b_L1_Y, b_M_Y = b_M_Y, b_AM_Y = b_AM_Y,
             mu_Y = mu_Y, c_male_Y = c_male_Y, c_soc_env_Y = c_soc_env_Y,
             c_A_Y = c_A_Y, c_L1_Y = c_L1_Y, c_M_Y = c_M_Y, c_AM_Y = c_AM_Y,
             sd_Y = sd_Y, A.M.inter = A.M.inter)

  return(coef)
}

gen.data.causal.model.2 <- function(N, A.M.inter) { # input parameters are the
  #   sample size N and the presence of A*M interaction with A.M.inter = 0 or 1

  b <- param.causal.model.2(A.M.interaction = A.M.inter)

  # baseline confounders: parent's educational level=L0_soc_env & sex=L0_male
  L0_male <- rbinom(N, size = 1, prob = b["p_L0_male"])
  L0_soc_env <- rbinom(N, size = 1, prob = b["p_L0_soc_env"])

  # exposure: A0_PM2.5
  A0_PM2.5 <- rbinom(N, size = 1, prob =  b["b_A"] +
                       b["b_male_A"] * L0_male +
                       b["b_soc_env_A"] * L0_soc_env )

  # intermediate confounder between M_diabetes and Y,
  L1 <- rbinom(N, size = 1, prob = b["b_L1"] +
                 b["b_male_L1"] * L0_male +
                 b["b_soc_env_L1"] * L0_soc_env +
                 b["b_A_L1"]* A0_PM2.5)

  # mediator: M_diabetes
  M_diabetes <- rbinom(N, size = 1, prob = b["b_M"] +
                         b["b_male_M"] * L0_male +
                         b["b_soc_env_M"] * L0_soc_env +
                         b["b_A_M"] * A0_PM2.5 +
                         b["b_L1_M"] * L1)

  # Y_death
  Y_death <- rbinom(N, size = 1, prob = b["b_Y"] +
                      b["b_male_Y"] * L0_male +
                      b["b_soc_env_Y"] * L0_soc_env +
                      b["b_A_Y"] * A0_PM2.5 +
                      b["b_L1_Y"] * L1 +
                      b["b_M_Y"] * M_diabetes +
                      b["b_AM_Y"] * A0_PM2.5 * M_diabetes * A.M.inter )

  # Y_qol
  Y_qol <- ( b["mu_Y"] +
               b["c_male_Y"] * L0_male +
               b["c_soc_env_Y"] * L0_soc_env +
               b["c_A_Y"] * A0_PM2.5 +
               b["c_L1_Y"] * L1 +
               b["c_M_Y"] * M_diabetes +
               b["c_AM_Y"] * A0_PM2.5 * M_diabetes * A.M.inter ) +
    rnorm(N, mean = 0, sd = b["sd_Y"])

  # data.frame
  data.sim <- data.frame(L0_male, L0_soc_env, A0_PM2.5, L1, M_diabetes,
                         Y_death, Y_qol)

  return( data.sim )
}

# create matrix to save simulation results
results.df2 <- matrix(NA, nrow = 1000, ncol = 16,
                      dimnames = list(c(1:1000),
                                      c("CRDE.death.gcomp", "CRIE.death.gcomp","total.death.gcomp",
                                        "CRDE.qol.gcomp", "CRIE.qol.gcomp","total.qol.gcomp",
                                        "CRDE.death.IPTW", "CRIE.death.IPTW","total.death.IPTW",
                                        "CRDE.qol.IPTW", "CRIE.qol.IPTW","total.qol.IPTW",
                                        "CRDE.death.cov", "CRIE.death.cov",
                                        "CRDE.qol.cov", "CRIE.qol.cov")))
set.seed(54321)
for(i in 1:1000) {
  print(paste0("simulation ",i))
  # generate data
  df2_int <- gen.data.causal.model.2(N=10000, A.M.inter=1)

  ### g-comp algorithm for CRDE and CRIE # # # # # # # # # # # # # # # # # # # #
  ## 1) Q^{a,a'}_R2 =  Y
  Q.death_R2 <- df2_int$Y_death
  Q.qol_R2 <- df2_int$Y_qol

  ## 2) Obtain Q^{a,a'}_L1, Q^{a,a'}_M1, Q^{a,a'}_R1
  ## 2.a) Obtain Q^{a,a'}_L1
  # Regress Q^{a,a'}_R2 on observed values (L(0),A,L(1),M)
  L1.model.death <- glm(Q.death_R2 ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                          M_diabetes + A0_PM2.5:M_diabetes,
                        family = "binomial", data = df2_int)
  L1.model.qol <- glm(Q.qol_R2 ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                        M_diabetes + A0_PM2.5:M_diabetes,
                      family = "gaussian", data = df2_int)

  # Evaluate the fitted function at the observed mediator M and covariate history L(1),L(0)
  # and the intervened exposure A = a
  data.Ais0 <- data.Ais1 <- df2_int
  data.Ais0$A0_PM2.5 <- 0
  data.Ais1$A0_PM2.5 <- 1

  # We will need 3 counterfactual quantities Q^{a,a'}: Q^{1,1}, Q^{1,0} and Q^{0,0}
  Q11.death.L1 <- predict(L1.model.death, newdata = data.Ais1, type = "response")
  Q10.death.L1 <- predict(L1.model.death, newdata = data.Ais1, type = "response")
  Q00.death.L1 <- predict(L1.model.death, newdata = data.Ais0, type = "response")

  Q11.qol.L1 <- predict(L1.model.qol, newdata = data.Ais1, type = "response")
  Q10.qol.L1 <- predict(L1.model.qol, newdata = data.Ais1, type = "response")
  Q00.qol.L1 <- predict(L1.model.qol, newdata = data.Ais0, type = "response")

  ## 2.b) Obtain Q^{a,a'}_M1
  # Regress Q^{a,a'}_L1 on observed values (L(0),A,L(1))
  M1.model.death.11 <- glm(Q11.death.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                           family = "quasibinomial", data = df2_int)
  M1.model.death.10 <- glm(Q10.death.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                           family = "quasibinomial", data = df2_int)
  M1.model.death.00 <- glm(Q00.death.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                           family = "quasibinomial", data = df2_int)

  M1.model.qol.11 <- glm(Q11.qol.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                         family = "gaussian", data = df2_int)
  M1.model.qol.10 <- glm(Q10.qol.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                         family = "gaussian", data = df2_int)
  M1.model.qol.00 <- glm(Q00.qol.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                         family = "gaussian", data = df2_int)

  # Evaluate the fitted function at the observed covariate history L(1),L(0)
  # and the intervened exposure A = a'
  Q11.death.M1 <- predict(M1.model.death.11, newdata = data.Ais1, type = "response")
  Q10.death.M1 <- predict(M1.model.death.10, newdata = data.Ais0, type = "response")
  Q00.death.M1 <- predict(M1.model.death.00, newdata = data.Ais0, type = "response")

  Q11.qol.M1 <- predict(M1.model.qol.11, newdata = data.Ais1, type = "response")
  Q10.qol.M1 <- predict(M1.model.qol.10, newdata = data.Ais0, type = "response")
  Q00.qol.M1 <- predict(M1.model.qol.00, newdata = data.Ais0, type = "response")

  ## 2.c) Obtain Q^{a,a'}_R1
  # Regress Q^{a,a'}_M1 on observed values (L(0),A)
  R1.model.death.11 <- glm(Q11.death.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)
  R1.model.death.10 <- glm(Q10.death.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)
  R1.model.death.00 <- glm(Q00.death.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)

  R1.model.qol.11 <- glm(Q11.qol.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "gaussian", data = df2_int)
  R1.model.qol.10 <- glm(Q10.qol.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "gaussian", data = df2_int)
  R1.model.qol.00 <- glm(Q00.qol.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "gaussian", data = df2_int)

  # Evaluate the fitted function at the observed covariate history L(0)
  # and the intervened exposure A = a
  Q11.death.R1 <- predict(R1.model.death.11, newdata = data.Ais1, type = "response")
  Q10.death.R1 <- predict(R1.model.death.10, newdata = data.Ais1, type = "response")
  Q00.death.R1 <- predict(R1.model.death.00, newdata = data.Ais0, type = "response")

  Q11.qol.R1 <- predict(R1.model.qol.11, newdata = data.Ais1, type = "response")
  Q10.qol.R1 <- predict(R1.model.qol.10, newdata = data.Ais1, type = "response")
  Q00.qol.R1 <- predict(R1.model.qol.00, newdata = data.Ais0, type = "response")

  ## 4) Estimate the Conditional Randomized Direct and Indirect Effects
  ## For death
  CRDE.death.gcomp <- mean(Q10.death.R1) - mean(Q00.death.R1)
  CRIE.death.gcomp <- mean(Q11.death.R1) - mean(Q10.death.R1)
  ## For quality of life
  CRDE.qol.gcomp <- mean(Q10.qol.R1) - mean(Q00.qol.R1)
  CRIE.qol.gcomp <- mean(Q11.qol.R1) - mean(Q10.qol.R1)

  ### save results
  results.df2[i,"CRDE.death.gcomp"] <- CRDE.death.gcomp - (0.078282) # direct 0.078282
  results.df2[i,"CRIE.death.gcomp"] <- CRIE.death.gcomp - (0.011) # indirect 0.011
  results.df2[i,"total.death.gcomp"] <- CRDE.death.gcomp + CRIE.death.gcomp - (0.089282) # direct 0.089282
  results.df2[i,"CRDE.qol.gcomp"] <- CRDE.qol.gcomp - (-7.207) # indirect -7.207
  results.df2[i,"CRIE.qol.gcomp"] <- CRIE.qol.gcomp - (-1.4) # ATE -1.4
  results.df2[i,"total.qol.gcomp"] <- CRDE.qol.gcomp + CRIE.qol.gcomp - (-8.607) # direct -8.607


  ### IPTW algorithm for CRDE and CRIE # # # # # # # # # # # # # # # # # # # #
  ## 1. Estimate gA and gM
  gA <- glm(A0_PM2.5 ~ L0_male + L0_soc_env,
            family = "binomial", data = df2_int)
  gM <- glm(M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
            family = "binomial", data = df2_int)

  ## 2. Predict each individual's probability of being exposed to her own exposure
  # predict the probabilities P(A0_PM2.5=1|L(0)) & P(A0_PM2.5=0|L(0))
  pred.gA1 <- predict(gA, type = "response")
  pred.gA0 <- 1 - pred.gA1
  # the predicted probability of the observed treatment A_i=a is :
  gAobs <- rep(NA, nrow(df2_int))
  gAobs[df2_int$A0_PM2.5 == 1] <- pred.gA1[df2_int$A0_PM2.5 == 1]
  gAobs[df2_int$A0_PM2.5 == 0] <- pred.gA0[df2_int$A0_PM2.5 == 0]

  # predict the probabilities
  # P(M=1|L(0),A=1,L(1)) & P(M=0|L(0),A=1,L(1)) setting A = 1
  # and P(M=1|L(0),A=0,L(1)) & P(M=0|L(0),A=0,L(1)) setting A = 0
  # data.Ais0 <- data.Ais1 <- df2_int
  # data.Ais0$A0_PM2.5 <- 0
  # data.Ais1$A0_PM2.5 <- 1

  pred.gM1.Ais1 <- predict(gM, newdata = data.Ais1, type = "response")
  pred.gM0.Ais1 <- 1 - pred.gM1.Ais1
  # the predicted probability of the observed treatment M_i=m is :
  gMobs.Ais1 <- rep(NA, nrow(df2_int))
  gMobs.Ais1[df2_int$M_diabetes == 1] <- pred.gM1.Ais1[df2_int$M_diabetes == 1]
  gMobs.Ais1[df2_int$M_diabetes == 0] <- pred.gM0.Ais1[df2_int$M_diabetes == 0]

  pred.gM1.Ais0 <- predict(gM, newdata = data.Ais0, type = "response")
  pred.gM0.Ais0 <- 1 - pred.gM1.Ais0
  # the predicted probability of the observed treatment M_i=m is :
  gMobs.Ais0 <- rep(NA, nrow(df2_int))
  gMobs.Ais0[df2_int$M_diabetes == 1] <- pred.gM1.Ais0[df2_int$M_diabetes == 1]
  gMobs.Ais0[df2_int$M_diabetes == 0] <- pred.gM0.Ais0[df2_int$M_diabetes == 0]

  ## 3. Calculate D^{a,a'}
  D.death.11 <- (df2_int$Y_death * (I(df2_int$A0_PM2.5 == 1) / gAobs) *
                   (gMobs.Ais1 / gMobs.Ais1))
  D.death.10 <- (df2_int$Y_death * (I(df2_int$A0_PM2.5 == 1) / gAobs) *
                   (gMobs.Ais0 / gMobs.Ais1))
  D.death.00 <- (df2_int$Y_death * (I(df2_int$A0_PM2.5 == 0) / gAobs) *
                   (gMobs.Ais0 / gMobs.Ais0))

  D.qol.11 <- (df2_int$Y_qol * (I(df2_int$A0_PM2.5 == 1) / gAobs) *
                 (gMobs.Ais1 / gMobs.Ais1))
  D.qol.10 <- (df2_int$Y_qol * (I(df2_int$A0_PM2.5 == 1) / gAobs) *
                 (gMobs.Ais0 / gMobs.Ais1))
  D.qol.00 <- (df2_int$Y_qol * (I(df2_int$A0_PM2.5 == 0) / gAobs) *
                 (gMobs.Ais0 / gMobs.Ais0))

  ## 4. Calculate CRDE and CRIE
  ## For death
  CRDE.death.IPTW <- mean(D.death.10) - mean(D.death.00)
  CRIE.death.IPTW <- mean(D.death.11) - mean(D.death.10)

  ## For quality of life
  CRDE.qol.IPTW <- mean(D.qol.10) - mean(D.qol.00)
  CRIE.qol.IPTW <- mean(D.qol.11) - mean(D.qol.10)

  ## 5. Calculate 95% CI based on the influence curve D^{a,a'}
  # the variance of the estimator Psi^{a,a'} = mean(D^{a,a'}) is var(D^{a,a'}) / n)
  N <- nrow(df2_int)
  se.CRDE.death <- sqrt(var(D.death.10 - D.death.00) / N)
  c(CRDE.death.IPTW - qnorm(0.975) * se.CRDE.death,
    CRDE.death.IPTW + qnorm(0.975) * se.CRDE.death)
  # [1] 0.04141535 0.10668602
  se.CRIE.death <- sqrt(var(D.death.11 - D.death.10) / N)
  c(CRIE.death.IPTW - qnorm(0.975) * se.CRIE.death,
    CRIE.death.IPTW + qnorm(0.975) * se.CRIE.death)
  # [1] 0.003380238 0.013017342

  se.CRDE.qol <- sqrt(var(D.qol.10 - D.qol.00) / N)
  c(CRDE.qol.IPTW - qnorm(0.975) * se.CRDE.qol,
    CRDE.qol.IPTW + qnorm(0.975) * se.CRDE.qol)
  # [1] -11.748095  -3.378137
  se.CRIE.qol <- sqrt(var(D.qol.11 - D.qol.10) / N)
  c(CRIE.qol.IPTW - qnorm(0.975) * se.CRIE.qol,
    CRIE.qol.IPTW + qnorm(0.975) * se.CRIE.qol)
  # [1] -1.4111404 -0.3362223

  ### save results
  results.df2[i,"CRDE.death.IPTW"] <- CRDE.death.IPTW - (0.078282) # direct 0.078282
  results.df2[i,"CRIE.death.IPTW"] <- CRIE.death.IPTW - (0.011) # indirect 0.011
  results.df2[i,"total.death.IPTW"] <- CRDE.death.IPTW + CRIE.death.IPTW - (0.089282) # direct 0.089282
  results.df2[i,"CRDE.qol.IPTW"] <- CRDE.qol.IPTW - (-7.207) # indirect -7.207
  results.df2[i,"CRIE.qol.IPTW"] <- CRIE.qol.IPTW - (-1.4) # ATE -1.4
  results.df2[i,"total.qol.IPTW"] <- CRDE.qol.IPTW + CRIE.qol.IPTW - (-8.607) # direct -8.607

  results.df2[i,"CRDE.death.cov"] <- ((0.078282 > CRDE.death.IPTW -
                                        qnorm(0.975) * se.CRDE.death) &
                                        (0.078282 < CRDE.death.IPTW +
                                           qnorm(0.975) * se.CRDE.death))
  results.df2[i,"CRIE.death.cov"] <- ((0.011 > CRIE.death.IPTW -
                                        qnorm(0.975) * se.CRIE.death) &
                                        (0.011 < CRIE.death.IPTW +
                                           qnorm(0.975) * se.CRIE.death))
  results.df2[i,"CRDE.qol.cov"] <- ((-7.207 > CRDE.qol.IPTW -
                                      qnorm(0.975) * se.CRDE.qol) &
                                      (-7.207 < CRDE.qol.IPTW +
                                         qnorm(0.975) * se.CRDE.qol))
  results.df2[i,"CRIE.qol.cov"] <- ((-1.4 > CRIE.qol.IPTW -
                                      qnorm(0.975) * se.CRIE.qol) &
                                      (-1.4 < CRIE.qol.IPTW +
                                         qnorm(0.975) * se.CRIE.qol))

}

sapply(data.frame(results.df2) , mean, na.rm = TRUE)
# CRDE.death.gcomp  CRIE.death.gcomp total.death.gcomp    CRDE.qol.gcomp    CRIE.qol.gcomp   total.qol.gcomp
#    -2.351832e-03     -6.237371e-04     -2.975569e-03     -1.949122e-03     -1.147990e-03     -3.097112e-03
# CRDE.death.IPTW   CRIE.death.IPTW  total.death.IPTW     CRDE.qol.IPTW     CRIE.qol.IPTW    total.qol.IPTW
#   -4.510467e-04     -4.272167e-05     -4.937684e-04     -2.741841e-01      4.604712e-02     -2.281370e-01
# CRDE.death.cov    CRIE.death.cov      CRDE.qol.cov      CRIE.qol.cov
#   9.870000e-01      9.470000e-01      1.000000e+00      9.780000e-01 => 95% CI are too wide but this is expected with IPTW

boxplot(results.df2[,1:12])
boxplot(subset(data.frame(results.df2), select = c("CRDE.death.gcomp", "CRDE.death.IPTW"))) # IPTW slightly better
boxplot(subset(data.frame(results.df2), select = c("CRIE.death.gcomp", "CRIE.death.IPTW"))) # IPTW slightly better
boxplot(subset(data.frame(results.df2), select = c("total.death.gcomp", "total.death.IPTW"))) # IPTW slightly better

boxplot(subset(data.frame(results.df2), select = c("CRDE.qol.gcomp", "CRDE.qol.IPTW"))) # gcomp sligthly better
boxplot(subset(data.frame(results.df2), select = c("CRIE.qol.gcomp", "CRIE.qol.IPTW"))) # gcomp sligthly better
boxplot(subset(data.frame(results.df2), select = c("total.qol.gcomp", "total.qol.IPTW"))) # gcomp sligthly better


results.df2.rel <- results.df2
results.df2.rel[,c("CRDE.death.gcomp","CRDE.death.IPTW")] <- results.df2.rel[,c("CRDE.death.gcomp","CRDE.death.IPTW")] / (0.078282)
results.df2.rel[,c("CRIE.death.gcomp","CRIE.death.IPTW")] <- results.df2.rel[,c("CRIE.death.gcomp","CRIE.death.IPTW")] / (0.011)
results.df2.rel[,c("total.death.gcomp","total.death.IPTW")] <- results.df2.rel[,c("total.death.gcomp","total.death.IPTW")] / (0.089282)
results.df2.rel[,c("CRDE.qol.gcomp","CRDE.qol.IPTW")] <- results.df2.rel[,c("CRDE.qol.gcomp","CRDE.qol.IPTW")] / (-7.207)
results.df2.rel[,c("CRIE.qol.gcomp","CRIE.qol.IPTW")] <- results.df2.rel[,c("CRIE.qol.gcomp","CRIE.qol.IPTW")] / (-1.4)
results.df2.rel[,c("total.qol.gcomp","total.qol.IPTW")] <- results.df2.rel[,c("total.qol.gcomp","total.qol.IPTW")] / (-8.607)
sapply(data.frame(results.df2.rel) , mean)
# CRDE.death.gcomp  CRIE.death.gcomp total.death.gcomp    CRDE.qol.gcomp    CRIE.qol.gcomp   total.qol.gcomp
#    -0.0300430786     -0.0567033737     -0.0333277636      0.0002704484      0.0008199927      0.0003598364
# CRDE.death.IPTW   CRIE.death.IPTW  total.death.IPTW     CRDE.qol.IPTW     CRIE.qol.IPTW    total.qol.IPTW
#   -0.0057618190     -0.0038837880     -0.0055304360      0.0380441388     -0.0328907991      0.0265059823
# CRDE.death.cov    CRIE.death.cov      CRDE.qol.cov      CRIE.qol.cov
#   0.9870000000      0.9470000000      1.0000000000      0.9780000000

# OK
