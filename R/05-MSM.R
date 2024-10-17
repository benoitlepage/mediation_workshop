### test programs for estimations by MSMs

rm(list=ls())

df1 <- read.csv(file = "data/df1.csv")
df1_int <- read.csv(file = "data/df1_int.csv")
df2 <- read.csv(file = "data/df2.csv")
df2_int <- read.csv(file = "data/df2_int.csv")

# ---------------------------------------------------------------------------- #
# 1) Estimation of the Average Total Effect (ATE) ----
# ---------------------------------------------------------------------------- #

rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")
### MSM of ATE|L_male, estimated by IPTW ---------------------------------------

## 1. Denominator of the weight
# 1a. Estimate g(A=a_i|L(0)) (denominator of the weight)
g.A.L <- glm(A0_PM2.5 ~ L0_male + L0_soc_env,
           family = "binomial", data = df1_int)

# 1b. Predict each individual's probability of being exposed to her own exposure
# predict the probabilities P(A0_PM2.5=1) & P(A0_PM2.5=0)
pred.g1.L <- predict(g.A.L, type="response")
pred.g0.L <- 1 - pred.g1.L
# the predicted probability of the observed treatment P(A = a_i | L(0)) is :
gAi.L <- rep(NA, nrow(df1_int))
gAi.L[df1_int$A0_PM2.5==1] <- pred.g1.L[df1_int$A0_PM2.5==1]
gAi.L[df1_int$A0_PM2.5==0] <- pred.g0.L[df1_int$A0_PM2.5==0]

## 2. Numerator of the weight
# The numerator of the weight can be 1 for simple weights,
# or g(A=a_i|V) to obtain stabilized weights which put less weight to individuals
# with less observation. Stabilized weights enable a weaker positivity assumption.

# 2a. Estimate g(A=a_i | sex) (numerator of the stabilized weight)
g.A.sex <- glm(A0_PM2.5 ~ L0_male,
               family = "binomial", data = df1_int)

# 2b. Predict each individual's probability of being exposed to her own exposure
# predict the probabilities P(A0_PM2.5=1 | sex) & P(A0_PM2.5=0 | sex)
pred.g1.sex <- predict(g.A.sex, type="response")
pred.g0.sex <- 1 - pred.g1.sex
# the predicted probability of the observed treatment P(A = a_i | sex) is :
gAi.sex <- rep(NA, nrow(df1_int))
gAi.sex[df1_int$A0_PM2.5==1] <- pred.g1.sex[df1_int$A0_PM2.5==1]
gAi.sex[df1_int$A0_PM2.5==0] <- pred.g0.sex[df1_int$A0_PM2.5==0]


## 3. Define individual weights:
# We can use simple weights w = 1 / g(A=a_i | L(0))
w <- 1 / gAi.L

# Or alternatively, we can use stabilized weights : sw = g(A=a_i) / g(A=a_i | L(0))
sw <- gAi.sex / gAi.L

par(mfcol = c(1,2))
boxplot(w ~ df1_int$A0_PM2.5)
boxplot(sw ~ df1_int$A0_PM2.5)
par(mfcol = c(1,1))

## applying these weights creates a pseudo-population were the baseline
## confounders are balanced, relative to the exposure:
## before applying weights to the individuals:
table(df1_int$L0_male, df1_int$A0_PM2.5, deparse.level = 2)
#                df1_int$A0_PM2.5
# df1_int$L0_male    0    1
#               0 4527  463
#               1 4349  661
prop.table(table(df1_int$L0_male, df1_int$A0_PM2.5, deparse.level = 2),
           margin = 2)
#                          df1_int$A0_PM2.5
# df1_int$L0_male         0         1
#               0 0.5100270 0.4119217
#               1 0.4899730 0.5880783

## after applying weights to the individuals:
library(questionr) # The questionr package enables to describe weighted populations
wtd.table(x = df1_int$L0_male, y = df1_int$A0_PM2.5,
          weights = w)
#          0        1
# 0 4989.425 4918.462
# 1 5010.862 5057.676
prop.table(wtd.table(x = df1_int$L0_male, y = df1_int$A0_PM2.5,
                     weights = w), margin = 2)
#           0         1
# 0 0.4989282 0.4930227
# 1 0.5010718 0.5069773


## 4. Estimate coefficients of the MSM using a weighted regression E(Y | A, sex)
# a GLM with gaussian family can be applied to estimate risk difference
msm1 <- glm(Y_death ~ A0_PM2.5 + L0_male + A0_PM2.5*L0_male,
           weights = w,
           family = "gaussian",
           data = df1_int)
coef(msm1)
# (Intercept)         A0_PM2.5          L0_male A0_PM2.5:L0_male
#  0.17573472       0.05883228       0.04598911       0.03077614

msm2 <- glm(Y_death ~ A0_PM2.5 + L0_male + A0_PM2.5*L0_male,
            weights = sw,
            family = "gaussian",
            data = df1_int)
coef(msm2)
# (Intercept)         A0_PM2.5          L0_male A0_PM2.5:L0_male
#  0.17573472       0.05883228       0.04598911       0.03077614

## 5. Estimate the ATE stratified by sex
# According to MSM1 (with simple weights)
ATE.msm1.male0 <- coef(msm1)["A0_PM2.5"]
# 0.05883228
ATE.msm1.male1 <- coef(msm1)["A0_PM2.5"] + coef(msm1)["A0_PM2.5:L0_male"]
# 0.08960842

# According to MSM2 (with stabilized weights)
ATE.msm2.male0 <- coef(msm2)["A0_PM2.5"]
# 0.05883228
ATE.msm2.male1 <- coef(msm2)["A0_PM2.5"] + coef(msm2)["A0_PM2.5:L0_male"]
# 0.08960842
# results are the same because there is no violation of the positivity assumption

# en pratique, il faudrait vérifier s'il y a une différence entre homme et femme
# à partir du système générateur de données ?

### Checking only - not shown in the document
## 6. Estimation of 95% confidence interval by bootstrap
set.seed(1234)
B <- 2000
bootstrap.estimates <- data.frame(matrix(NA, nrow = B, ncol = 4))
colnames(bootstrap.estimates) <- c("boot.msm1.male0", "boot.msm1.male1",
                                   "boot.msm2.male0", "boot.msm2.male1")
for (b in 1:B){
  # sample the indices 1 to n with replacement
  bootIndices <- sample(1:nrow(df1_int), replace=T)
  bootData <- df1_int[bootIndices,]

  if ( round(b/100, 0) == b/100 ) print(paste0("bootstrap number ",b))

  # 1a. Estimate g(A=a_i|L(0)) (denominator of the weight)
  g.A.L <- glm(A0_PM2.5 ~ L0_male + L0_soc_env,
               family = "binomial", data = bootData)

  # 1b. Predict each individual's probability of being exposed to her own exposure
  # the predicted probability of the observed treatment P(A = a_i | L(0)) is :
  gAi.L <- rep(NA, nrow(bootData))
  gAi.L[bootData$A0_PM2.5==1] <- predict(g.A.L, type="response")[bootData$A0_PM2.5==1]
  gAi.L[bootData$A0_PM2.5==0] <- (1 - predict(g.A.L, type="response"))[bootData$A0_PM2.5==0]

  # 2a. Estimate g(A=a_i | sex) (numerator of the stabilized weight)
  g.A.sex <- glm(A0_PM2.5 ~ L0_male,
                 family = "binomial", data = bootData)

  # 2b. Predict each individual's probability of being exposed to her own exposure
  # the predicted probability of the observed treatment P(A = a_i | sex) is :
  gAi.sex <- rep(NA, nrow(bootData))
  gAi.sex[bootData$A0_PM2.5==1] <- predict(g.A.sex, type="response")[bootData$A0_PM2.5==1]
  gAi.sex[bootData$A0_PM2.5==0] <- (1 - predict(g.A.sex, type="response"))[bootData$A0_PM2.5==0]

  # 3. Define individual weights:
  w <- 1 / gAi.L
  sw <- gAi.sex / gAi.L

  # 4. Estimate coefficients of the MSM using a weighted regression E(Y | A, sex)
  msm1 <- glm(Y_death ~ A0_PM2.5 + L0_male + A0_PM2.5*L0_male,
              weights = w,
              family = "gaussian",
              data = bootData)

  msm2 <- glm(Y_death ~ A0_PM2.5 + L0_male + A0_PM2.5*L0_male,
              weights = sw,
              family = "gaussian",
              data = bootData)

  ## 5. Estimate the ATE stratified by sex
  # Using coefficients of msm1
  bootstrap.estimates[b,"boot.msm1.male0"] <- coef(msm1)["A0_PM2.5"]
  bootstrap.estimates[b,"boot.msm1.male1"] <- (coef(msm1)["A0_PM2.5"] +
                                                 coef(msm1)["A0_PM2.5:L0_male"])

  # Using coefficients of msm2
  bootstrap.estimates[b,"boot.msm2.male0"] <- coef(msm2)["A0_PM2.5"]
  bootstrap.estimates[b,"boot.msm2.male1"] <- (coef(msm2)["A0_PM2.5"] +
                                                 coef(msm2)["A0_PM2.5:L0_male"])
}

IC95.ATE.msm1.male0 <- c(ATE.msm1.male0 -
                           qnorm(0.975)*sd(bootstrap.estimates[,"boot.msm1.male0"]),
                         ATE.msm1.male0 +
                           qnorm(0.975)*sd(bootstrap.estimates[,"boot.msm1.male0"]))
# 0.01878219 0.09888237

IC95.ATE.msm1.male1 <- c(ATE.msm1.male1 -
                           qnorm(0.975)*sd(bootstrap.estimates[,"boot.msm1.male1"]),
                         ATE.msm1.male1 +
                           qnorm(0.975)*sd(bootstrap.estimates[,"boot.msm1.male1"]))
# 0.05152628 0.12769057

IC95.ATE.msm2.male0 <- c(ATE.msm2.male0 -
                           qnorm(0.975)*sd(bootstrap.estimates[,"boot.msm2.male0"]),
                         ATE.msm2.male0 +
                           qnorm(0.975)*sd(bootstrap.estimates[,"boot.msm2.male0"]))
# 0.01878219 0.09888237

IC95.ATE.msm2.male1 <- c(ATE.msm2.male1 -
                           qnorm(0.975)*sd(bootstrap.estimates[,"boot.msm2.male1"]),
                         ATE.msm2.male1 +
                           qnorm(0.975)*sd(bootstrap.estimates[,"boot.msm2.male1"]))
# 0.05152628 0.12769057 # ne contient pas la vraie valeur ??

# encore une fois pas de différence entre poids simple et poids stabilisé car
# pas de problème de positivité ?
par(mfcol = c(1,1))
hist(bootstrap.estimates$boot.msm1.male0)
hist(bootstrap.estimates$boot.msm1.male1)
hist(bootstrap.estimates$boot.msm2.male0)
hist(bootstrap.estimates$boot.msm2.male1)

# Theoretical ATE in women and men:
param.causal.model.1 <- function(A.M.interaction = NULL) {
  # L0
  p_L0_male <- 0.5
  p_L0_soc_env <- 0.65

  # A: A0_PM2.5 <- rbinom( 0.05 + 0.04 * L0_male + 0.06 * L0_soc_env )
  b_A <- 0.05   # reference prevalence is 5%
  b_male_A <- 0.04  # + 0.04 for the effect of L0_male -> A0_PM2.5
  b_parent_educ_A <- 0.06  # +0.06 for the effect of L0_soc_env -> A0_PM2.5

  # L1: intermediate confounder between M and Y, not influenced by A
  p_L1 <- 0.3

  # M: M_diabetes <- rbinom( 0.2 + 0.05 * L0_male + 0.06 * L0_soc_env + 0.07 * L1 +
  #                         0.1 * A0_PM2.5 )
  b_M <- 0.2 # reference prevalence is 20%
  b_male_M <- 0.05 # +0.05 for the effect of L0_male -> M_diabetes
  b_parent_educ_M <- 0.06 # +0.06 for the effect of L0_soc_env -> M_diabetes
  b_L1_M <- 0.07 # +0.07 for the effect of L1 -> M_diabetes
  b_A_M <- 0.1 # +0.10 for the effect of A0_PM2.5 -> M_diabetes

  # Y binary: rbinom( 0.10 + 0.06 * L0_male + 0.04 * L0_soc_env + 0.05 * A0_PM2.5 +
  #                   0.07 * L1 + 0.08 * M_diabetes +
  #                   0.03 * A0_PM2.5 * M_diabetes * A.M.inter )
  b_Y <- 0.1 # reference prevalence is 10%
  b_male_Y <- 0.06 # +0.06 for the effect of L0_male -> Y
  b_parent_educ_Y <- 0.04 # +0.04 for the effect of L0_soc_env -> Y
  b_A_Y <- 0.05 # 0.05 for the effect of A0_PM2.5 -> Y
  b_L1_Y <- 0.07 # +0.07 for the effect of L1 -> Y
  b_M_Y <- 0.08 # 0.08 for the effect of M_diabetes -> Y
  b_AM_Y <- 0.03 # 0.03 for the interaction effect A0_PM2.5 * M_diabetes -> Y

  # Y continuous: (75 - 1 * L0_male - 3 * L0_soc_env - 4 * A0_PM2.5 -3.5 * L1 -
  #                9 * M_diabetes -5 * A0_PM2.5 * M_diabetes * A.M.inter ) +
  #                rnorm(N, mean = 0, sd = 10)
  mu_Y <- 75 # reference mean for QoL
  c_male_Y <- -1 # -1 for the effect of L0_male -> Y
  c_parent_educ_Y <- -3 # -3 for the effect of L0_soc_env -> Y
  c_A_Y <- -4 # -4 for the effect of A0_PM2.5 -> Y
  c_L1_Y <- -3.5 # -3.5 for the effect of L1 -> Y
  c_M_Y <- -9 # -9 for the effect of M_diabetes -> Y
  c_AM_Y <- -5  # - 5 for the interaction effect A0_PM2.5 * M_diabetes  -> Y
  sd_Y <- 10 # standard deviation of the residuals

  # A*M interaction ?
  A.M.inter <- A.M.interaction

  coef <- c( p_L0_male = p_L0_male, p_L0_soc_env = p_L0_soc_env,
             b_A = b_A, b_male_A = b_male_A, b_parent_educ_A = b_parent_educ_A,
             p_L1 = p_L1,
             b_M = b_M, b_male_M = b_male_M, b_parent_educ_M = b_parent_educ_M,
             b_L1_M = b_L1_M, b_A_M = b_A_M,
             b_Y = b_Y, b_male_Y = b_male_Y, b_parent_educ_Y = b_parent_educ_Y,
             b_A_Y = b_A_Y, b_L1_Y = b_L1_Y, b_M_Y = b_M_Y, b_AM_Y = b_AM_Y,
             mu_Y = mu_Y, c_male_Y = c_male_Y, c_parent_educ_Y = c_parent_educ_Y,
             c_A_Y = c_A_Y, c_L1_Y = c_L1_Y, c_M_Y = c_M_Y, c_AM_Y = c_AM_Y,
             sd_Y = sd_Y, A.M.inter = A.M.inter)

  return(coef)
}

true.ATE1 <- function(interaction = NULL) {
  b <- param.causal.model.1(A.M.interaction = interaction)

  # binary outcome (death)
  S1 <- S0 <- cbind(expand.grid(c(0,1),c(0,1),c(0,1)),rep(NA,n=2^3))
  colnames(S1) <- colnames(S0) <- list("parent_educ","L1","M","sum")
  for (n in 1:8) {
    S1[n,"sum"] <- ( ( ( b["b_Y"] +
                          b["b_male_Y"] * 1 +
                          b["b_parent_educ_Y"] * S1[n,"parent_educ"] +
                          b["b_A_Y"] * 1 +
                          b["b_L1_Y"] * S1[n,"L1"] +
                          b["b_M_Y"] * S1[n,"M"] +
                          b["b_AM_Y"] * 1 * S1[n,"M"] * b["A.M.inter"] ) *
                        (( b["b_M"] +
                             b["b_male_M"] * 1 +
                             b["b_parent_educ_M"] * S1[n,"parent_educ"] +
                             b["b_L1_M"] * S1[n,"L1"] +
                             b["b_A_M"] * 1 )^( S1[n,"M"] )) *
                        (( 1 - (b["b_M"] +
                                  b["b_male_M"] * 1 +
                                  b["b_parent_educ_M"] * S1[n,"parent_educ"] +
                                  b["b_L1_M"] * S1[n,"L1"] +
                                  b["b_A_M"] * 1) )^( 1 - S1[n,"M"] )) ) -
                      ( ( b["b_Y"] +
                            b["b_male_Y"] * 1 +
                            b["b_parent_educ_Y"] * S1[n,"parent_educ"] +
                            b["b_A_Y"] * 0 +
                            b["b_L1_Y"] * S1[n,"L1"] +
                            b["b_M_Y"] * S1[n,"M"] +
                            b["b_AM_Y"] * 0 * S1[n,"M"] * b["A.M.inter"] ) *
                          (( b["b_M"] +
                               b["b_male_M"] * 1 +
                               b["b_parent_educ_M"] * S1[n,"parent_educ"] +
                               b["b_L1_M"] * S1[n,"L1"] +
                               b["b_A_M"] * 0 )^( S1[n,"M"] )) *
                          (( 1 - (b["b_M"] +
                                    b["b_male_M"] * 1 +
                                    b["b_parent_educ_M"] * S1[n,"parent_educ"] +
                                    b["b_L1_M"] * S1[n,"L1"] +
                                    b["b_A_M"] * 0) )^( 1 - S1[n,"M"] )) ) ) *
      ((b["p_L1"])^(S1[n,"L1"])) *
      ((1 - b["p_L1"])^(1 - S1[n,"L1"])) *
      # ((b["p_L0_male"])^1) *
      # ((1 - b["p_L0_male"])^(1 - 1)) *
      ((b["p_L0_soc_env"])^(S1[n,"parent_educ"])) *
      ((1 - b["p_L0_soc_env"])^(1 - S1[n,"parent_educ"]))

    S0[n,"sum"] <- ( ( ( b["b_Y"] +
                           b["b_male_Y"] * 0 +
                           b["b_parent_educ_Y"] * S0[n,"parent_educ"] +
                           b["b_A_Y"] * 1 +
                           b["b_L1_Y"] * S0[n,"L1"] +
                           b["b_M_Y"] * S0[n,"M"] +
                           b["b_AM_Y"] * 1 * S0[n,"M"] * b["A.M.inter"] ) *
                         (( b["b_M"] +
                              b["b_male_M"] * 0 +
                              b["b_parent_educ_M"] * S0[n,"parent_educ"] +
                              b["b_L1_M"] * S0[n,"L1"] +
                              b["b_A_M"] * 1 )^( S0[n,"M"] )) *
                         (( 1 - (b["b_M"] +
                                   b["b_male_M"] * 0 +
                                   b["b_parent_educ_M"] * S0[n,"parent_educ"] +
                                   b["b_L1_M"] * S0[n,"L1"] +
                                   b["b_A_M"] * 1) )^( 1 - S0[n,"M"] )) ) -
                       ( ( b["b_Y"] +
                             b["b_male_Y"] * 0 +
                             b["b_parent_educ_Y"] * S0[n,"parent_educ"] +
                             b["b_A_Y"] * 0 +
                             b["b_L1_Y"] * S0[n,"L1"] +
                             b["b_M_Y"] * S0[n,"M"] +
                             b["b_AM_Y"] * 0 * S0[n,"M"] * b["A.M.inter"] ) *
                           (( b["b_M"] +
                                b["b_male_M"] * 0 +
                                b["b_parent_educ_M"] * S0[n,"parent_educ"] +
                                b["b_L1_M"] * S0[n,"L1"] +
                                b["b_A_M"] * 0 )^( S0[n,"M"] )) *
                           (( 1 - (b["b_M"] +
                                     b["b_male_M"] * 0 +
                                     b["b_parent_educ_M"] * S0[n,"parent_educ"] +
                                     b["b_L1_M"] * S0[n,"L1"] +
                                     b["b_A_M"] * 0) )^( 1 - S0[n,"M"] )) ) ) *
      ((b["p_L1"])^(S0[n,"L1"])) *
      ((1 - b["p_L1"])^(1 - S0[n,"L1"])) *
      # ((b["p_L0_male"])^(1 - 1)) *
      # ((1 - b["p_L0_male"])^(1)) *
      ((b["p_L0_soc_env"])^(S0[n,"parent_educ"])) *
      ((1 - b["p_L0_soc_env"])^(1 - S0[n,"parent_educ"]))
  }

  ATE.death.male0 <- sum(S0[,"sum"])
  ATE.death.male1 <- sum(S1[,"sum"])

  return(list(ATE.death.male0 = ATE.death.male0, ATE.death.male1 = ATE.death.male1))
}

true <- true.ATE1(interaction = 1)
true
# $ATE.death.male0
# [1] 0.0688
#
# $ATE.death.male1
# [1] 0.0703

# il y a une petite différence, mais faible !
# l'effet total est la moyenne des deux : 0.06955 => oui, ça a l'air d'être ça !

true2 <- true.ATE1(interaction = 0)
true2
# la petite différence précédente venait du terme d'interaction :
# sans terme d'interaction, les deux effets sont les même qq soit le sexe => OK


rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")
### MSM of ATE|L_male, estimated by G-computation ------------------------------
## 1. Estimate Qbar
Q.tot.death <- glm(Y_death ~ A0_PM2.5 + L0_male + L0_soc_env,
                   family = "gaussian", data = df1_int)
# The final result would be sligthly different if we applied a binomial family
# The Gaussian family corresponds to the true generating model in this example.

## 2. Predict an outcome for each subject, setting A=0 and A=1
# prepare data sets used to predict the outcome under the counterfactual
# scenarios setting A=0 and A=1
data.A1 <- data.A0 <- df1_int
data.A1$A0_PM2.5 <- 1
data.A0$A0_PM2.5 <- 0

# predict values under the same name in the corresponding counterfactual dataset
data.A1$Ya.death.pred <- predict(Q.tot.death, newdata = data.A1, type = "response")
data.A0$Ya.death.pred <- predict(Q.tot.death, newdata = data.A0, type = "response")

## 3. Append both counterfactual datasets in a single dataset
# number of row is twice the initial value (we have 2 counterfactual scenarios)
data.2scenarios <- rbind(data.A0, data.A1)

## 4. fit the MSM: E(Y_a|sex)
MSM.ATE.gcomp <- glm(Ya.death.pred ~ A0_PM2.5 + L0_male + A0_PM2.5:L0_male,
                     family = "gaussian",
                     data = data.2scenarios)
coef(MSM.ATE.gcomp)
#  (Intercept)         A0_PM2.5          L0_male A0_PM2.5:L0_male
# 1.743994e-01     7.720726e-02     4.874750e-02     5.087110e-16

## 5. Estimate the ATE stratified by sex
# According to MSM.ATE.gcomp
ATE.MSM.gcomp.male0 <- coef(MSM.ATE.gcomp)["A0_PM2.5"]
# 0.07720726
ATE.MSM.gcomp.male1 <- (coef(MSM.ATE.gcomp)["A0_PM2.5"] +
                          coef(MSM.ATE.gcomp)["A0_PM2.5:L0_male"])
# 0.07720726
# The results are the same in both strata, because in the first Qbar model,
# we did not include any (A * sex) interaction term

# Applying a binomial family for the first Qbar model would result in two
# different values of the ATE stratified by sex.
# => 0.06880798 in the L0_male = 0 strata
# => 0.08053575 in the L0_male = 1 strata

# Indeed, applying a gaussian family (additive model) with no interaction terms
# implies the assumption of some interaction terms in a multiplicative model (such
# as a glm with a binomial family).
# On the contrary, applying a binomial family (multiplicative model) with no
# interaction terms implies some interaction terms in an additive model (such as
# a glm with a gaussian family)

# Using the true data generating model used to simulate the illustrative datasets,
# the "true" value of the ATE stratified by sex can be calculated:
# (ATE | L0_male = 0) = 0.0688
# (ATE | L0_male = 1) = 0.0703


# ---------------------------------------------------------------------------- #
# 2) Estimation of the CDE ----
# ---------------------------------------------------------------------------- #
## estimation by IPTW ----
### without intermediate confounder affected by the exposure ----
rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")
### MSM of CDE, estimated by IPTW
## 1. Stabilized weight for the exposure sw_{A,i}
# 1a. Estimate g(A=a_i|L(0)) (denominator of the weight)
g.A.L <- glm(A0_PM2.5 ~ L0_male + L0_soc_env,
             family = "binomial", data = df1_int)
# 1b. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i | L(0)) is :
gAi.L <- rep(NA, nrow(df1_int))
gAi.L[df1_int$A0_PM2.5==1] <- predict(g.A.L, type="response")[df1_int$A0_PM2.5==1]
gAi.L[df1_int$A0_PM2.5==0] <- (1 - predict(g.A.L, type="response"))[df1_int$A0_PM2.5==0]

# 1c. Estimate g(A=a_i) (numerator of the weight)
g.A <- glm(A0_PM2.5 ~ 1, family = "binomial", data = df1_int)
# 1d. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i) is :
gAi <- rep(NA, nrow(df1_int))
gAi[df1_int$A0_PM2.5==1] <- predict(g.A, type="response")[df1_int$A0_PM2.5==1]
gAi[df1_int$A0_PM2.5==0] <- (1 - predict(g.A, type="response"))[df1_int$A0_PM2.5==0]

# 1e. Calculate sw_{A,i}
sw_Ai <- gAi / gAi.L

## 2. Stabilized weight for the mediator sw_{M,i}
# 2a. Estimate g(M=m_i|L(0),A,L(1)) (denominator of the weight)
g.M.L <- glm(M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
             family = "binomial", data = df1_int)
# 2b. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i | L(0)) is :
gMi.L <- rep(NA, nrow(df1_int))
gMi.L[df1_int$M_diabetes==1] <- predict(g.M.L, type="response")[df1_int$M_diabetes==1]
gMi.L[df1_int$M_diabetes==0] <- (1 - predict(g.M.L, type="response"))[df1_int$M_diabetes==0]

# 2c. Estimate g(M=m_i|A) (numerator of the weight)
g.M.A <- glm(M_diabetes ~ A0_PM2.5, family = "binomial", data = df1_int)
# 2d. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(M = m_i|A) is :
gMi.A <- rep(NA, nrow(df1_int))
gMi.A[df1_int$M_diabetes==1] <- predict(g.M.A, type="response")[df1_int$M_diabetes==1]
gMi.A[df1_int$M_diabetes==0] <- (1 - predict(g.M.A, type="response"))[df1_int$M_diabetes==0]
# 2e. Calculate sw_{M,i}
sw_Mi <- gMi.A / gMi.L

## 3. Define the individual stabilized weight for the CDE_m
sw_cde <- sw_Ai * sw_Mi

## 4. Estimate coefficients of the MSM using a weighted regression E(Y | A, sex)
# a GLM with gaussian family can be applied to estimate risk difference
msm_cde <- glm(Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5*M_diabetes,
               weights = sw_cde,
               family = "gaussian",
               data = df1_int)
coef(msm_cde)
# (Intercept)            A0_PM2.5          M_diabetes A0_PM2.5:M_diabetes
#  0.17891689          0.06798282          0.06729724         -0.00495314

## 5. Estimate CDE for m=0 and for m=1 using the MSM's coefficients
CDE_mis0 <- coef(msm_cde)["A0_PM2.5"]
# 0.06798282
CDE_mis1 <- coef(msm_cde)["A0_PM2.5"] + coef(msm_cde)["A0_PM2.5:M_diabetes"]
# 0.06302968


### with intermediate confounder affected by the exposure ----
### MSM of CDE, estimated by IPTW
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")

## 1. Stabilized weight for the exposure sw_{A,i}
# 1a. Estimate g(A=a_i|L(0)) (denominator of the weight)
g.A.L <- glm(A0_PM2.5 ~ L0_male + L0_soc_env,
             family = "binomial", data = df2_int)
# 1b. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i | L(0)) is :
gAi.L <- rep(NA, nrow(df2_int))
gAi.L[df2_int$A0_PM2.5==1] <- predict(g.A.L, type="response")[df2_int$A0_PM2.5==1]
gAi.L[df2_int$A0_PM2.5==0] <- (1 - predict(g.A.L, type="response"))[df2_int$A0_PM2.5==0]

# 1c. Estimate g(A=a_i) (numerator of the weight)
g.A <- glm(A0_PM2.5 ~ 1, family = "binomial", data = df2_int)
# 1d. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i) is :
gAi <- rep(NA, nrow(df2_int))
gAi[df2_int$A0_PM2.5==1] <- predict(g.A, type="response")[df2_int$A0_PM2.5==1]
gAi[df2_int$A0_PM2.5==0] <- (1 - predict(g.A, type="response"))[df2_int$A0_PM2.5==0]

# 1e. Calculate the weight for the exposure A: sw_{A,i}
sw_Ai <- gAi / gAi.L

## 2. Stabilized weight for the mediator sw_{M,i}
# 2a. Estimate g(M=m_i|L(0),A,L(1)) (denominator of the weight)
g.M.L <- glm(M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
             family = "binomial", data = df2_int)
# 2b. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i | L(0)) is :
gMi.L <- rep(NA, nrow(df2_int))
gMi.L[df2_int$M_diabetes==1] <- predict(g.M.L, type="response")[df2_int$M_diabetes==1]
gMi.L[df2_int$M_diabetes==0] <- (1 - predict(g.M.L, type="response"))[df2_int$M_diabetes==0]

# 2c. Estimate g(M=m_i|A) (numerator of the weight)
g.M.A <- glm(M_diabetes ~ A0_PM2.5, family = "binomial", data = df2_int)
# 2d. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(M = m_i|A) is :
gMi.A <- rep(NA, nrow(df2_int))
gMi.A[df2_int$M_diabetes==1] <- predict(g.M.A, type="response")[df2_int$M_diabetes==1]
gMi.A[df2_int$M_diabetes==0] <- (1 - predict(g.M.A, type="response"))[df2_int$M_diabetes==0]
# 2e. Calculate the weight for the mediator M: sw_{M,i}
sw_Mi <- gMi.A / gMi.L

## 3. Define the individual stabilized weight for the CDE_m
sw_cde <- sw_Ai * sw_Mi

## 4. Estimate coefficients of the MSM using a weighted regression E(Y | A, sex)
# a GLM with gaussian family can be applied to estimate risk differences
msm_cde <- glm(Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5*M_diabetes,
               weights = sw_cde,
               family = "gaussian",
               data = df2_int)
coef(msm_cde)
# (Intercept)            A0_PM2.5          M_diabetes A0_PM2.5:M_diabetes
#  0.17932146          0.06012920          0.07239661          0.03017266

## 5. Estimate CDE for m=0 and for m=1 using the MSM's coefficients
CDE_mis0 <- coef(msm_cde)["A0_PM2.5"]
# 0.0601292
CDE_mis1 <- coef(msm_cde)["A0_PM2.5"] + coef(msm_cde)["A0_PM2.5:M_diabetes"]
# 0.09030186


### example with ltmle
library(ltmle)
Qform <- c(L1="Q.kplus1 ~ L0_male + L0_soc_env + A0_PM2.5",
           Y_death="Q.kplus1 ~ L0_male + L0_soc_env + L1 +
                    A0_PM2.5 * M_diabetes")
gform <- c("A0_PM2.5 ~ L0_male + L0_soc_env",
           "M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5 + L1")
data_binary <- subset(df2_int, select = c(L0_male, L0_soc_env,
                                          A0_PM2.5, L1,
                                          M_diabetes, Y_death))
CDE_ltmle_M0_death <- ltmle(data = data_binary,
                            Anodes = c("A0_PM2.5", "M_diabetes"),
                            Lnodes = c("L1"), # intermediate confounders +/- baseline
                            Ynodes = c("Y_death"),
                            survivalOutcome = FALSE, # TRUE for time-to-event outcomes Y
                            Qform = Qform,
                            gform = gform,
                            abar = list(c(1,0), # counterfactual intervention do(A=1,M=0)
                                        c(0,0)), # counterfactual intervention do(A=0,M=0)
                            SL.library = NULL,
                            estimate.time = FALSE, # estimate computation time?
                            iptw.only = TRUE, # for only IPTW estimations (no TMLE)
                            gcomp = FALSE,
                            variance.method = "iptw") # IPTW influence curve
# CDE with M=0
summary(CDE_ltmle_M0_death, estimator = "iptw")$effect.measures$ATE
# $estimate
# [1] 0.0601292
# $CI
#            2.5%      97.5%
# [1,] 0.02424819 0.09601021

CDE_ltmle_M1_death <- ltmle(data = data_binary,
                            Anodes = c("A0_PM2.5", "M_diabetes"),
                            Lnodes = c("L1"), # intermediate confounders +/- baseline
                            Ynodes = c("Y_death"),
                            survivalOutcome = FALSE, # TRUE for time-to-event outcomes Y
                            Qform = Qform,
                            gform = gform,
                            abar = list(c(1,1), # counterfactual intervention do(A=1,M=0)
                                        c(0,1)), # counterfactual intervention do(A=0,M=0)
                            SL.library = NULL,
                            estimate.time = FALSE, # estimate computation time?
                            iptw.only = TRUE, # for only IPTW estimations (no TMLE)
                            gcomp = FALSE,
                            variance.method = "iptw") # IPTW influence curve
# CDE with M=1
summary(CDE_ltmle_M1_death, estimator = "iptw")$effect.measures$ATE
# $estimate
# [1] 0.09030186
# $CI
#            2.5%     97.5%
# [1,] 0.04084792 0.1397558

## we obtain the same results as the IPTW estimation by hand

library(CMAverse)
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
cmdag(outcome = "Y_death", exposure = "A0_PM2.5", mediator = "M_diabetes",
      basec = c("L0_male", "L0_soc_env"), postc = "L1", node = TRUE, text_col = "white")
# In this setting, we can use the marginal structural model and the $g$-formula approach. The results are shown below.

## The Marginal Structural Model
res_msm_RD <- cmest(data = df2_int,
                    model = "msm",
                    outcome = "Y_death",
                    exposure = "A0_PM2.5",
                    mediator = "M_diabetes",
                    basec = c("L0_male", "L0_soc_env"),
                    postc = "L1",
                    EMint = TRUE, # E*M interaction
                    ereg = "logistic",
                    yreg = "linear", # MSM is a linear regression (to get RD)
                    mreg = list("logistic"),
                    wmnomreg = list("logistic"),
                    wmdenomreg = list("logistic"),
                    astar = 0, #E(Y_{A=0,M=1})
                    a = 1,  #E(Y_{A=1,M=1})
                    mval = list(1), # for the CDE, set mediator to M=1
                    estimation = "imputation",
                    inference = "bootstrap",
                    nboot = 2)

summary(res_msm_RD)
# Causal Mediation Analysis
# # Outcome regression:
# Call:
#   glm(formula = Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5 * M_diabetes,
#       family = gaussian(), data = getCall(x$reg.output$yreg)$data,
#       weights = getCall(x$reg.output$yreg)$weights)
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)      # MSM coef
#   (Intercept)         0.179321   0.005242  34.210  < 2e-16 ***  # are the
#   A0_PM2.5            0.060129   0.017253   3.485 0.000494 ***  # same as
#   M_diabetes          0.072397   0.009242   7.834 5.21e-15 ***  # calculation
#   A0_PM2.5:M_diabetes 0.030173   0.025895   1.165 0.243960      # by hand
#
# # Mediator regressions:
# Call:
#   glm(formula = M_diabetes ~ A0_PM2.5, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -0.73432    0.02268 -32.384  < 2e-16 ***
#   A0_PM2.5     0.51531    0.06422   8.024 1.02e-15 ***
#
# # Mediator regressions for weighting (denominator):
# Call:
#   glm(formula = M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env +
#         L1, family = binomial(), data = getCall(x$reg.output$wmdenomreg[[1L]])$data,
#       weights = getCall(x$reg.output$wmdenomreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -1.36249    0.04783 -28.488  < 2e-16 ***
#   A0_PM2.5     0.30994    0.06668   4.648 3.35e-06 ***
#   L0_male      0.24661    0.04369   5.644 1.66e-08 ***
#   L0_soc_env   0.30628    0.04650   6.587 4.50e-11 ***
#   L1           0.86045    0.04493  19.152  < 2e-16 ***
#
# # Mediator regressions for weighting (nominator):
# Call:
#   glm(formula = M_diabetes ~ A0_PM2.5, family = binomial(), data = getCall(x$reg.output$wmnomreg[[1L]])$data,
#       weights = getCall(x$reg.output$wmnomreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -0.74205    0.02271 -32.680   <2e-16 ***
#   A0_PM2.5     0.55288    0.06408   8.628   <2e-16 ***
#
# # Exposure regression for weighting:
# Call:
#   glm(formula = A0_PM2.5 ~ L0_male + L0_soc_env, family = binomial(),
#       data = getCall(x$reg.output$ereg)$data, weights = getCall(x$reg.output$ereg)$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -2.73244    0.07425 -36.799  < 2e-16 ***
#   L0_male      0.40580    0.06447   6.294 3.09e-10 ***
#   L0_soc_env   0.64060    0.07350   8.716  < 2e-16 ***
#
# # Effect decomposition on the mean difference scale via the marginal structural model
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#                   Estimate  Std.error    95% CIL 95% CIU  P.val
#   cde            9.030e-02  1.601e-02  6.094e-02   0.082 <2e-16 ***
#   rpnde          7.003e-02  1.600e-02  5.838e-02   0.080 <2e-16 ***
#   rtnde          7.374e-02  1.597e-02  5.888e-02   0.080 <2e-16 ***
#   rpnie          8.903e-03  9.159e-04  8.284e-03   0.010 <2e-16 ***
#   rtnie          1.261e-02  9.523e-04  8.738e-03   0.010 <2e-16 ***
#   te             8.265e-02  1.505e-02  6.839e-02   0.089 <2e-16 ***
#   rintref       -2.027e-02  7.341e-06 -2.571e-03  -0.003 <2e-16 ***
#   rintmed        3.710e-03  3.634e-05  4.542e-04   0.001 <2e-16 ***
#   cde(prop)      1.093e+00  2.940e-02  8.907e-01   0.930 <2e-16 ***
#   rintref(prop) -2.452e-01  6.289e-03 -3.752e-02  -0.029 <2e-16 ***
#   rintmed(prop)  4.489e-02  1.662e-03  5.140e-03   0.007 <2e-16 ***
#   rpnie(prop)    1.077e-01  3.403e-02  9.376e-02   0.139 <2e-16 ***
#   rpm            1.526e-01  3.569e-02  9.890e-02   0.147 <2e-16 ***
#   rint          -2.004e-01  4.627e-03 -3.014e-02  -0.024 <2e-16 ***
#   rpe           -9.263e-02  2.940e-02  6.983e-02   0.109 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#  (cde: controlled direct effect;
#  rpnde: randomized analogue of pure natural direct effect;
#  rtnde: randomized analogue of total natural direct effect;
#  rpnie: randomized analogue of pure natural indirect effect;
#  rtnie: randomized analogue of total natural indirect effect;
#  te: total effect; rintref: randomized analogue of reference interaction;
#  rintmed: randomized analogue of mediated interaction;
#  cde(prop): proportion cde;
#  rintref(prop): proportion rintref;
#  rintmed(prop): proportion rintmed;
#  rpnie(prop): proportion rpnie;
#  rpm: randomized analogue of overall proportion mediated;
#  rint: randomized analogue of overall proportion attributable to interaction;
#  rpe: randomized analogue of overall proportion eliminated)


## MSM of CDE by G-computation
### MSM of CDE, estimated by G-computation (df1_int) ----
rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")
## 1. Estimate Qbar(A,M,L0,L1)
Q.cde.death <- glm(Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5:M_diabetes + L0_male + L0_soc_env + L1,
                   family = "gaussian", data = df1_int)
# The final result would be sligthly different if we applied a binomial family
# The Gaussian family corresponds to the true generating model in this example.

## 2. Predict an outcome for each subject, in each counterfactual scenario
# prepare data sets used to predict the outcome under the counterfactual
# 4 counterfactual scenarios setting (A=0,M=0), (A=1,M=0), (A=0,M=1) and (A=1,M=1)
data.A0M0 <- data.A1M0 <- data.A0M1 <- data.A1M1 <- df1_int
data.A0M0$A0_PM2.5 <- 0
data.A0M0$M_diabetes <- 0

data.A1M0$A0_PM2.5 <- 1
data.A1M0$M_diabetes <- 0

data.A0M1$A0_PM2.5 <- 0
data.A0M1$M_diabetes <- 1

data.A1M1$A0_PM2.5 <- 1
data.A1M1$M_diabetes <- 1

# predict values under the same name in the corresponding counterfactual dataset
data.A0M0$Yam.death.pred <- predict(Q.cde.death, newdata = data.A0M0, type = "response")
data.A1M0$Yam.death.pred <- predict(Q.cde.death, newdata = data.A1M0, type = "response")
data.A0M1$Yam.death.pred <- predict(Q.cde.death, newdata = data.A0M1, type = "response")
data.A1M1$Yam.death.pred <- predict(Q.cde.death, newdata = data.A1M1, type = "response")

## 3. Append both counterfactual datasets in a single dataset
# number of row is twice the initial value (we have 2 counterfactual scenarios)
data.4scenarios <- rbind(data.A0M0, data.A1M0,data.A0M1,data.A1M1)

## 4. fit the MSM: E(Y_a|sex)
MSM.CDE.gcomp <- glm(Yam.death.pred ~ A0_PM2.5 +  M_diabetes + A0_PM2.5:M_diabetes,
                     family = "gaussian", # gaussian family for risk differences
                     data = data.4scenarios)
coef(MSM.CDE.gcomp)
# (Intercept)            A0_PM2.5          M_diabetes A0_PM2.5:M_diabetes
#  0.17968603          0.06000138          0.06757214          0.01918153

## 5. Estimate the CDEm
# CDE(M=0) = E(Y_{A=1,M=0}) - E(Y_{A=0,M=0})
CDE_mis0_gcomp <- coef(MSM.CDE.gcomp)["A0_PM2.5"]
# 0.06000138

# CDE(M=1) = E(Y_{A=1,M=1}) - E(Y_{A=0,M=1})
CDE_mis1_gcomp <- (coef(MSM.CDE.gcomp)["A0_PM2.5"] +
                     coef(MSM.CDE.gcomp)["A0_PM2.5:M_diabetes"])
# 0.07918291

# Note: Applying a binomial family for the first Qbar model would result in two
# sligthly different values of the CDE(M=m)
# => 0.05934409 in setting M=0
# => 0.07537874 in setting M=1


rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
### MSM of CDE, estimated by G-computation (by ICE) (df2_int) ----
## 1a) Regress the outcome on L0, A, L1 and M (and the A*M interaction if appropriate)
Y.death.model <- glm(Y_death ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                       M_diabetes + A0_PM2.5:M_diabetes,
                     family = "binomial", data = df2_int)

## 1b) Generate predicted values by evaluating the regression setting the mediator
##    value to M=0 or to M=1
data.A0M0 <- data.A1M0 <- data.A0M1 <- data.A1M1 <- df2_int
data.A0M0$A0_PM2.5 <- 0
data.A0M0$M_diabetes <- 0

data.A1M0$A0_PM2.5 <- 1
data.A1M0$M_diabetes <- 0

data.A0M1$A0_PM2.5 <- 0
data.A0M1$M_diabetes <- 1

data.A1M1$A0_PM2.5 <- 1
data.A1M1$M_diabetes <- 1

Q.Y.death.A0M0 <- predict(Y.death.model, newdata = data.A0M0, type = "response")
Q.Y.death.A1M0 <- predict(Y.death.model, newdata = data.A1M0, type = "response")
Q.Y.death.A0M1 <- predict(Y.death.model, newdata = data.A0M1, type = "response")
Q.Y.death.A1M1 <- predict(Y.death.model, newdata = data.A1M1, type = "response")

## 2a) Regress the predicted values conditional on the observed exposure A
##    and baseline confounders L(0)
L1.death.A0M0.model <- glm(Q.Y.death.A0M0 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)
L1.death.A1M0.model <- glm(Q.Y.death.A1M0 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)
L1.death.A0M1.model <- glm(Q.Y.death.A0M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)
L1.death.A1M1.model <- glm(Q.Y.death.A1M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)


## 2b) generate predicted values by evaluating the regression at exposure
##    of interest: {A=0,M=0}, {A=1,M=0}, {A=0,M=1}, {A=1,M=1}
data.A0M0$Yam.death.pred <- predict(L1.death.A0M0.model,
                                    newdata = data.A0M0, type = "response")
data.A1M0$Yam.death.pred <- predict(L1.death.A1M0.model,
                                    newdata = data.A1M0, type = "response")
data.A0M1$Yam.death.pred <- predict(L1.death.A0M1.model,
                                    newdata = data.A0M1, type = "response")
data.A1M1$Yam.death.pred <- predict(L1.death.A1M1.model,
                                    newdata = data.A1M1, type = "response")


## 3. Append the 4 counterfactual datasets in a single long dataset
# number of row is 4 times the initial value (we have 4 counterfactual scenarios)
data.4scenarios <- rbind(data.A0M0, data.A1M0,data.A0M1,data.A1M1)

## 4. fit the MSM: E(Y_am) = alpha_0 + alpha_A a + alpha_M m + alpha_AM a:m
MSM.CDE.gcomp <- glm(Yam.death.pred ~ A0_PM2.5 +  M_diabetes + A0_PM2.5:M_diabetes,
                     family = "gaussian", # gaussian family for risk differences
                     data = data.4scenarios)
coef(MSM.CDE.gcomp)
# (Intercept)            A0_PM2.5          M_diabetes A0_PM2.5:M_diabetes
#  0.17974947          0.06342833          0.07366466          0.02469485

## 5. Estimate the CDE(M=m)
# CDE(M=0) = E(Y_{A=1,M=0}) - E(Y_{A=0,M=0})
CDE_mis0_gcomp_ice <- coef(MSM.CDE.gcomp)["A0_PM2.5"]
# 0.06342833

# CDE(M=1) = E(Y_{A=1,M=1}) - E(Y_{A=0,M=1})
CDE_mis1_gcomp_ice <- (coef(MSM.CDE.gcomp)["A0_PM2.5"] +
                         coef(MSM.CDE.gcomp)["A0_PM2.5:M_diabetes"])
# 0.08812318


##### Alternative --------------------------------------------------------------
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
### MSM of CDE, estimated by G-computation (by ICE) ----------------------------
## 1a) Regress the outcome on L0, A, L1 and M (and the A*M interaction if appropriate)
Y.death.model <- glm(Y_death ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                       M_diabetes + A0_PM2.5:M_diabetes,
                     family = "binomial", data = df2_int)

## 1b) Generate predicted values by evaluating the regression setting the mediator
##    value to M=0 or to M=1
data.A0M0 <- data.A1M0 <- data.A0M1 <- data.A1M1 <- df2_int
# data.A0M0$A0_PM2.5 <- 0
data.A0M0$M_diabetes <- 0

# data.A1M0$A0_PM2.5 <- 1
data.A1M0$M_diabetes <- 0

# data.A0M1$A0_PM2.5 <- 0
data.A0M1$M_diabetes <- 1

# data.A1M1$A0_PM2.5 <- 1
data.A1M1$M_diabetes <- 1

data.A0M0$Q.Y.death <- predict(Y.death.model, newdata = data.A0M0, type = "response")
data.A1M0$Q.Y.death <- predict(Y.death.model, newdata = data.A1M0, type = "response")
data.A0M1$Q.Y.death <- predict(Y.death.model, newdata = data.A0M1, type = "response")
data.A1M1$Q.Y.death <- predict(Y.death.model, newdata = data.A1M1, type = "response")

## 2a) Regress the predicted values conditional on the exposure A
##    and baseline confounders L(0)
L1.death.A0M0.model <- glm(Q.Y.death ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = data.A0M0)
L1.death.A1M0.model <- glm(Q.Y.death ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = data.A1M0)
L1.death.A0M1.model <- glm(Q.Y.death ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = data.A0M1)
L1.death.A1M1.model <- glm(Q.Y.death ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = data.A1M1)


## 2b) generate predicted values by evaluating the regression at exposure
##    of interest: {A=0,M=0}, {A=1,M=0}, {A=0,M=1}, {A=1,M=1}
data.A0M0$A0_PM2.5 <- 0
data.A1M0$A0_PM2.5 <- 1
data.A0M1$A0_PM2.5 <- 0
data.A1M1$A0_PM2.5 <- 1
data.A0M0$Yam.death.pred <- predict(L1.death.A0M0.model,
                                    newdata = data.A0M0, type = "response")
data.A1M0$Yam.death.pred <- predict(L1.death.A1M0.model,
                                    newdata = data.A1M0, type = "response")
data.A0M1$Yam.death.pred <- predict(L1.death.A0M1.model,
                                    newdata = data.A0M1, type = "response")
data.A1M1$Yam.death.pred <- predict(L1.death.A1M1.model,
                                    newdata = data.A1M1, type = "response")


## 3. Append the 4 counterfactual datasets in a single long dataset
# number of row is 4 times the initial value (we have 4 counterfactual scenarios)
data.4scenarios <- rbind(data.A0M0, data.A1M0,data.A0M1,data.A1M1)

## 4. fit the MSM: E(Y_am) = alpha_0 + alpha_A a + alpha_M m + alpha_AM a:m
MSM.CDE.gcomp <- glm(Yam.death.pred ~ A0_PM2.5 +  M_diabetes + A0_PM2.5:M_diabetes,
                     family = "gaussian", # gaussian family for risk differences
                     data = data.4scenarios)
coef(MSM.CDE.gcomp)
# (Intercept)           A0_PM2.5        M_diabetes A0_PM2.5:M_diabetes
#  0.17974978       0.06341297       0.07366464       0.02469211

## 5. Estimate the CDE(M=m)
# CDE(M=0) = E(Y_{A=1,M=0}) - E(Y_{A=0,M=0})
CDE_mis0_gcomp_ice <- coef(MSM.CDE.gcomp)["A0_PM2.5"]
# 0.06341297

# CDE(M=1) = E(Y_{A=1,M=1}) - E(Y_{A=0,M=1})
CDE_mis1_gcomp_ice <- (coef(MSM.CDE.gcomp)["A0_PM2.5"] +
                         coef(MSM.CDE.gcomp)["A0_PM2.5:M_diabetes"])
# 0.08810508


##### ltmle  --------------------------------------------------------------
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
library(ltmle)
Qform <- c(L1="Q.kplus1 ~ L0_male + L0_soc_env + A0_PM2.5",
           Y_death="Q.kplus1 ~ L0_male + L0_soc_env + L1 +
                    A0_PM2.5 * M_diabetes")
gform <- c("A0_PM2.5 ~ L0_male + L0_soc_env",
           "M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5 + L1")
data_binary <- subset(df2_int, select = c(L0_male, L0_soc_env,
                                          A0_PM2.5, L1,
                                          M_diabetes, Y_death))
CDE_ltmle_M0_death <- ltmle(data = data_binary,
                            Anodes = c("A0_PM2.5", "M_diabetes"),
                            Lnodes = c("L1"), # intermediate confounders +/- baseline
                            Ynodes = c("Y_death"),
                            survivalOutcome = FALSE, # TRUE for time-to-event outcomes Y
                            Qform = Qform,
                            gform = gform,
                            abar = list(c(1,0), # counterfactual intervention do(A=1,M=0)
                                        c(0,0)), # counterfactual intervention do(A=0,M=0)
                            SL.library = NULL,
                            estimate.time = FALSE, # estimate computation time?
                            gcomp = TRUE,
                            variance.method = "ic")
# CDE with M=0
summary(CDE_ltmle_M0_death)$effect.measures$ATE$estimate
#   Parameter Estimate:  0.06342833

CDE_ltmle_M1_death <- ltmle(data = data_binary,
                            Anodes = c("A0_PM2.5", "M_diabetes"),
                            Lnodes = c("L1"), # intermediate confounders +/- baseline
                            Ynodes = c("Y_death"),
                            survivalOutcome = FALSE, # TRUE for time-to-event outcomes Y
                            Qform = Qform,
                            gform = gform,
                            abar = list(c(1,1), # counterfactual intervention do(A=1,M=0)
                                        c(0,1)), # counterfactual intervention do(A=0,M=0)
                            SL.library = NULL,
                            estimate.time = FALSE, # estimate computation time?
                            gcomp = TRUE,
                            variance.method = "ic")
# CDE with M=1
summary(CDE_ltmle_M1_death)$effect.measures$ATE$estimate
#   Parameter Estimate:  0.08812318


rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")


# ---------------------------------------------------------------------------- #
# 3) MSM of NDE & NIE, estimated by IPTW ----
# ---------------------------------------------------------------------------- #
rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")

## 1. Stabilized weight for the MSM1
# 1a. sw_Ai = g(A=a_i | L(0)) / g(A=a_i | L(0)) = 1
sw_Ai <- rep(1, nrow(df1_int))

# 1b. sw_Mi = g(M=m_i | A,L(0)) / g(M=m_i | A,L(0),L(1))
g.M.AL0 <- glm(M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env,
               family = "binomial", data = df1_int)
g.Mis1.AL0 <- predict(g.M.AL0, type = "response")
sw_M.num <- rep(NA, nrow(df1_int))
sw_M.num[df1_int$M_diabetes==1] <- g.Mis1.AL0[df1_int$M_diabetes==1]
sw_M.num[df1_int$M_diabetes==0] <- (1 - g.Mis1.AL0[df1_int$M_diabetes==0])

g.M.AL0L1 <- glm(M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env + L1,
                 family = "binomial", data = df1_int)
g.Mis1.AL0L1 <- predict(g.M.AL0L1, type = "response")
sw_M.denom <- rep(NA, nrow(df1_int))
sw_M.denom[df1_int$M_diabetes==1] <- g.Mis1.AL0L1[df1_int$M_diabetes==1]
sw_M.denom[df1_int$M_diabetes==0] <- (1 - g.Mis1.AL0L1[df1_int$M_diabetes==0])

sw_msm1 <- sw_Ai * sw_M.num / sw_M.denom

## 2. Estimate coefficients of the MSM1
MSM1 <- glm(Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5:M_diabetes +
              L0_male + L0_soc_env,
            weights = sw_msm1,
            family = "gaussian",
            data = df1_int)
coef(MSM1)
# (Intercept)                A0_PM2.5             M_diabetes
#  0.12033221            0.06381257            0.06691712
#     L0_male L0_soc_env      A0_PM2.5:M_diabetes
#  0.04671886            0.05521263            0.01652446

## 3. Stabilized weight for the MSM2
# 3a. sw_A = g(A=a_i) / g(A=a_i | L(0))
# numerator
g.A <- glm(A0_PM2.5 ~ 1, family = "binomial", data = df1_int)
g.Ais1 <- predict(g.A, type = "response")
sw_msm2.num <- rep(NA, nrow(df1_int))
sw_msm2.num[df1_int$A0_PM2.5==1] <- g.Ais1[df1_int$A0_PM2.5==1]
sw_msm2.num[df1_int$A0_PM2.5==0] <- (1 - g.Ais1[df1_int$A0_PM2.5==0])

# denominator
g.A.L0 <- glm(A0_PM2.5 ~ L0_male + L0_soc_env,
              family = "binomial", data = df1_int)
g.Ais1.L0 <- predict(g.A.L0, type = "response")
sw_msm2.denom <- rep(NA, nrow(df1_int))
sw_msm2.denom[df1_int$A0_PM2.5==1] <- g.Ais1.L0[df1_int$A0_PM2.5==1]
sw_msm2.denom[df1_int$A0_PM2.5==0] <- (1 - g.Ais1.L0[df1_int$A0_PM2.5==0])

# stabilized weight
sw_msm2 <- sw_msm2.num / sw_msm2.denom

## 3. Estimate coefficients of the MSM2
MSM2 <- glm(M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env,
            weights = sw_msm2,
            family = "binomial",
            data = df1_int)
coef(MSM2)
# (Intercept)                A0_PM2.5               L0_male L0_soc_env
#  -1.2723106             0.5883720             0.2566129             0.3270087

## 4. Estimate PNDE conditional on L(0), and the marginal value of PNDE
# a = 1 and a* = 0
# PNDE|L(0) = (a - a*)[alpha_A + alpha_AM.g^-1(a^*,l(0))]
g.minus1.A0 <- plogis(coef(MSM2)["(Intercept)"] + coef(MSM2)["A0_PM2.5"] * 0 +
                      coef(MSM2)["L0_male"] * df1_int$L0_male +
                      coef(MSM2)["L0_soc_env"] * df1_int$L0_soc_env)

# PNDE conditional on L(0)
PNDE_L0 <- (1 - 0) * (coef(MSM1)["A0_PM2.5"] +
                        coef(MSM1)["A0_PM2.5:M_diabetes"] * g.minus1.A0)
# marginal PNDE
PNDE <- mean(PNDE_L0)
# [1] 0.06850657

## 4. Estimate TNIE conditional on L(0), and the marginal value of TNIE
# TNIE|L(0) = [g^-1(a,l(0)) - g^-1(a^*,l(0))] * (alpha_M + alpha_AM * a)
g.minus1.A1 <- plogis(coef(MSM2)["(Intercept)"] + coef(MSM2)["A0_PM2.5"] * 1 +
                        coef(MSM2)["L0_male"] * df1_int$L0_male +
                        coef(MSM2)["L0_soc_env"] * df1_int$L0_soc_env)

# TNIE conditional on L(0)
TNIE_L0 <- (g.minus1.A1 - g.minus1.A0) * (coef(MSM1)["M_diabetes"] +
                                            coef(MSM1)["A0_PM2.5:M_diabetes"] * 1)
# marginal PNDE
TNIE <- mean(TNIE_L0)
# [1] 0.01096799


# ---------------------------------------------------------------------------- #
# IV) CMAverse -----------
# ---------------------------------------------------------------------------- #
## IV.1) without intermediate confounders affected by the exposure ----
# ---------------------------------------------------------------------------- #
rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")
library(CMAverse)
## The Marginal Structural Model
res_msm_df1 <- cmest(data = df1_int,
                     model = "msm",
                     outcome = "Y_death",
                     exposure = "A0_PM2.5",
                     mediator = "M_diabetes",
                     basec = c("L0_male", "L0_soc_env","L1"),
                     postc = NULL,
                     EMint = TRUE, # E*M interaction
                     ereg = "logistic", # exposure regression model g(A=1|L(0))
                     yreg = "linear", # to get risk difference
                     mreg = list("logistic"), # mediation model g(M=1|L1,A,L0)
                     wmnomreg = list("logistic"), #g(M=1|A) wgt nominator
                     wmdenomreg = list("logistic"), # g(M=1|L1,A,L(0)) wgt denom
                     astar = 0, #E(Y_{A=0,M=1})
                     a = 1,  #E(Y_{A=1,M=1})
                     mval = list(1), # mediator value at which the variable is controlled
                     estimation = "imputation",
                     inference = "bootstrap",
                     nboot = 2)
summary(res_msm_df1)
# Causal Mediation Analysis
#
# # Outcome regression:
# Call:
#   glm(formula = Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5 * M_diabetes,
#       family = gaussian(), data = getCall(x$reg.output$yreg)$data,
#       weights = getCall(x$reg.output$yreg)$weights)
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
#   (Intercept)          0.178967   0.005049  35.445  < 2e-16 ***
#   A0_PM2.5             0.067151   0.016679   4.026 5.72e-05 ***
#   M_diabetes           0.067321   0.009508   7.080 1.54e-12 ***
#   A0_PM2.5:M_diabetes -0.004491   0.026027  -0.173    0.863
#
# # Mediator regressions:
# Call:
#   glm(formula = M_diabetes ~ A0_PM2.5, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -0.92405    0.02353 -39.264   <2e-16 ***
#   A0_PM2.5     0.58413    0.06501   8.986   <2e-16 ***
#
# # Mediator regressions for weighting (denominator):
# Call:
#   glm(formula = M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env +
#         L1, family = binomial(), data = getCall(x$reg.output$wmdenomreg[[1L]])$data,
#       weights = getCall(x$reg.output$wmdenomreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -1.37880    0.04872 -28.303  < 2e-16 ***
#   A0_PM2.5     0.56260    0.06555   8.583  < 2e-16 ***
#   L0_male      0.25861    0.04437   5.829 5.57e-09 ***
#   L0_soc_env   0.33050    0.04744   6.967 3.23e-12 ***
#   L1           0.33462    0.04744   7.054 1.74e-12 ***
#
# # Mediator regressions for weighting (nominator):
# Call:
#   glm(formula = M_diabetes ~ A0_PM2.5, family = binomial(), data = getCall(x$reg.output$wmnomreg[[1L]])$data,
#       weights = getCall(x$reg.output$wmnomreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -0.93236    0.02358 -39.544   <2e-16 ***
#   A0_PM2.5     0.62388    0.06481   9.627   <2e-16 ***
#
# # Exposure regression for weighting:
# Call:
#   glm(formula = A0_PM2.5 ~ L0_male + L0_soc_env + L1, family = binomial(),
#       data = getCall(x$reg.output$ereg)$data, weights = getCall(x$reg.output$ereg)$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.74236    0.07729 -35.484  < 2e-16 ***
#   L0_male      0.40610    0.06448   6.298 3.01e-10 ***
#   L0_soc_env   0.64079    0.07350   8.718  < 2e-16 ***
#   L1           0.03257    0.06968   0.467     0.64
#
# # Effect decomposition on the mean difference scale via the marginal structural model
#
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
# Estimate  Std.error    95% CIL 95% CIU  P.val
#   cde           0.0626596  0.0241140  0.0283293   0.061 <2e-16 ***
#   pnde          0.0658736  0.0172598  0.0327941   0.056 <2e-16 ***
#   tnde          0.0652848  0.0183434  0.0319353   0.057 <2e-16 ***
#   pnie          0.0088254  0.0033214  0.0057592   0.010 <2e-16 ***
#   tnie          0.0082366  0.0022378  0.0063562   0.009 <2e-16 ***
#   te            0.0741102  0.0150221  0.0421568   0.062 <2e-16 ***
#   intref        0.0032140  0.0068542 -0.0047438   0.004      1
#   intmed       -0.0005888  0.0010836 -0.0008588   0.001      1
#   cde(prop)     0.8454925  0.2258321  0.6680758   0.971 <2e-16 ***
#   intref(prop)  0.0433679  0.1360406 -0.0744993   0.108      1
#   intmed(prop) -0.0079448  0.0223848 -0.0207595   0.009      1
#   pnie(prop)    0.1190844  0.1121762  0.0937028   0.244 <2e-16 ***
#   pm            0.1111396  0.0897914  0.1030174   0.224 <2e-16 ***
#   int           0.0354231  0.1136558 -0.0651847   0.088      1
#   pe            0.1545075  0.2258321  0.0285181   0.332 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (cde: controlled direct effect;
#  pnde: pure natural direct effect;
#  tnde: total natural direct effect;
#  pnie: pure natural indirect effect;
#  tnie: total natural indirect effect;
#  te: total effect;
#  intref: reference interaction;
#  intmed: mediated interaction;
#  cde(prop): proportion cde;
#  intref(prop): proportion intref;
#  intmed(prop): proportion intmed;
#  pnie(prop): proportion pnie;
#  pm: overall proportion mediated;
#  int: overall proportion attributable to interaction;
#  pe: overall proportion eliminated)










