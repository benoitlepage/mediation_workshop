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
# 3) MSM (conditional on L0,L1) of NDE & NIE, estimated by IPTW ----
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
# 4) CMAverse -----------
# ---------------------------------------------------------------------------- #
## 4.0) Using the mediation formula and MSM of the outcome and the mediator ----

# E(Y_{a,M_a*}) = sum_m E(Y_{am}) * P(M_a* = m)

# We can estimate 2 MSMs, using weighted regressions:
#  - E(Y_am) = alpha_0 + alpha_1 a + alpha_2 m + alpha_3 (a*m)
#  - P(M_a*) = expit(beta_0 + beta_1 a*)

rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")

## 1. Stabilized weight for the exposure sw_{A,i}
# 1a. Estimate g(A=a_i|L(0)) (denominator of the weight)
g.A.L <- glm(A0_PM2.5 ~ L0_male + L0_soc_env + L1,
             family = "binomial", data = df1_int)
summary(g.A.L)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -2.74236    0.07729 -35.484  < 2e-16 ***
#   L0_male      0.40610    0.06448   6.298 3.01e-10 ***
#   L0_soc_env   0.64079    0.07350   8.718  < 2e-16 ***
#   L1           0.03257    0.06968   0.467     0.64

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
summary(g.M.L)
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -1.37880    0.04872 -28.303  < 2e-16 ***
#   L0_male      0.25861    0.04437   5.829 5.57e-09 ***
#   L0_soc_env   0.33050    0.04744   6.967 3.23e-12 ***
#   A0_PM2.5     0.56260    0.06555   8.583  < 2e-16 ***
#   L1           0.33462    0.04744   7.054 1.74e-12 ***

# 2b. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i | L(0)) is :
gMi.L <- rep(NA, nrow(df1_int))
gMi.L[df1_int$M_diabetes==1] <- predict(g.M.L, type="response")[df1_int$M_diabetes==1]
gMi.L[df1_int$M_diabetes==0] <- (1 - predict(g.M.L, type="response"))[df1_int$M_diabetes==0]

# 2c. Estimate g(M=m_i|A) (numerator of the weight)
g.M.A <- glm(M_diabetes ~ A0_PM2.5, family = "binomial", data = df1_int)
summary(g.M.A)
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -0.93236    0.02358 -39.544   <2e-16 ***
#   A0_PM2.5     0.62388    0.06481   9.627   <2e-16 ***

# 2d. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(M = m_i|A) is :
gMi.A <- rep(NA, nrow(df1_int))
gMi.A[df1_int$M_diabetes==1] <- predict(g.M.A, type="response")[df1_int$M_diabetes==1]
gMi.A[df1_int$M_diabetes==0] <- (1 - predict(g.M.A, type="response"))[df1_int$M_diabetes==0]
# 2e. Calculate sw_{M,i}
sw_Mi <- gMi.A / gMi.L

## 3. Estimate marginal model for the counterfactual outcome Y_{am}
model.Yam.death <- glm(Y_death ~ A0_PM2.5 * M_diabetes,
                       weights = sw_Ai * sw_Mi,
                       family = "gaussian",
                       data = df1_int)
summary(model.Yam.death)
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
#   (Intercept)          0.178967   0.005049  35.445  < 2e-16 ***
#   A0_PM2.5             0.067151   0.016679   4.026 5.72e-05 ***
#   M_diabetes           0.067321   0.009508   7.080 1.54e-12 ***
#   A0_PM2.5:M_diabetes -0.004491   0.026027  -0.173    0.863

## 4) Estimate marginal model for the counterfactual mediator P(M_a* = 1)
model.Ma <- glm(M_diabetes ~ A0_PM2.5,
                weights = sw_Ai,
                family = "binomial",
                data = df1_int)
summary(model.Ma)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -0.92405    0.02353 -39.264   <2e-16 ***
#   A0_PM2.5     0.58413    0.06501   8.986   <2e-16 ***

## 5) Estimate population average counterfactuals using the "mediation formula"
##    E(Y_{a,M_a*}) = sum_m E(Y_{am}) * P(M_a^* = m)
E.Y0M0 <- ((coef(model.Yam.death)["(Intercept)"] + # sum_m E(Y_{a0}) * P(M_a^* = 0)
             coef(model.Yam.death)["A0_PM2.5"] * 0 +
             coef(model.Yam.death)["M_diabetes"] * 0 +
             coef(model.Yam.death)["A0_PM2.5:M_diabetes"] * 0 * 0) *
             (1 - plogis(coef(model.Ma)["(Intercept)"] +
                           coef(model.Ma)["A0_PM2.5"] * 0))) +
  ((coef(model.Yam.death)["(Intercept)"] +  # sum_m E(Y_{a1}) * P(M_a^* = 1)
      coef(model.Yam.death)["A0_PM2.5"] * 0 +
      coef(model.Yam.death)["M_diabetes"] * 1 +
      coef(model.Yam.death)["A0_PM2.5:M_diabetes"] * 0 * 1) *
     (plogis(coef(model.Ma)["(Intercept)"] +
               coef(model.Ma)["A0_PM2.5"] * 0)))

E.Y1M0 <- ((coef(model.Yam.death)["(Intercept)"] + # sum_m E(Y_{a,M=0}) * P(M_0 = 0)
              coef(model.Yam.death)["A0_PM2.5"] * 1 +
              coef(model.Yam.death)["M_diabetes"] * 0 +
              coef(model.Yam.death)["A0_PM2.5:M_diabetes"] * 1 * 0) *
             (1 - plogis(coef(model.Ma)["(Intercept)"] +
                           coef(model.Ma)["A0_PM2.5"] * 0))) +
  ((coef(model.Yam.death)["(Intercept)"] +  # sum_m E(Y_{a,M=1}) * P(M_0 = 1)
      coef(model.Yam.death)["A0_PM2.5"] * 1 +
      coef(model.Yam.death)["M_diabetes"] * 1 +
      coef(model.Yam.death)["A0_PM2.5:M_diabetes"] * 1 * 1) *
     (plogis(coef(model.Ma)["(Intercept)"] +
               coef(model.Ma)["A0_PM2.5"] * 0)))

E.Y1M1 <- ((coef(model.Yam.death)["(Intercept)"] + #  E(Y_{a,M=0}) * P(M_1 = 0)
              coef(model.Yam.death)["A0_PM2.5"] * 1 +
              coef(model.Yam.death)["M_diabetes"] * 0 +
              coef(model.Yam.death)["A0_PM2.5:M_diabetes"] * 1 * 0) *
             (1 - plogis(coef(model.Ma)["(Intercept)"] +
                           coef(model.Ma)["A0_PM2.5"] * 1))) +
  ((coef(model.Yam.death)["(Intercept)"] +  #  E(Y_{a,M=1}) * P(M_1 = 1)
      coef(model.Yam.death)["A0_PM2.5"] * 1 +
      coef(model.Yam.death)["M_diabetes"] * 1 +
      coef(model.Yam.death)["A0_PM2.5:M_diabetes"] * 1 * 1) *
     (plogis(coef(model.Ma)["(Intercept)"] +
               coef(model.Ma)["A0_PM2.5"] * 1)))

PNDE.death <- E.Y1M0 - E.Y0M0
# 0.06587481
TNIE.death <- E.Y1M1 - E.Y1M0
# 0.008274351


## 6) Estimate population average counterfactuals using simulated
##    values of the mediator
PNDE.death.sim <- rep(NA,5)
TNIE.death.sim <- rep(NA,5)

set.seed(1234)
for(k in 1:15){
  M0 <- rbinom(n = nrow(df1_int), size = 1,
               prob = plogis(coef(model.Ma)["(Intercept)"] +
                               coef(model.Ma)["A0_PM2.5"] * 0))
  M1 <- rbinom(n = nrow(df1_int), size = 1,
               prob = plogis(coef(model.Ma)["(Intercept)"] +
                               coef(model.Ma)["A0_PM2.5"] * 1))
  E.Y0M0 <- mean((coef(model.Yam.death)["(Intercept)"] +
                    coef(model.Yam.death)["A0_PM2.5"] * 0 +
                    coef(model.Yam.death)["M_diabetes"] * M0 +
                    coef(model.Yam.death)["A0_PM2.5:M_diabetes"] * 0 * M0))
  E.Y1M0 <- mean((coef(model.Yam.death)["(Intercept)"] +
                    coef(model.Yam.death)["A0_PM2.5"] * 1 +
                    coef(model.Yam.death)["M_diabetes"] * M0 +
                    coef(model.Yam.death)["A0_PM2.5:M_diabetes"] * 1 * M0))
  E.Y1M1 <- mean((coef(model.Yam.death)["(Intercept)"] +
                coef(model.Yam.death)["A0_PM2.5"] * 1 +
                coef(model.Yam.death)["M_diabetes"] * M1 +
                coef(model.Yam.death)["A0_PM2.5:M_diabetes"] * 1 * M1))
  PNDE.death.sim[k] <- E.Y1M0  - E.Y0M0
  TNIE.death.sim[k] <- E.Y1M1 - E.Y1M0
}
mean(PNDE.death.sim)
# [1] 0.06587304
mean(TNIE.death.sim)
# [1] 0.00824998

## 4.1) Estimation of NDE and NIE using CMAverse (by IPTW or gcomputation) ----
##      without intermediate confounders affected by the exposure
# ---------------------------------------------------------------------------- #
rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")
library(CMAverse)

### CMAverse MSM, by IPTW ----
## Using the CMAverse to estimate MSM estimated by IPTW
set.seed(1234) # note: there is randomness, even with estimation by IPTW
               #       => counterfactuals are imputed (estimation = "imputation")
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
                     mval = list(0), # mediator value at which the variable is controlled
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
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -2.74236    0.07729 -35.484  < 2e-16 ***
#   L0_male      0.40610    0.06448   6.298 3.01e-10 ***
#   L0_soc_env   0.64079    0.07350   8.718  < 2e-16 ***
#   L1           0.03257    0.06968   0.467     0.64
#
# # Effect decomposition on the mean difference scale via the marginal structural model
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
#                  Estimate  Std.error    95% CIL 95% CIU  P.val
#   cde           0.0671509  0.0157808  0.0761499   0.097 <2e-16 ***
#   pnde          0.0658727  0.0169713  0.0643983   0.087 <2e-16 ***
#   tnde          0.0652889  0.0172698  0.0594343   0.083 <2e-16 ***
#   pnie          0.0087514  0.0013389  0.0070189   0.009 <2e-16 ***
#   tnie          0.0081675  0.0010405  0.0024559   0.004 <2e-16 ***
#   te            0.0740402  0.0159309  0.0682520   0.090 <2e-16 ***
#   intref       -0.0012782  0.0011906 -0.0117517  -0.010 <2e-16 ***
#   intmed       -0.0005839  0.0002984 -0.0049640  -0.005 <2e-16 ***
#   cde(prop)     0.9069522  0.0222813  1.0860362   1.116 <2e-16 ***
#   intref(prop) -0.0172642  0.0439631 -0.1726807  -0.114 <2e-16 ***
#   intmed(prop) -0.0078856  0.0162851 -0.0729153  -0.051 <2e-16 ***
#   pnie(prop)    0.1181976  0.0379668  0.0786163   0.130 <2e-16 ***
#   pm            0.1103120  0.0216818  0.0275800   0.057 <2e-16 ***
#   int          -0.0251498  0.0602482 -0.2455960  -0.165 <2e-16 ***
#   pe            0.0930478  0.0222813 -0.1159712  -0.086 <2e-16 ***
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

plot(res_msm_df1$reg.output$yreg$weights, sw_Ai * sw_Mi)
abline(0,1)

plot(res_msm_df1$reg.output$mreg[[1]]$weights, sw_Ai) # not the same weights for the mediator MSM ?
abline(0,1)
# what are those weights ??
head(data.frame(df1_int, wA = res_msm_df1$reg.output$mreg[[1]]$weights), 20)
head(data.frame(df1_int[,1:4],
                wA = res_msm_df1$reg.output$mreg[[1]]$weights,
                test = predict(g.A.L, type="response")), 20) #predict(g.A, type="response") / predict(g.A.L, type="response")
# Should check in the code what are the weights for the mreg model

# We could have used direct the "mediation formula" (without imputing counterfactuals)
E.Y0M0 <- ((coef(res_msm_df1$reg.output$yreg)["(Intercept)"] + # sum_m E(Y_{a0}) * P(M_a^* = 0)
              coef(res_msm_df1$reg.output$yreg)["A0_PM2.5"] * 0 +
              coef(res_msm_df1$reg.output$yreg)["M_diabetes"] * 0 +
              coef(res_msm_df1$reg.output$yreg)["A0_PM2.5:M_diabetes"] * 0 * 0) *
             (1 - plogis(coef(res_msm_df1$reg.output$mreg[[1]])["(Intercept)"] +
                           coef(res_msm_df1$reg.output$mreg[[1]])["A0_PM2.5"] * 0))) +
  ((coef(res_msm_df1$reg.output$yreg)["(Intercept)"] +  # sum_m E(Y_{a1}) * P(M_a^* = 1)
      coef(res_msm_df1$reg.output$yreg)["A0_PM2.5"] * 0 +
      coef(res_msm_df1$reg.output$yreg)["M_diabetes"] * 1 +
      coef(res_msm_df1$reg.output$yreg)["A0_PM2.5:M_diabetes"] * 0 * 1) *
     (plogis(coef(res_msm_df1$reg.output$mreg[[1]])["(Intercept)"] +
               coef(res_msm_df1$reg.output$mreg[[1]])["A0_PM2.5"] * 0)))

E.Y1M0 <- ((coef(res_msm_df1$reg.output$yreg)["(Intercept)"] + # sum_m E(Y_{a,M=0}) * P(M_0 = 0)
              coef(res_msm_df1$reg.output$yreg)["A0_PM2.5"] * 1 +
              coef(res_msm_df1$reg.output$yreg)["M_diabetes"] * 0 +
              coef(res_msm_df1$reg.output$yreg)["A0_PM2.5:M_diabetes"] * 1 * 0) *
             (1 - plogis(coef(res_msm_df1$reg.output$mreg[[1]])["(Intercept)"] +
                           coef(res_msm_df1$reg.output$mreg[[1]])["A0_PM2.5"] * 0))) +
  ((coef(res_msm_df1$reg.output$yreg)["(Intercept)"] +  # sum_m E(Y_{a,M=1}) * P(M_0 = 1)
      coef(res_msm_df1$reg.output$yreg)["A0_PM2.5"] * 1 +
      coef(res_msm_df1$reg.output$yreg)["M_diabetes"] * 1 +
      coef(res_msm_df1$reg.output$yreg)["A0_PM2.5:M_diabetes"] * 1 * 1) *
     (plogis(coef(res_msm_df1$reg.output$mreg[[1]])["(Intercept)"] +
               coef(res_msm_df1$reg.output$mreg[[1]])["A0_PM2.5"] * 0)))

E.Y1M1 <- ((coef(res_msm_df1$reg.output$yreg)["(Intercept)"] + #  E(Y_{a,M=0}) * P(M_1 = 0)
              coef(res_msm_df1$reg.output$yreg)["A0_PM2.5"] * 1 +
              coef(res_msm_df1$reg.output$yreg)["M_diabetes"] * 0 +
              coef(res_msm_df1$reg.output$yreg)["A0_PM2.5:M_diabetes"] * 1 * 0) *
             (1 - plogis(coef(res_msm_df1$reg.output$mreg[[1]])["(Intercept)"] +
                           coef(res_msm_df1$reg.output$mreg[[1]])["A0_PM2.5"] * 1))) +
  ((coef(res_msm_df1$reg.output$yreg)["(Intercept)"] +  #  E(Y_{a,M=1}) * P(M_1 = 1)
      coef(res_msm_df1$reg.output$yreg)["A0_PM2.5"] * 1 +
      coef(res_msm_df1$reg.output$yreg)["M_diabetes"] * 1 +
      coef(res_msm_df1$reg.output$yreg)["A0_PM2.5:M_diabetes"] * 1 * 1) *
     (plogis(coef(res_msm_df1$reg.output$mreg[[1]])["(Intercept)"] +
               coef(res_msm_df1$reg.output$mreg[[1]])["A0_PM2.5"] * 1)))

PNDE.death <- E.Y1M0 - E.Y0M0
# 0.06587481
TNIE.death <- E.Y1M1 - E.Y1M0
# 0.008274351

### CMAverse MSM, by g-comp ----
## Using the CMAverse to estimate MSM estimated by g-comp
set.seed(1234) # note: there is randomness, even with estimation by IPTW
#       => counterfactuals are imputed (estimation = "imputation")
res_msm_df1.qol <- cmest(data = df1_int,
                     model = "gformula",
                     outcome = "Y_qol",
                     exposure = "A0_PM2.5",
                     mediator = "M_diabetes",
                     basec = c("L0_male", "L0_soc_env","L1"),
                     postc = NULL,
                     EMint = TRUE, # E*M interaction
                     # ereg = "logistic", # # not needed with gcomp
                     yreg = "linear", # to get risk difference
                     mreg = list("logistic"), # mediation model g(M=1|L1,A,L0)
                     # wmnomreg = list("logistic"), # not needed with gcomp
                     # wmdenomreg = list("logistic"), ## not needed with gcomp
                     astar = 0, #E(Y_{A=0,M=1})
                     a = 1,  #E(Y_{A=1,M=1})
                     mval = list(0), # mediator value at which the variable is controlled
                     estimation = "imputation",
                     inference = "bootstrap",
                     nboot = 2)
summary(res_msm_df1.qol)
# Causal Mediation Analysis
#
# # Outcome regression:
# Call:
#   glm(formula = Y_qol ~ A0_PM2.5 + M_diabetes + A0_PM2.5 * M_diabetes +
#         L0_male + L0_soc_env + L1, family = gaussian(), data = getCall(x$reg.output$yreg)$data,
#       weights = getCall(x$reg.output$yreg)$weights)
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)
#   (Intercept)          74.7669     0.2139 349.557  < 2e-16 ***
#   A0_PM2.5             -3.7153     0.4160  -8.931  < 2e-16 ***
#   M_diabetes           -8.6317     0.2385 -36.197  < 2e-16 ***
#   L0_male              -0.7235     0.2018  -3.586 0.000337 ***
#   L0_soc_env           -2.8899     0.2112 -13.684  < 2e-16 ***
#   L1                   -3.4280     0.2212 -15.494  < 2e-16 ***
#   A0_PM2.5:M_diabetes  -5.6154     0.6514  -8.621  < 2e-16 ***
#
# # Mediator regressions:
# Call:
#   glm(formula = M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env +
#         L1, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -1.37880    0.04872 -28.303  < 2e-16 ***
#   A0_PM2.5     0.56260    0.06555   8.583  < 2e-16 ***
#   L0_male      0.25861    0.04437   5.829 5.57e-09 ***
#   L0_soc_env   0.33050    0.04744   6.967 3.23e-12 ***
#   L1           0.33462    0.04744   7.054 1.74e-12 ***
#
# # Effect decomposition on the mean difference scale via the g-formula approach
#
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
#                 Estimate Std.error   95% CIL 95% CIU  P.val
#   cde          -3.715265  0.114707 -3.084527  -2.930 <2e-16 ***
#   pnde         -5.103390  0.008919 -4.820137  -4.808 <2e-16 ***
#   tnde         -5.659875  0.006095 -5.579563  -5.571 <2e-16 ***
#   pnie         -0.855400  0.040149 -1.006846  -0.953 <2e-16 ***
#   tnie         -1.411886  0.025135 -1.758083  -1.724 <2e-16 ***
#   te           -6.515276  0.034054 -6.578220  -6.532 <2e-16 ***
#   intref       -1.388125  0.105788 -1.877736  -1.736 <2e-16 ***
#   intmed       -0.556485  0.015014 -0.771408  -0.751 <2e-16 ***
#   cde(prop)     0.570239  0.015115  0.448589   0.469 <2e-16 ***
#   intref(prop)  0.213057  0.017570  0.263846   0.287 <2e-16 ***
#   intmed(prop)  0.085412  0.002894  0.114201   0.118 <2e-16 ***
#   pnie(prop)    0.131292  0.005348  0.145871   0.153 <2e-16 ***
#   pm            0.216704  0.002454  0.263960   0.267 <2e-16 ***
#   int           0.298469  0.020463  0.378048   0.406 <2e-16 ***
#   pe            0.429761  0.015115  0.531104   0.551 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (cde: controlled direct effect;
#  pnde: pure natural direct effect;
#  tnde: total natural direct effect;
#  pnie: pure natural indirect effect;
#  tnie: total natural indirect effect;
#  te: total effect; intref:
#  reference interaction;
#  intmed: mediated interaction;
#  cde(prop): proportion cde;
#  intref(prop): proportion intref;
#  intmed(prop): proportion intmed;
#  pnie(prop): proportion pnie;
#  pm: overall proportion mediated;
#  int: overall proportion attributable to interaction;
#  pe: overall proportion eliminated)


# ---------------------------------------------------------------------------- #
# 5) Unified approach based on MSM, estimated by IPTW ----
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
## 5.1) manual calculation ----
# ---------------------------------------------------------------------------- #
# from Theis Lange, Stijn Vansteelandt, Maarten Bekaert, A Simple Unified Approach
# for Estimating Natural Direct and Indirect Effects. Am J Epidemiol 2012;176,(3)190–195.
# https://doi-org.proxy.insermbiblio.inist.fr/10.1093/aje/kwr525
rm(list=ls())
# this method is used to estimate Natural direct and indirect effects when they are
# identifiable
df1_int <- read.csv(file = "data/df1_int.csv")

## Step 1. Estimate a suitable model for the exposure conditional on confounders,
##         using the original data set
# for A weight numerator
g.A <- glm(A0_PM2.5 ~ 1,
           family = "binomial", data = df1_int)

# for A weight denominator
g.A.L0 <- glm(A0_PM2.5 ~ L0_male + L0_soc_env + L1, # L1 not necessary in our example
              family = "binomial", data = df1_int)

## Step 2. Estimate a suitable model for the mediator conditional on confounders,
##         using the original data set
df1_int$A0_PM2.5_temp <- df1_int$A0_PM2.5
g.M.AL0 <- glm(M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5_temp + L1,
               family = "binomial", data = df1_int)

## Step 3. Construct a new data set by repeating the original data set twice,
##         and including an additional variable A*, equal to the original exposure
##         in the 1st data set, and equal to 1-A in the 2nd data set
##         (note for multicategorical exposure, see supplementary material of Lange et al)

# create identifier
df1_int$id <- 1:nrow(df1_int)

# dupplicate the original data set
df1_int.1 <- df1_int.2 <- df1_int

# Add an A* variable
df1_int.1$A0_PM2.5_star <- df1_int$A0_PM2.5
df1_int.2$A0_PM2.5_star <- 1 - df1_int$A0_PM2.5

# stack the 2 new datasets
df1_int.new <- rbind(df1_int.1,df1_int.2)
# sort the data set by the id-variable to use geeglm
# (necessary to use geepack and the geeglm() later)
df1_int.new <- df1_int.new[order(df1_int.new$id), ]

head(df1_int.new, 12)
#       L0_male L0_soc_env A0_PM2.5 L1 M_diabetes Y_death Y_qol A0_PM2.5_temp id A0_PM2.5_star
# 1           0          1        0  1          0       0  93.4             0  1             0
# 10001       0          1        0  1          0       0  93.4             0  1             1
# 2           1          1        1  1          0       1  64.0             1  2             1
# 10002       1          1        1  1          0       1  64.0             1  2             0
# 3           1          1        0  0          0       0  75.6             0  3             0
# 10003       1          1        0  0          0       0  75.6             0  3             1
# 4           1          0        0  0          0       0  89.8             0  4             0
# 10004       1          0        0  0          0       0  89.8             0  4             1
# 5           1          1        0  0          0       0  77.2             0  5             0
# 10005       1          1        0  0          0       0  77.2             0  5             1
# 6           1          1        0  0          1       0  73.9             0  6             0
# 10006       1          1        0  0          1       0  73.9             0  6             1


## Step 4. Compute weights by applying the fitted models from steps 1 and 2
# weights for A
g.Ais1 <- predict(g.A, newdata = df1_int.new, type = "response")
sw_A.num <- rep(NA, nrow(df1_int.new))
sw_A.num[df1_int.new$A0_PM2.5==1] <- g.Ais1[df1_int.new$A0_PM2.5==1]
sw_A.num[df1_int.new$A0_PM2.5==0] <- (1 - g.Ais1[df1_int.new$A0_PM2.5==0])

g.Ais1.L0 <- predict(g.A.L0, newdata = df1_int.new, type = "response")
sw_A.denom <- rep(NA, nrow(df1_int.new))
sw_A.denom[df1_int.new$A0_PM2.5==1] <- g.Ais1.L0[df1_int.new$A0_PM2.5==1]
sw_A.denom[df1_int.new$A0_PM2.5==0] <- (1 - g.Ais1.L0[df1_int.new$A0_PM2.5==0])

w_A <- 1 / sw_A.denom
sw_A <- sw_A.num / sw_A.denom

# weights for M|A*
df1_int.new$A0_PM2.5_temp <- df1_int.new$A0_PM2.5_star
g.Mis1.Astar <- predict(g.M.AL0, newdata = df1_int.new, type = "response")
sw_M.Astar <- rep(NA, nrow(df1_int.new))
sw_M.Astar[df1_int.new$M_diabetes==1] <- g.Mis1.Astar[df1_int.new$M_diabetes==1]
sw_M.Astar[df1_int.new$M_diabetes==0] <- 1 - g.Mis1.Astar[df1_int.new$M_diabetes==0]

# weight for M|A
df1_int.new$A0_PM2.5_temp <- df1_int.new$A0_PM2.5
g.Mis1.A <- predict(g.M.AL0, newdata = df1_int.new, type = "response")
sw_M.A <- rep(NA, nrow(df1_int.new))
sw_M.A[df1_int.new$M_diabetes==1] <- g.Mis1.A[df1_int.new$M_diabetes==1]
sw_M.A[df1_int.new$M_diabetes==0] <- 1 - g.Mis1.A[df1_int.new$M_diabetes==0]

# non-stabilized weight
w <- w_A * (sw_M.Astar / sw_M.A)
# stabilized weight
sw <- sw_A * (sw_M.Astar / sw_M.A)

boxplot(data.frame(w, sw))

## Step 5. Fit a suitable model to the outcome including only A and A* (and their
##         interaction)
# we will use a GEE model in order to get "robust" standard errors
library(geepack)
# note: be careful that individuals in the "df1_int.new" data set
#       and in the "w" vector have the same order

## for a risk difference in death
MSM.model.death <- geeglm(Y_death ~ A0_PM2.5 * A0_PM2.5_star,
                          weights = sw,
                          family = "gaussian",
                          id = df1_int.new$id,
                          data = df1_int.new,
                          scale.fix = TRUE)
summary(MSM.model.death)
# Coefficients:
#                           Estimate   Std.err    Wald Pr(>|W|)
#   (Intercept)             0.198840  0.004251 2187.61  < 2e-16 ***
#   A0_PM2.5                0.065950  0.014539   20.58  5.7e-06 ***
#   A0_PM2.5_star           0.008534  0.001242   47.25  6.2e-12 ***
#   A0_PM2.5:A0_PM2.5_star -0.000406  0.003762    0.01     0.91
unif.PNDE.death <- coef(MSM.model.death)["A0_PM2.5"]
# 0.066
unif.TNIE.death <- (coef(MSM.model.death)["A0_PM2.5_star"] +
                      coef(MSM.model.death)["A0_PM2.5:A0_PM2.5_star"])
# 0.00813

## For Quality of life
MSM.model.qol <- geeglm(Y_qol ~ A0_PM2.5 * A0_PM2.5_star,
                        weights = sw,
                        family = "gaussian",
                        id = df1_int.new$id,
                        data = df1_int.new,
                        scale.fix = TRUE)
summary(MSM.model.qol)
# Coefficients:
#                          Estimate Std.err     Wald Pr(>|W|)
#   (Intercept)             69.0853  0.1174 346240.4  < 2e-16 ***
#   A0_PM2.5                -5.3771  0.4010    179.8  < 2e-16 ***
#   A0_PM2.5_star           -1.0769  0.0311   1197.2  < 2e-16 ***
#   A0_PM2.5:A0_PM2.5_star  -0.6947  0.0909     58.4  2.2e-14 ***
unif.PNDE.qol <- coef(MSM.model.qol)["A0_PM2.5"]
# -5.38
unif.TNIE.qol <- (coef(MSM.model.qol)["A0_PM2.5_star"] +
                    coef(MSM.model.qol)["A0_PM2.5:A0_PM2.5_star"])
# -1.77

# 95% confidence intervals can also be computed by bootstrap

# ---------------------------------------------------------------------------- #
## 5.2) using medflex package ----
# ---------------------------------------------------------------------------- #
library(medflex)
# this package can used to estimate Natural direct and indirect effects when they are
# identifiable

# ---------------------------------------------------------------------------- #
### 5.2.1) estimate PNDE and TNDE, conditional on L0 and L1 ----
# ---------------------------------------------------------------------------- #
## Weighting based approach of medflex
# rm(list=ls())
# get back to the original data set
df1_int <- subset(df1_int, select = -c(A0_PM2.5_temp, id))

## 1) Expand the dataset and calculate weights using the neWeight() function
g.M <- glm(M_diabetes ~ factor(A0_PM2.5) + # the 1st variable should be the exposure
             L0_male + L0_soc_env + L1, # then add baseline confounders
           family = "binomial", data = df1_int)

# expend the data set and add a second column A* = 1 - A
# and calculate the weights = sw_M.Astar / sw_M.A
exp.Data <- neWeight(g.M)

# it is also possible to use the neWeight function directly:
# exp.Data <- neWeight(M_diabetes ~ factor(A0_PM2.5) + # the 1st variable should be the exposure
#                        L0_male + L0_soc_env + L1, # then add baseline confounders
#                      family = "binomial", data = df1_int)

class(exp.Data)
# [1] "data.frame" "expData"    "weightData"
head(exp.Data)
#   id A0_PM2.50 A0_PM2.51 L0_male L0_soc_env L1 M_diabetes Y_death Y_qol
# 1  1         0         0       0          1  1          0       0  93.4
# 2  1         0         1       0          1  1          0       0  93.4
# 3  2         1         1       1          1  1          0       1  64.0
# 4  2         1         0       1          1  1          0       1  64.0
# 5  3         0         0       1          1  0          0       0  75.6
# 6  3         0         1       1          1  0          0       0  75.6

# the new variables A0_PM2.50 and A0_PM2.51 correspond to A and A*

# check the weights:
w <- weights(exp.Data)
boxplot(w)

# Note that we estimated the same weights by hand previously:
head(data.frame(w, sw_M.Astar / sw_M.A), 6)
#       w sw_M.Astar.sw_M.A
# 1 1.000             1.000
# 2 0.801             0.801
# 3 1.000             1.000
# 4 1.293             1.293
# 5 1.000             1.000
# 6 0.809             0.809
plot(w, sw_M.Astar / sw_M.A)
abline(0,1)

## 2) Fit the natural effect model using the neModel() function
##    (adjusted for baseline confounders)

## for death
set.seed(1234)
neMod.death <- neModel(Y_death ~ A0_PM2.50 * A0_PM2.51 + L0_male + L0_soc_env + L1,
                       family = "gaussian", # to estimate risk difference
                       expData = exp.Data,
                       se = c("bootstrap"), # or "robust"
                       nBoot = 10, # use >= 1000 samples, if se = bootstrap
                       )
summary(neMod.death)
#                         Estimate Std. Error z value Pr(>|z|)
#   (Intercept)            0.10982    0.00807   13.61  < 2e-16 ***
#   A0_PM2.501             0.06576    0.01642    4.01  6.2e-05 ***
#   A0_PM2.511             0.00836    0.00214    3.91  9.1e-05 ***
#   L0_male                0.05041    0.00741    6.80  1.0e-11 ***
#   L0_soc_env             0.06113    0.00609   10.03  < 2e-16 ***
#   L1                     0.08341    0.01281    6.51  7.5e-11 ***
#   A0_PM2.501:A0_PM2.511  0.00238    0.00446    0.53     0.59

## estimate the conditional PNDE and TNDE, given L(0)=0 and L(1)=0
medflex.PNDE.death.L0 <- coef(neMod.death)["A0_PM2.501"]
# 0.0658
medflex.TNIE.death.L0 <- (coef(neMod.death)["A0_PM2.511"] +
                            coef(neMod.death)["A0_PM2.501:A0_PM2.511"])
# 0.0107

# 95% CI of the TNIE (which combines 2 coefficients) ?
# cov(b1,b2) = var(b1) + var(b2) + 2 * cov(b1,b2)
c(medflex.TNIE.death.L0 -
    qnorm(0.975) * sqrt(var(neMod.death$bootRes$t[,3]) + # column of A0_PM2.511
                          var(neMod.death$bootRes$t[,7]) + # column of A0_PM2.501:A0_PM2.511
                          2 * var(neMod.death$bootRes$t[,c(3,7)])[1,2]),
  medflex.TNIE.death.L0 +
    qnorm(0.975) * sqrt(var(neMod.death$bootRes$t[,3]) +
                          var(neMod.death$bootRes$t[,7]) +
                          2 * var(neMod.death$bootRes$t[,c(3,7)])[1,2]))
# A0_PM2.511 A0_PM2.511
#     0.0027     0.0188

## for quality of life
neMod.qol <- neModel(Y_qol ~ A0_PM2.50 * A0_PM2.51 + L0_male + L0_soc_env + L1,
                     family = "gaussian", # to estimate risk difference
                     expData = exp.Data,
                     se = c("robust") # or "bootstrap"
                     # nBoot = 1000, # if se = bootstrap
                     )
summary(neMod.qol)
#                         Estimate Std. Error z value Pr(>|z|)
#   (Intercept)             73.187      0.226  324.09  < 2e-16 ***
#   A0_PM2.501              -5.336      0.338  -15.77  < 2e-16 ***
#   A0_PM2.511              -1.069      0.136   -7.87  3.4e-15 ***
#   L0_male                 -1.231      0.223   -5.53  3.3e-08 ***
#   L0_soc_env              -3.546      0.230  -15.41  < 2e-16 ***
#   L1                      -4.101      0.243  -16.87  < 2e-16 ***
#   A0_PM2.501:A0_PM2.511   -0.765      0.126   -6.09  1.1e-09 ***

## estimate the conditional PNDE and TNDE, given L(0)=0 and L(1)=0
medflex.PNDE.qol.L0 <- coef(neMod.qol)["A0_PM2.501"]
# -5.34
medflex.TNIE.qol.L0 <- (coef(neMod.qol)["A0_PM2.511"] +
                            coef(neMod.qol)["A0_PM2.501:A0_PM2.511"])
# -1.83

# 95% CI of the TNIE (which combines 2 coefficients) ?
# cov(b1,b2) = var(b1) + var(b2) + 2cov(b1,b2)
c(medflex.TNIE.qol.L0 -
    qnorm(0.975) * sqrt(neMod.qol$vcov["A0_PM2.511","A0_PM2.511"] +
                          neMod.qol$vcov["A0_PM2.501:A0_PM2.511","A0_PM2.501:A0_PM2.511"] +
                          2 * neMod.qol$vcov["A0_PM2.511","A0_PM2.501:A0_PM2.511"]),
  medflex.TNIE.qol.L0 +
    qnorm(0.975) * sqrt(neMod.qol$vcov["A0_PM2.511","A0_PM2.511"] +
                          neMod.qol$vcov["A0_PM2.501:A0_PM2.511","A0_PM2.501:A0_PM2.511"] +
                          2 * neMod.qol$vcov["A0_PM2.511","A0_PM2.501:A0_PM2.511"]))
# A0_PM2.511 A0_PM2.511
#     -2.30      -1.37

# plot results
plot(neMod.qol)

# ---------------------------------------------------------------------------- #
### 5.2.2) estimate marginal PNDE and TNDE (population average) ----
# ---------------------------------------------------------------------------- #

## 1) Expand the dataset and calculate weights using the neWeight() function

## Estimate a model of the exposure
g.A.L0 <- glm(A0_PM2.5 ~ L0_male + L0_soc_env + L1, # will no work without L(1) !
              family = "binomial", data = df1_int)

## Estimate a model of the mediator
g.M <- glm(M_diabetes ~ factor(A0_PM2.5) + # the 1st variable should be the exposure
             L0_male + L0_soc_env + L1, # then add baseline confounders
           family = "binomial", data = df1_int)

# expend the data set and add a second column A* = 1 - A
# and calculate the weights = sw_M.Astar / sw_M.A
exp.Data <- neWeight(g.M)

## 2) Fit the natural effect model using the neModel() function
## for death
neMod.death.pop <- neModel(Y_death ~ A0_PM2.50 * A0_PM2.51, # no need for L(0),L(1)
                       family = "gaussian", # to estimate risk difference
                       expData = exp.Data,
                       xFit = g.A.L0, # model of the exposure
                       se = c("robust")) # or "bootstrap"
summary(neMod.death.pop)
# Parameter estimates:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)            0.198840   0.004247   46.82  < 2e-16 ***
#   A0_PM2.501             0.065950   0.014346    4.60  4.3e-06 ***
#   A0_PM2.511             0.008534   0.001619    5.27  1.4e-07 ***
#   A0_PM2.501:A0_PM2.511 -0.000406   0.003744   -0.11     0.91

## estimate the population average PNDE and TNDE,
medflex.PNDE.death <- coef(neMod.death.pop)["A0_PM2.501"]
# 0.066
medflex.TNIE.death <- (coef(neMod.death.pop)["A0_PM2.511"] +
                          coef(neMod.death.pop)["A0_PM2.501:A0_PM2.511"])
# 0.00813
# the results are the same than results calculated by hand

## for quality of life
neMod.qol.pop <- neModel(Y_qol ~ A0_PM2.50 * A0_PM2.51, # no need for L(0),L(1)
                         family = "gaussian", # to estimate risk difference
                         expData = exp.Data,
                         xFit = g.A.L0, # model of the exposure
                         se = c("robust")) # or "bootstrap"
summary(neMod.qol.pop)
# Parameter estimates:
#                         Estimate Std. Error z value Pr(>|z|)
#   (Intercept)             69.085      0.117  590.36  < 2e-16 ***
#   A0_PM2.501              -5.377      0.350  -15.38  < 2e-16 ***
#   A0_PM2.511              -1.077      0.137   -7.88  3.3e-15 ***
#   A0_PM2.501:A0_PM2.511   -0.695      0.120   -5.80  6.7e-09 ***

## estimate the population average PNDE and TNDE,
medflex.PNDE.qol <- coef(neMod.qol.pop)["A0_PM2.501"]
# -5.38
medflex.TNIE.qol <- (coef(neMod.qol.pop)["A0_PM2.511"] +
                       coef(neMod.qol.pop)["A0_PM2.501:A0_PM2.511"])
# -1.77
# the results are the same than results calculated by hand




# ---------------------------------------------------------------------------- #
# 6) Unified approach based on MSM, estimated by g-computation ----
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
## 6.1) manual calculation ----
# ---------------------------------------------------------------------------- #
# from Theis Lange, Stijn Vansteelandt, Maarten Bekaert, A Simple Unified Approach
# for Estimating Natural Direct and Indirect Effects. Am J Epidemiol 2012;176,(3)190–195.
# https://doi-org.proxy.insermbiblio.inist.fr/10.1093/aje/kwr525
rm(list=ls())
# this method is used to estimate Natural direct and indirect effects when they are
# identifiable
df1_int <- read.csv(file = "data/df1_int.csv")

## Step 1. Estimate a suitable model for the outcome conditional on the exposure,
##         mediator and confounders, using the original data set
Q.Y.death <- glm(Y_death ~ L0_male + L0_soc_env + L1  + A0_PM2.5 * M_diabetes,
                 family = "binomial",
                 data = df1_int)

Q.Y.qol <- glm(Y_qol ~ L0_male + L0_soc_env + L1 + A0_PM2.5 * M_diabetes,
               family = "gaussian",
               data = df1_int)

## Step 2. Expand data and impute counterfactuals E(Y_{a,M_a*}|L(0),L(1))
## Expand data
# create identifier
df1_int$id <- 1:nrow(df1_int)

# dupplicate the original data set
df1_int.1 <- df1_int.2 <- df1_int

# Add an A* variable where A* = A
df1_int.1$A0_PM2.5_star <- df1_int$A0_PM2.5
df1_int.2$A0_PM2.5_star <- df1_int$A0_PM2.5

# for the 2nd data set, change the exposure A = 1 - A
df1_int.2$A0_PM2.5 <- 1 - df1_int$A0_PM2.5

# stack the 2 new datasets
df1_int.new <- rbind(df1_int.1,df1_int.2)
# sort the data set by the id-variable to use geeglm
# (necessary to use geepack and the geeglm() later)
df1_int.new <- df1_int.new[order(df1_int.new$id), ]

head(df1_int.new, 12)
#       L0_male L0_soc_env A0_PM2.5 L1 M_diabetes Y_death Y_qol A0_PM2.5_temp id A0_PM2.5_star
# 1           0          1        0  1          0       0  93.4             0  1             0
# 10001       0          1        1  1          0       0  93.4             0  1             0
# 2           1          1        1  1          0       1  64.0             1  2             1
# 10002       1          1        0  1          0       1  64.0             1  2             1
# 3           1          1        0  0          0       0  75.6             0  3             0
# 10003       1          1        1  0          0       0  75.6             0  3             0
# 4           1          0        0  0          0       0  89.8             0  4             0
# 10004       1          0        1  0          0       0  89.8             0  4             0
# 5           1          1        0  0          0       0  77.2             0  5             0
# 10005       1          1        1  0          0       0  77.2             0  5             0
# 6           1          1        0  0          1       0  73.9             0  6             0
# 10006       1          1        1  0          1       0  73.9             0  6             0

## Impute counterfactuals
# death - we predict the probability P(Y=1|A,M,L(0),L(1))
Y.pred.death <- predict(Q.Y.death, newdata = df1_int.new, type = "response")

# qol - we can predict the mean E(Y|A,M,L(0),L(1))
Y.pred.qol <- predict(Q.Y.qol, newdata = df1_int.new, type = "response")

## Step 3. Fit a suitable model to the outcome incluing only A and A* (and their
##         interaction)
library(geepack)
# for death (conditional)
MSM.model.death.1 <- geeglm(Y.pred.death ~ A0_PM2.5 * A0_PM2.5_star +
                              L0_male + L0_soc_env + L1,
                            family = "gaussian",
                            id = df1_int.new$id,
                            data = df1_int.new,
                            scale.fix = TRUE)
summary(MSM.model.death.1)
#                          Estimate  Std.err   Wald Pr(>|W|)
#   (Intercept)            0.104153 0.000611  29104  < 2e-16 ***
#   A0_PM2.5               0.063543 0.000138 210638  < 2e-16 ***
#   A0_PM2.5_star          0.007032 0.001112     40  2.5e-10 ***
#   L0_male                0.054291 0.000706   5922  < 2e-16 ***
#   L0_soc_env             0.065692 0.000692   9004  < 2e-16 ***
#   L1                     0.086446 0.000860  10100  < 2e-16 ***
#   A0_PM2.5:A0_PM2.5_star 0.004918 0.000397    154  < 2e-16 ***

## estimate the conditional PNDE and TNDE, given L(0)=0 and L(1)=0
PNDE.death.L0 <- coef(MSM.model.death.1)["A0_PM2.5"]
# 0.0635
TNIE.death.L0 <- (coef(MSM.model.death.1)["A0_PM2.5_star"] +
                            coef(MSM.model.death.1)["A0_PM2.5:A0_PM2.5_star"])
# 0.012

# For quality of life (conditional)
MSM.model.qol.1 <- geeglm(Y.pred.qol ~ A0_PM2.5 * A0_PM2.5_star +
                              L0_male + L0_soc_env + L1,
                            family = "gaussian",
                            id = df1_int.new$id,
                            data = df1_int.new,
                            scale.fix = TRUE)
summary(MSM.model.qol.1)
#                          Estimate Std.err     Wald Pr(>|W|)
#   (Intercept)             73.3371  0.0949 597466.3  < 2e-16 ***
#   A0_PM2.5                -5.3013  0.0268  39033.2  < 2e-16 ***
#   A0_PM2.5_star           -1.0535  0.1353     60.6  6.9e-15 ***
#   L0_male                 -1.3289  0.1038    164.0  < 2e-16 ***
#   L0_soc_env              -3.6467  0.1061   1180.6  < 2e-16 ***
#   L1                      -4.2325  0.1168   1312.8  < 2e-16 ***
#   A0_PM2.5:A0_PM2.5_star  -0.7920  0.0870     82.9  < 2e-16 ***

## estimate the conditional PNDE and TNDE, given L(0)=0 and L(1)=0
PNDE.qol.L0 <- coef(MSM.model.qol.1)["A0_PM2.5"]
# -5.3
TNIE.qol.L0 <- (coef(MSM.model.qol.1)["A0_PM2.5_star"] +
                  coef(MSM.model.qol.1)["A0_PM2.5:A0_PM2.5_star"])
# -1.85

## For population average estimations
## Estimate a model of the exposure
g.A.L0 <- glm(A0_PM2.5 ~ L0_male + L0_soc_env + L1, # will no work without L(1) !
              family = "binomial", data = df1_int)
g.Ais1.L0 <- predict(g.A.L0, newdata = df1_int.new, type = "response")

w <- rep(NA, nrow(df1_int.new))
w[df1_int.new$A0_PM2.5_star == 1] <- 1 / g.Ais1.L0[df1_int.new$A0_PM2.5_star == 1]
w[df1_int.new$A0_PM2.5_star == 0] <- 1 / (1 - g.Ais1.L0[df1_int.new$A0_PM2.5_star == 0])

# for death (population average)
MSM.model.death.2 <- geeglm(Y.pred.death ~ A0_PM2.5 * A0_PM2.5_star,
                            family = "gaussian",
                            weights = w,
                            id = df1_int.new$id,
                            data = df1_int.new,
                            scale.fix = TRUE)
summary(MSM.model.death.2)
#                          Estimate  Std.err    Wald Pr(>|W|)
#   (Intercept)            0.198882 0.000644  95250.3  < 2e-16 ***
#   A0_PM2.5               0.063867 0.000138 213552.3  < 2e-16 ***
#   A0_PM2.5_star          0.009213 0.002009     21.0  4.5e-06 ***
#   A0_PM2.5:A0_PM2.5_star 0.002330 0.000445     27.4  1.7e-07 ***

## estimate the conditional PNDE and TNDE, given L(0)=0 and L(1)=0
PNDE.death <- coef(MSM.model.death.2)["A0_PM2.5"]
# 0.0639
TNIE.death <- (coef(MSM.model.death.2)["A0_PM2.5_star"] +
                 coef(MSM.model.death.2)["A0_PM2.5:A0_PM2.5_star"])
# 0.0115

# For quality of life (conditional)
MSM.model.qol.2 <- geeglm(Y.pred.qol ~ A0_PM2.5 * A0_PM2.5_star,
                          family = "gaussian",
                          weights = w,
                          id = df1_int.new$id,
                          data = df1_int.new,
                          scale.fix = TRUE)
summary(MSM.model.qol.2)
#                          Estimate Std.err     Wald Pr(>|W|)
#   (Intercept)             69.0849  0.0492 1.97e+06  < 2e-16 ***
#   A0_PM2.5                -5.3108  0.0269 3.88e+04  < 2e-16 ***
#   A0_PM2.5_star           -1.1690  0.1624 5.18e+01  6.1e-13 ***
#   A0_PM2.5:A0_PM2.5_star  -0.7395  0.0906 6.66e+01  3.3e-16 ***

## estimate the conditional PNDE and TNDE, given L(0)=0 and L(1)=0
PNDE.qol.L0 <- coef(MSM.model.qol.2)["A0_PM2.5"]
# -5.31
TNIE.qol.L0 <- (coef(MSM.model.qol.2)["A0_PM2.5_star"] +
                  coef(MSM.model.qol.2)["A0_PM2.5:A0_PM2.5_star"])
# -1.91


## Other solution: estimate coefficients of the MSM directly for population level
## estimates of PNDE and TNIE
MSM.model.death.3 <- geeglm(Y.pred.death ~ A0_PM2.5 * A0_PM2.5_star,
                          family = "gaussian",
                          id = df1_int.new$id,
                          data = df1_int.new,
                          scale.fix = TRUE)
summary(MSM.model.death.3)
# Coefficients:
#                          Estimate  Std.err   Wald Pr(>|W|)
#   (Intercept)            0.197386 0.000639  95446   <2e-16 ***
#   A0_PM2.5               0.063543 0.000138 210638   <2e-16 ***
#   A0_PM2.5_star          0.021519 0.001923    125   <2e-16 ***
#   A0_PM2.5:A0_PM2.5_star 0.004918 0.000397    154   <2e-16 ***
unif.PNDE.death <- coef(MSM.model.death.3)["A0_PM2.5"]
# 0.0635
unif.TNIE.death <- (coef(MSM.model.death.3)["A0_PM2.5_star"] +
                      coef(MSM.model.death.3)["A0_PM2.5:A0_PM2.5_star"])
# 0.0264

## For Quality of life
MSM.model.qol.3 <- geeglm(Y.pred.qol ~ A0_PM2.5 * A0_PM2.5_star,
                          family = "gaussian",
                          id = df1_int.new$id,
                          data = df1_int.new,
                          scale.fix = TRUE)
summary(MSM.model.qol.3)
# Coefficients:
#                          Estimate Std.err     Wald Pr(>|W|)
#   (Intercept)             69.1529  0.0491 1.98e+06   <2e-16 ***
#   A0_PM2.5                -5.3013  0.0268 3.90e+04   <2e-16 ***
#   A0_PM2.5_star           -1.6894  0.1516 1.24e+02   <2e-16 ***
#   A0_PM2.5:A0_PM2.5_star  -0.7920  0.0870 8.29e+01   <2e-16 ***
unif.PNDE.qol <- coef(MSM.model.qol.3)["A0_PM2.5"]
# -5.3
unif.TNIE.qol <- (coef(MSM.model.qol.3)["A0_PM2.5_star"] +
                    coef(MSM.model.qol.3)["A0_PM2.5:A0_PM2.5_star"])
# -2.48


# ---------------------------------------------------------------------------- #
## 6.2) using medflex package ----
# ---------------------------------------------------------------------------- #
library(medflex)
# this package can used to estimate Natural direct and indirect effects when they are
# identifiable

# ---------------------------------------------------------------------------- #
### 6.2.1) estimate PNDE and TNDE, conditional on L0 and L1 ----
# ---------------------------------------------------------------------------- #
## Weighting based approach of medflex
# rm(list=ls())
# get back to the original data set
df1_int <- subset(df1_int, select = -c(id))
df1_int$A0_PM2.5 <- factor(df1_int$A0_PM2.5)

## 1) Fit a model for the outcome
Q.Y.death <- glm(Y_death ~ A0_PM2.5 + # start with the exposure
                   M_diabetes + A0_PM2.5:M_diabetes + # then the mediator
                   L0_male + L0_soc_env + L1, # then baseline confounders
                 family = "binomial",
                 data = df1_int)

Q.Y.qol <- glm(Y_qol ~ A0_PM2.5 + # start with the exposure
                 M_diabetes + A0_PM2.5:M_diabetes + # then the mediator
                 L0_male + L0_soc_env + L1, # then baseline confounders
               family = "gaussian",
               data = df1_int)


## 2) expend the data set and add a second column A* = 1 - A
exp.Data.death <- neImpute(Q.Y.death)
head(exp.Data.death)
#   id A0_PM2.50 A0_PM2.51 L0_male L0_soc_env L1 M_diabetes Y_death Y_qol
# 1  1         0         0       0          1  1          0   0.222  93.4
# 2  1         1         0       0          1  1          0   0.292  93.4
# 3  2         1         1       1          1  1          0   0.356  64.0
# 4  2         0         1       1          1  1          0   0.277  64.0
# 5  3         0         0       1          1  0          0   0.197  75.6
# 6  3         1         0       1          1  0          0   0.261  75.6

exp.Data.qol <- neImpute(Q.Y.qol)
head(exp.Data.qol)
#   id A0_PM2.50 A0_PM2.51 L0_male L0_soc_env L1 M_diabetes Y_death Y_qol
# 1  1         0         0       0          1  1          0       0  68.4
# 2  1         1         0       0          1  1          0       0  64.7
# 3  2         1         1       1          1  1          0       1  64.0
# 4  2         0         1       1          1  1          0       1  67.7
# 5  3         0         0       1          1  0          0       0  71.2
# 6  3         1         0       1          1  0          0       0  67.4

# the predicted outcomes by medflex are the same than our previous predictions
plot(exp.Data.death$Y_death, Y.pred.death)
abline(0,1)

plot(exp.Data.qol$Y_qol, Y.pred.qol)
abline(0,1)

## 3) Fit the natural effect model using the neModel() function
##    (adjusted for baseline confounders)

## for death
neMod.death <- neModel(Y_death ~ A0_PM2.50 * A0_PM2.51 + L0_male + L0_soc_env + L1,
                       family = "gaussian", # to estimate risk difference
                       expData = exp.Data.death,
                       se = c("robust") # or "bootstrap"
                       )
summary(neMod.death)
#                         Estimate Std. Error z value Pr(>|z|)
#   (Intercept)            0.10415    0.00829   12.56  < 2e-16 ***
#   A0_PM2.501             0.06354    0.01373    4.63  3.7e-06 ***
#   A0_PM2.511             0.00703    0.00175    4.01  6.1e-05 ***
#   L0_male                0.05429    0.00867    6.26  3.8e-10 ***
#   L0_soc_env             0.06569    0.00882    7.45  9.3e-14 ***
#   L1                     0.08645    0.00999    8.66  < 2e-16 ***
#   A0_PM2.501:A0_PM2.511  0.00492    0.00391    1.26     0.21

## estimate the conditional PNDE and TNDE, given L(0)=0 and L(1)=0
medflex.PNDE.death.L0 <- coef(neMod.death)["A0_PM2.501"]
# 0.0635
medflex.TNIE.death.L0 <- (coef(neMod.death)["A0_PM2.511"] +
                            coef(neMod.death)["A0_PM2.501:A0_PM2.511"])
# 0.012

## for quality of life
neMod.qol <- neModel(Y_qol ~ A0_PM2.50 * A0_PM2.51 + L0_male + L0_soc_env + L1,
                     family = "gaussian", # to estimate risk difference
                     expData = exp.Data.qol,
                     se = c("robust")) # or "bootstrap"
summary(neMod.qol)
#                         Estimate Std. Error z value Pr(>|z|)
#   (Intercept)             73.337      0.230  319.17  < 2e-16 ***
#   A0_PM2.501              -5.301      0.340  -15.61  < 2e-16 ***
#   A0_PM2.511              -1.054      0.139   -7.58  3.6e-14 ***
#   L0_male                 -1.329      0.227   -5.85  5.0e-09 ***
#   L0_soc_env              -3.647      0.235  -15.54  < 2e-16 ***
#   L1                      -4.233      0.249  -17.01  < 2e-16 ***
#   A0_PM2.501:A0_PM2.511   -0.792      0.127   -6.21  5.1e-10 ***

## estimate the conditional PNDE and TNDE, given L(0)=0 and L(1)=0
medflex.PNDE.qol.L0 <- coef(neMod.qol)["A0_PM2.501"]
# -5.3
medflex.TNIE.qol.L0 <- (coef(neMod.qol)["A0_PM2.511"] +
                          coef(neMod.qol)["A0_PM2.501:A0_PM2.511"])
# -1.85

# ---------------------------------------------------------------------------- #
### 6.2.2) estimate marginal PNDE and TNDE (population average) ----
# ---------------------------------------------------------------------------- #
## Estimate a model of the exposure
g.A.L0 <- glm(A0_PM2.5 ~ L0_male + L0_soc_env + L1, # will no work without L(1) !
              family = "binomial", data = df1_int)

## 3) Fit the natural effect model using the neModel() function
## for death
neMod.death.pop <- neModel(Y_death ~ A0_PM2.50 * A0_PM2.51, # no need for L(0),L(1)
                           family = "gaussian", # to estimate risk difference
                           expData = exp.Data.death,
                           xFit = g.A.L0, # model of the exposure
                           se = c("robust")) # or "bootstrap"
summary(neMod.death.pop)
# Parameter estimates:
#                         Estimate Std. Error z value Pr(>|z|)
#   (Intercept)            0.19888    0.00425   46.83  < 2e-16 ***
#   A0_PM2.501             0.06387    0.01378    4.64  3.6e-06 ***
#   A0_PM2.511             0.00921    0.00169    5.45  5.1e-08 ***
#   A0_PM2.501:A0_PM2.511  0.00233    0.00367    0.64     0.53

## estimate the population average PNDE and TNDE,
medflex.PNDE.death <- coef(neMod.death.pop)["A0_PM2.501"]
# 0.0639
medflex.TNIE.death <- (coef(neMod.death.pop)["A0_PM2.511"] +
                         coef(neMod.death.pop)["A0_PM2.501:A0_PM2.511"])
# 0.0115
# the results are the same than results calculated by hand

## for quality of life
neMod.qol.pop <- neModel(Y_qol ~ A0_PM2.50 * A0_PM2.51, # no need for L(0),L(1)
                         family = "gaussian", # to estimate risk difference
                         expData = exp.Data.qol,
                         xFit = g.A.L0, # model of the exposure
                         se = c("robust")) # or "bootstrap"
summary(neMod.qol.pop)
# Parameter estimates:
#                         Estimate Std. Error z value Pr(>|z|)
#   (Intercept)             69.085      0.117  590.50  < 2e-16 ***
#   A0_PM2.501              -5.311      0.339  -15.65  < 2e-16 ***
#   A0_PM2.511              -1.169      0.145   -8.09  6.1e-16 ***
#   A0_PM2.501:A0_PM2.511   -0.740      0.125   -5.92  3.2e-09 ***

## estimate the population average PNDE and TNDE,
medflex.PNDE.qol <- coef(neMod.qol.pop)["A0_PM2.501"]
# -5.31
medflex.TNIE.qol <- (coef(neMod.qol.pop)["A0_PM2.511"] +
                       coef(neMod.qol.pop)["A0_PM2.501:A0_PM2.511"])
# -1.91
# the results are the same than results calculated by hand





