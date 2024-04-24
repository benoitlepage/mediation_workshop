### test programs for estimations by MSMs

rm(list=ls())

df1 <- read.csv(file = "data/df1.csv")
df1_int <- read.csv(file = "data/df1_int.csv")
df2 <- read.csv(file = "data/df2.csv")
df2_int <- read.csv(file = "data/df2_int.csv")

################################################################################
######################### Estimation of the Average Total Effect (ATE)
################################################################################

rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")
### MSM of ATE|L_male, estimated by IPTW ---------------------------------------

## 1. Denominator of the weight
# 1a. Estimate g(A=a_i|L(0)) (denominator of the weight)
g.A.L <- glm(A0_ace ~ L0_male + L0_parent_low_educ_lv,
           family = "binomial", data = df1_int)

# 1b. Predict each individual's probability of being exposed to her own exposure
# predict the probabilities P(A0_ace=1 & P(A0_ace=0)
pred.g1.L <- predict(g.A.L, type="response")
pred.g0.L <- 1 - pred.g1.L
# the predicted probability of the observed treatment P(A = a_i | L(0)) is :
gAi.L <- rep(NA, nrow(df1_int))
gAi.L[df1_int$A0_ace==1] <- pred.g1.L[df1_int$A0_ace==1]
gAi.L[df1_int$A0_ace==0] <- pred.g0.L[df1_int$A0_ace==0]

## 2. Numerator of the weight
# The numerator of the weight can be 1 for simple weights,
# or g(A=a_i|V) to obtain stabilized weights which put less weight to individuals
# with less observation. Stabilized weights enable a weaker positivity assumption.

# 2a. Estimate g(A=a_i | sex) (numerator of the stabilized weight)
g.A.sex <- glm(A0_ace ~ L0_male,
           family = "binomial", data = df1_int)

# 2b. Predict each individual's probability of being exposed to her own exposure
# predict the probabilities P(A0_ace=1 | sex) & P(A0_ace=0 | sex)
pred.g1.sex <- predict(g.A.sex, type="response")
pred.g0.sex <- 1 - pred.g1.sex
# the predicted probability of the observed treatment P(A = a_i | sex) is :
gAi.sex <- rep(NA, nrow(df1_int))
gAi.sex[df1_int$A0_ace==1] <- pred.g1.sex[df1_int$A0_ace==1]
gAi.sex[df1_int$A0_ace==0] <- pred.g0.sex[df1_int$A0_ace==0]


## 3. Define individual weights:
# We can use simple weights w = 1 / g(A=a_i | L(0))
w <- 1 / gAi.L

# Or alternatively, we can use stabilized weights : sw = g(A=a_i) / g(A=a_i | L(0))
sw <- gAi.sex / gAi.L

par(mfcol = c(1,2))
boxplot(w ~ df1_int$A0_ace)
boxplot(sw ~ df1_int$A0_ace)


## 4. Estimate coefficients of the MSM using a weighted regression E(Y | A, sex)
# a GLM with gaussian family can be applied to estimate risk difference
msm1 <- glm(Y_death ~ A0_ace + L0_male + A0_ace*L0_male,
           weights = w,
           family = "gaussian",
           data = df1_int)
coef(msm1)
# (Intercept)         A0_ace        L0_male A0_ace:L0_male
# 0.17573472     0.03589627     0.04598911     0.04136896

msm2 <- glm(Y_death ~ A0_ace + L0_male + A0_ace*L0_male,
            weights = sw,
            family = "gaussian",
            data = df1_int)
coef(msm2)
# (Intercept)         A0_ace        L0_male A0_ace:L0_male
# 0.17573472     0.03589627     0.04598911     0.04136896

## 5. Estimate the ATE stratified by sex
# According to MSM1 (with simple weights)
ATE.msm1.male0 <- coef(msm1)["A0_ace"]
# 0.03589627
ATE.msm1.male1 <- coef(msm1)["A0_ace"] + coef(msm1)["A0_ace:L0_male"]
# 0.07726522

# According to MSM2 (with stabilized weights)
# According to MSM1 (with simple weights)
ATE.msm2.male0 <- coef(msm2)["A0_ace"]
# 0.03589627
ATE.msm2.male1 <- coef(msm2)["A0_ace"] + coef(msm2)["A0_ace:L0_male"]
# 0.07726522
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
  g.A.L <- glm(A0_ace ~ L0_male + L0_parent_low_educ_lv,
               family = "binomial", data = bootData)

  # 1b. Predict each individual's probability of being exposed to her own exposure
  # the predicted probability of the observed treatment P(A = a_i | L(0)) is :
  gAi.L <- rep(NA, nrow(bootData))
  gAi.L[bootData$A0_ace==1] <- predict(g.A.L, type="response")[bootData$A0_ace==1]
  gAi.L[bootData$A0_ace==0] <- (1 - predict(g.A.L, type="response"))[bootData$A0_ace==0]

  # 2a. Estimate g(A=a_i | sex) (numerator of the stabilized weight)
  g.A.sex <- glm(A0_ace ~ L0_male,
                 family = "binomial", data = bootData)

  # 2b. Predict each individual's probability of being exposed to her own exposure
  # the predicted probability of the observed treatment P(A = a_i | sex) is :
  gAi.sex <- rep(NA, nrow(bootData))
  gAi.sex[bootData$A0_ace==1] <- predict(g.A.sex, type="response")[bootData$A0_ace==1]
  gAi.sex[bootData$A0_ace==0] <- (1 - predict(g.A.sex, type="response"))[bootData$A0_ace==0]

  # 3. Define individual weights:
  w <- 1 / gAi.L
  sw <- gAi.sex / gAi.L

  # 4. Estimate coefficients of the MSM using a weighted regression E(Y | A, sex)
  msm1 <- glm(Y_death ~ A0_ace + L0_male + A0_ace*L0_male,
              weights = w,
              family = "gaussian",
              data = bootData)

  msm2 <- glm(Y_death ~ A0_ace + L0_male + A0_ace*L0_male,
              weights = sw,
              family = "gaussian",
              data = bootData)

  ## 5. Estimate the ATE stratified by sex
  # Using coefficients of msm1
  bootstrap.estimates[b,"boot.msm1.male0"] <- coef(msm1)["A0_ace"]
  bootstrap.estimates[b,"boot.msm1.male1"] <- (coef(msm1)["A0_ace"] +
                                                 coef(msm1)["A0_ace:L0_male"])

  # Using coefficients of msm2
  bootstrap.estimates[b,"boot.msm2.male0"] <- coef(msm2)["A0_ace"]
  bootstrap.estimates[b,"boot.msm2.male1"] <- (coef(msm2)["A0_ace"] +
                                                 coef(msm2)["A0_ace:L0_male"])
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
  p_L0_parent_low_educ_lv <- 0.65

  # A: A0_ace <- rbinom( 0.05 + 0.04 * L0_male + 0.06 * L0_parent_low_educ_lv )
  b_A <- 0.05   # reference prevalence is 5%
  b_male_A <- 0.04  # + 0.04 for the effect of L0_male -> A0_ace
  b_parent_educ_A <- 0.06  # +0.06 for the effect of L0_parent_low_educ_lv -> A0_ace

  # L1: intermediate confounder between M and Y, not influenced by A
  p_L1 <- 0.3

  # M: M_smoking <- rbinom( 0.2 + 0.05 * L0_male + 0.06 * L0_parent_low_educ_lv + 0.07 * L1 +
  #                         0.1 * A0_ace )
  b_M <- 0.2 # reference prevalence is 20%
  b_male_M <- 0.05 # +0.05 for the effect of L0_male -> M_smoking
  b_parent_educ_M <- 0.06 # +0.06 for the effect of L0_parent_low_educ_lv -> M_smoking
  b_L1_M <- 0.07 # +0.07 for the effect of L1 -> M_smoking
  b_A_M <- 0.1 # +0.10 for the effect of A0_ace -> M_smoking

  # Y binary: rbinom( 0.10 + 0.06 * L0_male + 0.04 * L0_parent_low_educ_lv + 0.05 * A0_ace +
  #                   0.07 * L1 + 0.08 * M_smoking +
  #                   0.03 * A0_ace * M_smoking * A.M.inter )
  b_Y <- 0.1 # reference prevalence is 10%
  b_male_Y <- 0.06 # +0.06 for the effect of L0_male -> Y
  b_parent_educ_Y <- 0.04 # +0.04 for the effect of L0_parent_low_educ_lv -> Y
  b_A_Y <- 0.05 # 0.05 for the effect of A0_ace -> Y
  b_L1_Y <- 0.07 # +0.07 for the effect of L1 -> Y
  b_M_Y <- 0.08 # 0.08 for the effect of M_smoking -> Y
  b_AM_Y <- 0.03 # 0.03 for the interaction effect A0_ace * M_smoking -> Y

  # Y continuous: (75 - 1 * L0_male - 3 * L0_parent_low_educ_lv - 4 * A0_ace -3.5 * L1 -
  #                9 * M_smoking -5 * A0_ace * M_smoking * A.M.inter ) +
  #                rnorm(N, mean = 0, sd = 10)
  mu_Y <- 75 # reference mean for QoL
  c_male_Y <- -1 # -1 for the effect of L0_male -> Y
  c_parent_educ_Y <- -3 # -3 for the effect of L0_parent_low_educ_lv -> Y
  c_A_Y <- -4 # -4 for the effect of A0_ace -> Y
  c_L1_Y <- -3.5 # -3.5 for the effect of L1 -> Y
  c_M_Y <- -9 # -9 for the effect of M_smoking -> Y
  c_AM_Y <- -5  # - 5 for the interaction effect A0_ace * M_smoking  -> Y
  sd_Y <- 10 # standard deviation of the residuals

  # A*M interaction ?
  A.M.inter <- A.M.interaction

  coef <- c( p_L0_male = p_L0_male, p_L0_parent_low_educ_lv = p_L0_parent_low_educ_lv,
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
      ((b["p_L0_parent_low_educ_lv"])^(S1[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S1[n,"parent_educ"]))

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
      ((b["p_L0_parent_low_educ_lv"])^(S0[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S0[n,"parent_educ"]))
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
Q.tot.death <- glm(Y_death ~ A0_ace + L0_male + L0_parent_low_educ_lv,
                   family = "gaussian", data = df1_int)
# The final result would be sligthly different if we applied a binomial family
# The Gaussian family corresponds to the true generating model in this example.

## 2. Predict an outcome for each subject, setting A=0 and A=1
# prepare data sets used to predict the outcome under the counterfactual
# scenarios setting A=0 and A=1
data.A1 <- data.A0 <- df1_int
data.A1$A0_ace <- 1
data.A0$A0_ace <- 0

# predict values under the same name in the corresponding counterfactual dataset
data.A1$Ya.death.pred <- predict(Q.tot.death, newdata = data.A1, type = "response")
data.A0$Ya.death.pred <- predict(Q.tot.death, newdata = data.A0, type = "response")

## 3. Append both counterfactual datasets in a single dataset
# number of row is twice the initial value (we have 2 counterfactual scenarios)
data.2scenarios <- rbind(data.A0, data.A1)

## 4. fit the MSM: E(Y_a|sex)
MSM.ATE.gcomp <- glm(Ya.death.pred ~ A0_ace + L0_male + A0_ace:L0_male,
                     family = "gaussian",
                     data = data.2scenarios)
coef(MSM.ATE.gcomp)
#  (Intercept)         A0_ace        L0_male A0_ace:L0_male
# 1.743994e-01   7.720726e-02   4.874750e-02   4.802530e-16

## 5. Estimate the ATE stratified by sex
# According to MSM.ATE.gcomp
ATE.MSM.gcomp.male0 <- coef(MSM.ATE.gcomp)["A0_ace"]
# 0.07720726
ATE.MSM.gcomp.male1 <- (coef(MSM.ATE.gcomp)["A0_ace"] +
                          coef(MSM.ATE.gcomp)["A0_ace:L0_male"])
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


rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")
### MSM of CDE, estimated by IPTW ----------------------------------------------
## 1. Stabilized weight for the exposure sw_{A,i}
# 1a. Estimate g(A=a_i|L(0)) (denominator of the weight)
g.A.L <- glm(A0_ace ~ L0_male + L0_parent_low_educ_lv,
             family = "binomial", data = df1_int)
# 1b. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i | L(0)) is :
gAi.L <- rep(NA, nrow(df1_int))
gAi.L[df1_int$A0_ace==1] <- predict(g.A.L, type="response")[df1_int$A0_ace==1]
gAi.L[df1_int$A0_ace==0] <- (1 - predict(g.A.L, type="response"))[df1_int$A0_ace==0]

# 1c. Estimate g(A=a_i) (numerator of the weight)
g.A <- glm(A0_ace ~ 1, family = "binomial", data = df1_int)
# 1d. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i) is :
gAi <- rep(NA, nrow(df1_int))
gAi[df1_int$A0_ace==1] <- predict(g.A, type="response")[df1_int$A0_ace==1]
gAi[df1_int$A0_ace==0] <- (1 - predict(g.A, type="response"))[df1_int$A0_ace==0]

# 1e. Calculate sw_{A,i}
sw_Ai <- gAi / gAi.L

## 2. Stabilized weight for the mediator sw_{M,i}
# 2a. Estimate g(M=m_i|L(0),A,L(1)) (denominator of the weight)
g.M.L <- glm(M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1,
             family = "binomial", data = df1_int)
# 2b. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(A = a_i | L(0)) is :
gMi.L <- rep(NA, nrow(df1_int))
gMi.L[df1_int$M_smoking==1] <- predict(g.M.L, type="response")[df1_int$M_smoking==1]
gMi.L[df1_int$M_smoking==0] <- (1 - predict(g.M.L, type="response"))[df1_int$M_smoking==0]

# 2c. Estimate g(M=m_i|A) (numerator of the weight)
g.M.A <- glm(M_smoking ~ A0_ace, family = "binomial", data = df1_int)
# 2d. Predict each individual's probability of being exposed to her own exposure
# the predicted probability of the observed treatment g(M = m_i|A) is :
gMi.A <- rep(NA, nrow(df1_int))
gMi.A[df1_int$M_smoking==1] <- predict(g.M.A, type="response")[df1_int$M_smoking==1]
gMi.A[df1_int$M_smoking==0] <- (1 - predict(g.M.A, type="response"))[df1_int$M_smoking==0]
# 2e. Calculate sw_{M,i}
sw_Mi <- gMi.A / gMi.L

## 3. Define the individual stabilized weight for the CDE_m
sw_cde <- sw_Ai * sw_Mi

## 4. Estimate coefficients of the MSM using a weighted regression E(Y | A, sex)
# a GLM with gaussian family can be applied to estimate risk difference
msm_cde <- glm(Y_death ~ A0_ace + M_smoking + A0_ace*M_smoking,
               weights = sw_cde,
               family = "gaussian",
               data = df1_int)
coef(msm_cde)
# (Intercept)           A0_ace        M_smoking A0_ace:M_smoking
#  0.17891689       0.06798282       0.06729724      -0.00495314

## 5. Estimate CDE for m=0 and for m=1 using the MSM's coefficients
CDE_mis0 <- coef(msm_cde)["A0_ace"]
# 0.06798282
CDE_mis1 <- coef(msm_cde)["A0_ace"] + coef(msm_cde)["A0_ace:M_smoking"]
# 0.06302968


rm(list=ls())
df1_int <- read.csv(file = "data/df1_int.csv")
### MSM of CDE, estimated by G-computation -------------------------------------
## 1. Estimate Qbar(A,M,L0,L1)
Q.cde.death <- glm(Y_death ~ A0_ace + M_smoking + A0_ace:M_smoking + L0_male + L0_parent_low_educ_lv + L1,
                   family = "gaussian", data = df1_int)
# The final result would be sligthly different if we applied a binomial family
# The Gaussian family corresponds to the true generating model in this example.

## 2. Predict an outcome for each subject, in each counterfactual scenario
# prepare data sets used to predict the outcome under the counterfactual
# 4 counterfactual scenarios setting (A=0,M=0), (A=1,M=0), (A=0,M=1) and (A=1,M=1)
data.A0M0 <- data.A1M0 <- data.A0M1 <- data.A1M1 <- df1_int
data.A0M0$A0_ace <- 0
data.A0M0$M_smoking <- 0

data.A1M0$A0_ace <- 1
data.A1M0$M_smoking <- 0

data.A0M1$A0_ace <- 0
data.A0M1$M_smoking <- 1

data.A1M1$A0_ace <- 1
data.A1M1$M_smoking <- 1

# predict values under the same name in the corresponding counterfactual dataset
data.A0M0$Yam.death.pred <- predict(Q.cde.death, newdata = data.A0M0, type = "response")
data.A1M0$Yam.death.pred <- predict(Q.cde.death, newdata = data.A1M0, type = "response")
data.A0M1$Yam.death.pred <- predict(Q.cde.death, newdata = data.A0M1, type = "response")
data.A1M1$Yam.death.pred <- predict(Q.cde.death, newdata = data.A1M1, type = "response")

## 3. Append both counterfactual datasets in a single dataset
# number of row is twice the initial value (we have 2 counterfactual scenarios)
data.4scenarios <- rbind(data.A0M0, data.A1M0,data.A0M1,data.A1M1)

## 4. fit the MSM: E(Y_a|sex)
MSM.CDE.gcomp <- glm(Yam.death.pred ~ A0_ace +  M_smoking + A0_ace:M_smoking,
                     family = "gaussian", # gaussian family for risk differences
                     data = data.4scenarios)
coef(MSM.CDE.gcomp)
# (Intercept)           A0_ace        M_smoking A0_ace:M_smoking
#  0.17968603       0.06000138       0.06757214       0.01918153

## 5. Estimate the CDEm
# CDE(M=0) = E(Y_{A=1,M=0}) - E(Y_{A=0,M=0})
CDE_mis0_gcomp <- coef(MSM.CDE.gcomp)["A0_ace"]
# 0.06000138

# CDE(M=1) = E(Y_{A=1,M=1}) - E(Y_{A=0,M=1})
CDE_mis1_gcomp <- (coef(MSM.CDE.gcomp)["A0_ace"] +
                     coef(MSM.CDE.gcomp)["A0_ace:M_smoking"])
# 0.07918291

# Note: Applying a binomial family for the first Qbar model would result in two
# sligthly different values of the CDE(M=m)
# => 0.05934409 in setting M=0
# => 0.07537874 in setting M=1


rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
### MSM of CDE, estimated by G-computation (by ICE) ----------------------------
## 1a) Regress the outcome on L0, A, L1 and M (and the A*M interaction if appropriate)
Y.death.model <- glm(Y_death ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 +
                       M_smoking + A0_ace:M_smoking,
                     family = "binomial", data = df2_int)

## 1b) Generate predicted values by evaluating the regression setting the mediator
##    value to M=0 or to M=1
data.A0M0 <- data.A1M0 <- data.A0M1 <- data.A1M1 <- df2_int
data.A0M0$A0_ace <- 0
data.A0M0$M_smoking <- 0

data.A1M0$A0_ace <- 1
data.A1M0$M_smoking <- 0

data.A0M1$A0_ace <- 0
data.A0M1$M_smoking <- 1

data.A1M1$A0_ace <- 1
data.A1M1$M_smoking <- 1

Q.Y.death.A0M0 <- predict(Y.death.model, newdata = data.A0M0, type = "response")
Q.Y.death.A1M0 <- predict(Y.death.model, newdata = data.A1M0, type = "response")
Q.Y.death.A0M1 <- predict(Y.death.model, newdata = data.A0M1, type = "response")
Q.Y.death.A1M1 <- predict(Y.death.model, newdata = data.A1M1, type = "response")

## 2a) Regress the predicted values conditional on the observed exposure A
##    and baseline confounders L(0)
L1.death.A0M0.model <- glm(Q.Y.death.A0M0 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                           family = "quasibinomial", data = df2_int)
L1.death.A1M0.model <- glm(Q.Y.death.A1M0 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                           family = "quasibinomial", data = df2_int)
L1.death.A0M1.model <- glm(Q.Y.death.A0M1 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                           family = "quasibinomial", data = df2_int)
L1.death.A1M1.model <- glm(Q.Y.death.A1M1 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
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
MSM.CDE.gcomp <- glm(Yam.death.pred ~ A0_ace +  M_smoking + A0_ace:M_smoking,
                     family = "gaussian", # gaussian family for risk differences
                     data = data.4scenarios)
coef(MSM.CDE.gcomp)
# (Intercept)           A0_ace        M_smoking A0_ace:M_smoking
#  0.17974947       0.06342833       0.07366466       0.02469485

## 5. Estimate the CDE(M=m)
# CDE(M=0) = E(Y_{A=1,M=0}) - E(Y_{A=0,M=0})
CDE_mis0_gcomp_ice <- coef(MSM.CDE.gcomp)["A0_ace"]
# 0.06342833

# CDE(M=1) = E(Y_{A=1,M=1}) - E(Y_{A=0,M=1})
CDE_mis1_gcomp_ice <- (coef(MSM.CDE.gcomp)["A0_ace"] +
                         coef(MSM.CDE.gcomp)["A0_ace:M_smoking"])
# 0.08812318


##### Alternative --------------------------------------------------------------
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
### MSM of CDE, estimated by G-computation (by ICE) ----------------------------
## 1a) Regress the outcome on L0, A, L1 and M (and the A*M interaction if appropriate)
Y.death.model <- glm(Y_death ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 +
                       M_smoking + A0_ace:M_smoking,
                     family = "binomial", data = df2_int)

## 1b) Generate predicted values by evaluating the regression setting the mediator
##    value to M=0 or to M=1
data.A0M0 <- data.A1M0 <- data.A0M1 <- data.A1M1 <- df2_int
# data.A0M0$A0_ace <- 0
data.A0M0$M_smoking <- 0

# data.A1M0$A0_ace <- 1
data.A1M0$M_smoking <- 0

# data.A0M1$A0_ace <- 0
data.A0M1$M_smoking <- 1

# data.A1M1$A0_ace <- 1
data.A1M1$M_smoking <- 1

data.A0M0$Q.Y.death <- predict(Y.death.model, newdata = data.A0M0, type = "response")
data.A1M0$Q.Y.death <- predict(Y.death.model, newdata = data.A1M0, type = "response")
data.A0M1$Q.Y.death <- predict(Y.death.model, newdata = data.A0M1, type = "response")
data.A1M1$Q.Y.death <- predict(Y.death.model, newdata = data.A1M1, type = "response")

## 2a) Regress the predicted values conditional on the exposure A
##    and baseline confounders L(0)
L1.death.A0M0.model <- glm(Q.Y.death ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                           family = "quasibinomial", data = data.A0M0)
L1.death.A1M0.model <- glm(Q.Y.death ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                           family = "quasibinomial", data = data.A1M0)
L1.death.A0M1.model <- glm(Q.Y.death ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                           family = "quasibinomial", data = data.A0M1)
L1.death.A1M1.model <- glm(Q.Y.death ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                           family = "quasibinomial", data = data.A1M1)


## 2b) generate predicted values by evaluating the regression at exposure
##    of interest: {A=0,M=0}, {A=1,M=0}, {A=0,M=1}, {A=1,M=1}
data.A0M0$A0_ace <- 0
data.A1M0$A0_ace <- 1
data.A0M1$A0_ace <- 0
data.A1M1$A0_ace <- 1
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
MSM.CDE.gcomp <- glm(Yam.death.pred ~ A0_ace +  M_smoking + A0_ace:M_smoking,
                     family = "gaussian", # gaussian family for risk differences
                     data = data.4scenarios)
coef(MSM.CDE.gcomp)
# (Intercept)           A0_ace        M_smoking A0_ace:M_smoking
#  0.17974978       0.06341297       0.07366464       0.02469211

## 5. Estimate the CDE(M=m)
# CDE(M=0) = E(Y_{A=1,M=0}) - E(Y_{A=0,M=0})
CDE_mis0_gcomp_ice <- coef(MSM.CDE.gcomp)["A0_ace"]
# 0.06341297

# CDE(M=1) = E(Y_{A=1,M=1}) - E(Y_{A=0,M=1})
CDE_mis1_gcomp_ice <- (coef(MSM.CDE.gcomp)["A0_ace"] +
                         coef(MSM.CDE.gcomp)["A0_ace:M_smoking"])
# 0.08810508


##### ltmle  --------------------------------------------------------------
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
library(ltmle)
Qform <- c(L1="Q.kplus1 ~ L0_male + L0_parent_low_educ_lv + A0_ace",
           Y_death="Q.kplus1 ~ L0_male + L0_parent_low_educ_lv + L1 +
                    A0_ace * M_smoking")
gform <- c("A0_ace ~ L0_male + L0_parent_low_educ_lv",
           "M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1")
data_binary <- subset(df2_int, select = c(L0_male, L0_parent_low_educ_lv,
                                          A0_ace, L1,
                                          M_smoking, Y_death))
CDE_ltmle_M0_death <- ltmle(data = data_binary,
                            Anodes = c("A0_ace", "M_smoking"),
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
                            Anodes = c("A0_ace", "M_smoking"),
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
### MSM of NDE & NIE, estimated by IPTW ----------------------------------------
## 1. Stabilized weight for the MSM1
# 1a. sw_Ai = g(A=a_i | L(0)) / g(A=a_i | L(0)) = 1
sw_Ai <- rep(1, nrow(df1_int))

# 1b. sw_Mi = g(M=m_i | A,L(0)) / g(M=m_i | A,L(0),L(1))
g.M.AL0 <- glm(M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv,
               family = "binomial", data = df1_int)
g.Mis1.AL0 <- predict(g.M.AL0, type = "response")
sw_M.num <- rep(NA, nrow(df1_int))
sw_M.num[df1_int$M_smoking==1] <- g.Mis1.AL0[df1_int$M_smoking==1]
sw_M.num[df1_int$M_smoking==0] <- (1 - g.Mis1.AL0[df1_int$M_smoking==0])

g.M.AL0L1 <- glm(M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv + L1,
                 family = "binomial", data = df1_int)
g.Mis1.AL0L1 <- predict(g.M.AL0L1, type = "response")
sw_M.denom <- rep(NA, nrow(df1_int))
sw_M.denom[df1_int$M_smoking==1] <- g.Mis1.AL0L1[df1_int$M_smoking==1]
sw_M.denom[df1_int$M_smoking==0] <- (1 - g.Mis1.AL0L1[df1_int$M_smoking==0])

sw_msm1 <- sw_Ai * sw_M.num / sw_M.denom

## 2. Estimate coefficients of the MSM1
MSM1 <- glm(Y_death ~ A0_ace + M_smoking + A0_ace:M_smoking +
              L0_male + L0_parent_low_educ_lv,
            weights = sw_msm1,
            family = "gaussian",
            data = df1_int)
coef(MSM1)
# (Intercept)                A0_ace             M_smoking
#  0.12033221            0.06381257            0.06691712
#     L0_male L0_parent_low_educ_lv      A0_ace:M_smoking
#  0.04671886            0.05521263            0.01652446

## 3. Stabilized weight for the MSM2
# 3a. sw_A = g(A=a_i) / g(A=a_i | L(0))
# numerator
g.A <- glm(A0_ace ~ 1, family = "binomial", data = df1_int)
g.Ais1 <- predict(g.A, type = "response")
sw_msm2.num <- rep(NA, nrow(df1_int))
sw_msm2.num[df1_int$A0_ace==1] <- g.Ais1[df1_int$A0_ace==1]
sw_msm2.num[df1_int$A0_ace==0] <- (1 - g.Ais1[df1_int$A0_ace==0])

# denominator
g.A.L0 <- glm(A0_ace ~ L0_male + L0_parent_low_educ_lv,
              family = "binomial", data = df1_int)
g.Ais1.L0 <- predict(g.A.L0, type = "response")
sw_msm2.denom <- rep(NA, nrow(df1_int))
sw_msm2.denom[df1_int$A0_ace==1] <- g.Ais1.L0[df1_int$A0_ace==1]
sw_msm2.denom[df1_int$A0_ace==0] <- (1 - g.Ais1.L0[df1_int$A0_ace==0])

# stabilized weight
sw_msm2 <- sw_msm2.num / sw_msm2.denom

## 3. Estimate coefficients of the MSM2
MSM2 <- glm(M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv,
            weights = sw_msm2,
            family = "binomial",
            data = df1_int)
coef(MSM2)
# (Intercept)                A0_ace               L0_male L0_parent_low_educ_lv
#  -1.2723106             0.5883720             0.2566129             0.3270087

## 4. Estimate PNDE conditional on L(0), and the marginal value of PNDE
# a = 1 and a* = 0
# PNDE|L(0) = (a - a*)[alpha_A + alpha_AM.g^-1(a^*,l(0))]
g.minus1.A0 <- plogis(coef(MSM2)["(Intercept)"] + coef(MSM2)["A0_ace"] * 0 +
                      coef(MSM2)["L0_male"] * df1_int$L0_male +
                      coef(MSM2)["L0_parent_low_educ_lv"] * df1_int$L0_parent_low_educ_lv)

# PNDE conditional on L(0)
PNDE_L0 <- (1 - 0) * (coef(MSM1)["A0_ace"] +
                        coef(MSM1)["A0_ace:M_smoking"] * g.minus1.A0)
# marginal PNDE
PNDE <- mean(PNDE_L0)
# [1] 0.06850657

## 4. Estimate TNIE conditional on L(0), and the marginal value of TNIE
# TNIE|L(0) = [g^-1(a,l(0)) - g^-1(a^*,l(0))] * (alpha_M + alpha_AM * a)
g.minus1.A1 <- plogis(coef(MSM2)["(Intercept)"] + coef(MSM2)["A0_ace"] * 1 +
                        coef(MSM2)["L0_male"] * df1_int$L0_male +
                        coef(MSM2)["L0_parent_low_educ_lv"] * df1_int$L0_parent_low_educ_lv)

# TNIE conditional on L(0)
TNIE_L0 <- (g.minus1.A1 - g.minus1.A0) * (coef(MSM1)["M_smoking"] +
                                            coef(MSM1)["A0_ace:M_smoking"] * 1)
# marginal PNDE
TNIE <- mean(TNIE_L0)
# [1] 0.01096799


# ---------------------------------------------------------------------------- #
# IV) CMAverse -----------
# ---------------------------------------------------------------------------- #
# The DAG for this scientific setting is:
# devtools::install_github("BS1125/CMAverse")
library(CMAverse)
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
cmdag(outcome = "Y_death", exposure = "A0_ace", mediator = "M_smoking",
      basec = c("L0_male", "L0_parent_low_educ_lv"), postc = "L1", node = TRUE, text_col = "white")
# In this setting, we can use the marginal structural model and the $g$-formula approach. The results are shown below.

## The Marginal Structural Model
res_msm <- cmest(data = df2_int,
                 model = "msm",
                 outcome = "Y_death",
                 exposure = "A0_ace",
                 mediator = "M_smoking",
                 basec = c("L0_male", "L0_parent_low_educ_lv"),
                 postc = "L1",
                 EMint = TRUE, # E*M interaction
                 ereg = "logistic", # exposure regression model g(A=1|L(0))
                 # yreg = "logistic", not necessary with MSM
                 mreg = list("logistic"), # mediation regresion model g(M=1|L1,A,L0)
                 wmnomreg = list("logistic"), #g(M=1|A)
                 wmdenomreg = list("logistic"), # g(M=1|L1,A,L(0))
                 astar = 0, #E(Y_{A=0,M=1})
                 a = 1,  #E(Y_{A=1,M=1})
                 mval = list(1), # mediator value at which the variable is controlled
                 estimation = "imputation",
                 inference = "bootstrap",
                 nboot = 2)
summary(res_msm)
# Causal Mediation Analysis
#
# # Outcome regression:
#
# Call:
#   glm(formula = Y_death ~ A0_ace + M_smoking + A0_ace * M_smoking,
#       family = binomial(), data = getCall(x$reg.output$yreg)$data,
#       weights = getCall(x$reg.output$yreg)$weights)
#
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -1.6420  -0.6975  -0.5999  -0.5638   2.6288
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)      -1.52095    0.03362 -45.238  < 2e-16 ***
#   A0_ace            0.36526    0.10056   3.632 0.000281 ***
#   M_smoking         0.43148    0.05470   7.888 3.08e-15 ***
#   A0_ace:M_smoking  0.06991    0.14299   0.489 0.624927
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
# Null deviance: 10335  on 9999  degrees of freedom
# Residual deviance: 10219  on 9996  degrees of freedom
# AIC: 10340
#
# Number of Fisher Scoring iterations: 4
#
#
# # Mediator regressions:
#
# Call:
#   glm(formula = M_smoking ~ A0_ace, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
#
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -1.4730  -0.8841  -0.8739   1.4592   1.7251
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.73432    0.02268 -32.384  < 2e-16 ***
#   A0_ace       0.51531    0.06422   8.024 1.02e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
# Null deviance: 12789  on 9999  degrees of freedom
# Residual deviance: 12726  on 9998  degrees of freedom
# AIC: 12847
#
# Number of Fisher Scoring iterations: 4
#
#
# # Mediator regressions for weighting (denominator):
#
# Call:
#   glm(formula = M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv +
#         L1, family = binomial(), data = getCall(x$reg.output$wmdenomreg[[1L]])$data,
#       weights = getCall(x$reg.output$wmdenomreg[[1L]])$weights)
#
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -1.3340  -0.8581  -0.7529   1.2616   1.7835
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -1.36249    0.04783 -28.488  < 2e-16 ***
#   A0_ace                 0.30994    0.06668   4.648 3.35e-06 ***
#   L0_male                0.24661    0.04369   5.644 1.66e-08 ***
#   L0_parent_low_educ_lv  0.30628    0.04650   6.587 4.50e-11 ***
#   L1                     0.86045    0.04493  19.152  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
# Null deviance: 12783  on 9999  degrees of freedom
# Residual deviance: 12260  on 9995  degrees of freedom
# AIC: 12270
#
# Number of Fisher Scoring iterations: 4
#
#
# # Mediator regressions for weighting (nominator):
#
# Call:
#   glm(formula = M_smoking ~ A0_ace, family = binomial(), data = getCall(x$reg.output$wmnomreg[[1L]])$data,
#       weights = getCall(x$reg.output$wmnomreg[[1L]])$weights)
#
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -1.0982  -0.8825  -0.8825   1.5043   1.5043
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.74205    0.02271 -32.680   <2e-16 ***
#   A0_ace       0.55288    0.06408   8.628   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
# Null deviance: 12783  on 9999  degrees of freedom
# Residual deviance: 12710  on 9998  degrees of freedom
# AIC: 12714
#
# Number of Fisher Scoring iterations: 4
#
#
# # Exposure regression for weighting:
#
# Call:
#   glm(formula = A0_ace ~ L0_male + L0_parent_low_educ_lv, family = binomial(),
#       data = getCall(x$reg.output$ereg)$data, weights = getCall(x$reg.output$ereg)$weights)
#
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -0.5830  -0.5830  -0.4825  -0.3550   2.3645
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -2.73244    0.07425 -36.799  < 2e-16 ***
#   L0_male                0.40580    0.06447   6.294 3.09e-10 ***
#   L0_parent_low_educ_lv  0.64060    0.07350   8.716  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
# Null deviance: 7030.1  on 9999  degrees of freedom
# Residual deviance: 6909.0  on 9997  degrees of freedom
# AIC: 6915
#
# Number of Fisher Scoring iterations: 5
#
#
# # Effect decomposition on the odds ratio scale via the marginal structural model
#
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
#                   Estimate Std.error  95% CIL 95% CIU  P.val
#   Rcde            1.545217  0.055313 1.629922   1.704 <2e-16 ***
#   rRpnde          1.473363  0.103755 1.627216   1.767 <2e-16 ***
#   rRtnde          1.486793  0.071520 1.640039   1.736 <2e-16 ***
#   rRpnie          1.059492  0.012056 1.047202   1.063 <2e-16 ***
#   rRtnie          1.069150  0.007743 1.045052   1.055 <2e-16 ***
#   Rte             1.575246  0.095828 1.717453   1.846 <2e-16 ***
#   ERcde           0.446294  0.020216 0.520256   0.547 <2e-16 ***
#   rERintref       0.027069  0.123971 0.079953   0.247 <2e-16 ***
#   rERintmed       0.042391  0.019983 0.016156   0.043 <2e-16 ***
#   rERpnie         0.059492  0.012056 0.047206   0.063 <2e-16 ***
#   ERcde(prop)     0.775832  0.110347 0.615334   0.764 <2e-16 ***
#   rERintref(prop) 0.047056  0.133964 0.110548   0.291 <2e-16 ***
#   rERintmed(prop) 0.073692  0.030419 0.019257   0.060 <2e-16 ***
#   rERpnie(prop)   0.103420  0.006802 0.065741   0.075 <2e-16 ***
#   rpm             0.177112  0.023618 0.094136   0.126 <2e-16 ***
#   rint            0.120748  0.103545 0.170674   0.310 <2e-16 ***
#   rpe             0.224168  0.110347 0.236415   0.385 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Rcde: controlled direct effect odds ratio;
#  rRpnde: randomized analogue of pure natural direct effect odds ratio;
#  rRtnde: randomized analogue of total natural direct effect odds ratio;
#  rRpnie: randomized analogue of pure natural indirect effect odds ratio;
#  rRtnie: randomized analogue of total natural indirect effect odds ratio;
#  Rte: total effect odds ratio;
#  ERcde: excess relative risk due to controlled direct effect;
#  rERintref: randomized analogue of excess relative risk due to reference interaction;
#  rERintmed: randomized analogue of excess relative risk due to mediated interaction;
#  rERpnie: randomized analogue of excess relative risk due to pure natural indirect effect;
#  ERcde(prop): proportion ERcde;
#  rERintref(prop): proportion rERintref;
#  rERintmed(prop): proportion rERintmed;
#  rERpnie(prop): proportion rERpnie;
#  rpm: randomized analogue of overall proportion mediated;
#  rint: randomized analogue of overall proportion attributable to interaction;
#  rpe: randomized analogue of overall proportion eliminated)
#
# Relevant variable values:
#   $a
# [1] 1
#
# $astar
# [1] 0
#
# $yval
# [1] "1"
#
# $mval
# $mval[[1]]
# [1] 1

# can we get risk differences ?
res_msm_RD <- cmest(data = df2_int,
                    model = "msm",
                    outcome = "Y_death",
                    exposure = "A0_ace",
                    mediator = "M_smoking",
                    basec = c("L0_male", "L0_parent_low_educ_lv"),
                    postc = "L1",
                    EMint = TRUE, # E*M interaction
                    ereg = "logistic",
                    yreg = "linear", # estimate the binary outcome by linear regression
                    mreg = list("logistic"),
                    wmnomreg = list("logistic"),
                    wmdenomreg = list("logistic"),
                    astar = 0, #E(Y_{A=0,M=1})
                    a = 1,  #E(Y_{A=1,M=1})
                    mval = list(1), # mediator value at which the variable is controlled (for the CDE)
                    estimation = "imputation",
                    inference = "bootstrap",
                    nboot = 2)

summary(res_msm_RD)
# Estimate  Std.error    95% CIL 95% CIU  P.val
#   cde            0.0903019  0.0311309  0.0303279   0.072 <2e-16 *** # 0.09023883 if we apply a gaussian regression for the Q.L2 quantity in the ICE method for CDE
#   rpnde          0.0695928  0.0292988  0.0495508   0.089 <2e-16 ***
#   rtnde          0.0736457  0.0299728  0.0458806   0.086 <2e-16 ***
#   rpnie          0.0097246  0.0009589  0.0084005   0.010 <2e-16 ***
#   rtnie          0.0137775  0.0002849  0.0056358   0.006 <2e-16 ***
#   te             0.0833703  0.0290139  0.0555693   0.095 <2e-16 ***
#   rintref       -0.0207091  0.0018322  0.0167614   0.019 <2e-16 ***
#   rintmed        0.0040529  0.0006740 -0.0036703  -0.003 <2e-16 ***
#   cde(prop)      1.0831421  0.1630518  0.5416174   0.761 <2e-16 ***
#   rintref(prop) -0.2483987  0.1265182  0.1791686   0.349 <2e-16 ***
#   rintmed(prop)  0.0486133  0.0276120 -0.0667509  -0.030 <2e-16 ***
#   rpnie(prop)    0.1166433  0.0641455  0.0898076   0.176 <2e-16 ***
#   rpm            0.1652566  0.0365335  0.0601534   0.109 <2e-16 ***
#   rint          -0.1997854  0.0989062  0.1495145   0.282 <2e-16 ***
#   rpe           -0.0831421  0.1630518  0.2393220   0.458 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# cde: controlled direct effect;
# rpnde: randomized analogue of pure natural direct effect;
# rtnde: randomized analogue of total natural direct effect;
# rpnie: randomized analogue of pure natural indirect effect;
# rtnie: randomized analogue of total natural indirect effect;
# te: total effect;
# rintref: randomized analogue of reference interaction;
# rintmed: randomized analogue of mediated interaction;
# cde(prop): proportion cde;
# rintref(prop): proportion rintref;
# rintmed(prop): proportion rintmed;
# rpnie(prop): proportion rpnie;
# rpm: randomized analogue of overall proportion mediated;
# rint: randomized analogue of overall proportion attributable to interaction;
# rpe: randomized analogue of overall proportion eliminated

res_msm_RD_M0 <- cmest(data = df2_int,
                       model = "msm",
                       outcome = "Y_death",
                       exposure = "A0_ace",
                       mediator = "M_smoking",
                       basec = c("L0_male", "L0_parent_low_educ_lv"),
                       postc = "L1",
                       EMint = TRUE, # E*M interaction
                       ereg = "logistic",
                       yreg = "linear", # estimate the binary outcome by linear regression
                       mreg = list("logistic"),
                       # wmnomreg = list("logistic"),
                       # wmdenomreg = list("logistic"),
                       astar = 0, #E(Y_{A=0,M=1})
                       a = 1,  #E(Y_{A=1,M=1})
                       mval = list(0), # set do(M=0) for the CDE
                       estimation = "imputation",
                       inference = "bootstrap",
                       nboot = 2)

summary(res_msm_RD_M0)
#                  Estimate Std.error   95% CIL 95% CIU  P.val
#   cde           0.0601292 0.0130764 0.0227365   0.040 <2e-16 ***
#   rpnde         0.0698219 0.0157661 0.0412164   0.062 <2e-16 ***
#   rtnde         0.0735433 0.0164906 0.0481480   0.070 <2e-16 ***
#   rpnie         0.0089292 0.0007660 0.0076369   0.009 <2e-16 ***
#   rtnie         0.0126505 0.0000416 0.0155418   0.016 <2e-16 ***
#   te            0.0824725 0.0157245 0.0568141   0.078 <2e-16 ***
#   rintref       0.0096927 0.0026897 0.0184798   0.022 <2e-16 ***
#   rintmed       0.0037214 0.0007244 0.0069316   0.008 <2e-16 ***
#   cde(prop)     0.7290822 0.0872741 0.3990139   0.516 <2e-16 ***
#   rintref(prop) 0.1175267 0.0311988 0.2837739   0.326 <2e-16 ***
#   rintmed(prop) 0.0451228 0.0153616 0.1015740   0.122 <2e-16 ***
#   rpnie(prop)   0.1082682 0.0407137 0.0983853   0.153 <2e-16 ***
#   rpm           0.1533911 0.0560753 0.1999592   0.275 <2e-16 ***
#   rint          0.1626495 0.0465604 0.3853478   0.448 <2e-16 ***
#   rpe           0.2709178 0.0872741 0.4837331   0.601 <2e-16 ***
