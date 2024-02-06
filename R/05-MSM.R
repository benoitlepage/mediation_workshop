### test programs for estimations by MSMs

rm(list=ls())

df1 <- read.csv(file = "data/df1.csv")
df1_int <- read.csv(file = "data/df1_int.csv")

################################################################################
######################### Estimation of the Average Total Effect (ATE)
################################################################################

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


# 3. Define individual weights:
# We can use simple weights w = 1 / g(A=a_i | L(0))
w <- 1 / gAi.L

# Or alternatively, we can use stabilized weights : sw = g(A=a_i) / g(A=a_i | L(0))
sw <- gAi.sex / gAi.L

par(mfcol = c(1,2))
boxplot(w ~ df1_int$A0_ace)
boxplot(sw ~ df1_int$A0_ace)


# 4. Estimate coefficients of the MSM using a weighted regression E(Y | A, sex)
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

