### test programs for estimations based on TMLE

rm(list=ls())

df2_int <- read.csv(file = "data/df2_int.csv")

# ---------------------------------------------------------------------------- #
# Estimation of the Average Total Effect (ATE) ----
# ---------------------------------------------------------------------------- #
library(ltmle)

Qform <- c(Y_death="Q.kplus1 ~ L0_male + L0_parent_low_educ_lv + A0_ace")
gform <- c("A0_ace ~ L0_male + L0_parent_low_educ_lv")

Psi_Ais1 <- ltmle(data = subset(df2_int, select = c(L0_male,
                                                    L0_parent_low_educ_lv,
                                                    A0_ace,
                                                    Y_death)),
                  Anodes = "A0_ace",
                  Ynodes = "Y_death",
                  Qform = Qform,
                  gform = gform,
                  abar = 1,
                  SL.library = "glm",
                  variance.method = "ic")
summary(Psi_Ais1, "tmle")
# Parameter Estimate:  0.28714
#  Estimated Std Err:  0.013838
#            p-value:  <2e-16
#  95% Conf Interval: (0.26002, 0.31426)

head(Psi_Ais1$cum.g, 10)
# [,1]
# [1,] 0.10989220
# [2,] 0.15629749
# [3,] 0.15629749
# [4,] 0.08894074
# [5,] 0.15629749
# [6,] 0.15629749
# [7,] 0.10989220
# [8,] 0.10989220
# [9,] 0.15629749
# [10,] 0.15629749

table(Psi_Ais1$cum.g.used)
# FALSE  TRUE
# 8876  1124

Psi_Ais1$fit$Qstar
# Coefficients:
#   S1
# -0.001667
# Degrees of Freedom: 1124 Total (i.e. Null);  1123 Residual

head(Psi_Ais1$IC$tmle)
# [1]  0.0001003559  4.2532581791  0.0569935644 -0.0280052148  0.0569935644  0.0569935644


# ---------------------------------------------------------------------------- #
## manual calculation: ----
# ---------------------------------------------------------------------------- #
## 1) Estimate Qbar and predict Qbar when A0_ace is set to 1
Q.fit <- glm(Y_death ~ A0_ace + L0_male + L0_parent_low_educ_lv,
             family = "binomial", data = df2_int)
data.A1 <- df2_int
data.A1$A0_ace <- 1
# Qbar_Ais1 <- predict(Q.fit, newdata = data.A1, type = "response")
logitQ <- predict(Q.fit, newdata = data.A1, type = "link")

## 2) Estimate the treatment mechanism
g.L <- glm(A0_ace ~ L0_male + L0_parent_low_educ_lv,
           family = "binomial", data = df2_int)

# predict the probabilities P(A0_ace=1|L(0))
g1.L <- predict(g.L, type="response")
head(g1.L)
#          1          2          3          4          5          6
# 0.10989220 0.15629749 0.15629749 0.08894074 0.15629749 0.15629749

# pred.g1.L <- predict(g.L, type="response")
# pred.g0.L <- 1 - pred.g1.L

# the predicted probability of the observed treatment A_i=a is :
# gA.L <- rep(NA, nrow(df2_int))
# gA.L[df2_int$A0_ace==1] <- pred.g1.L[df2_int$A0_ace==1]
# gA.L[df2_int$A0_ace==0] <- pred.g0.L[df2_int$A0_ace==0]

# its useful to check the distribution of gA.L, as values close to 0 or 1
# (< 0.01 or > 0.999)
# are indicators of near positivity violation
summary(g1.L)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.06109 0.84370 0.89011 0.80279 0.91106 0.93891
# there is no positivity issues in this example.
# head(gA.L, 10)
# # [1] 0.8901078 0.1562975 0.8437025 0.9110593 0.8437025 0.8437025 0.8901078 0.8901078 0.1562975 0.8437025
#
#
# plot(gA.L, Psi_Ais1$cum.g)
# plot(gA.L[Psi_Ais1$cum.g.used == TRUE], Psi_Ais1$cum.g[Psi_Ais1$cum.g.used == TRUE])
# en fait, on a pas besoin de calculer la probabilité d'observer sa propre observation, avec l'indicatrice, on retombe sur nos pieds !


## 3) Determine a parametric family of fluctuations of Qbar,
#     defined by parametric regression including a clever covarate chosen so the
#     loss function L(Q(epsilon)) holds with epsilon playing the role of the coefficient
#     in front of the clever covariate.

# The clever covariate H(A,L(0)) depends on g:
H <- (df2_int$A0_ace == 1) / g1.L


## 3) Update the initial fit Qbar from step 1.
# This is achieved by holding Qbar fixed (as intercept) while estimating the
# coefficient epsilon for H

# for example:
update.fit <- glm(df2_int$Y_death ~ -1 + offset(logitQ) + H,
                  family = "quasibinomial")
# Coefficients:
#   H
# -0.0001756
# Degrees of Freedom: 10000 Total (i.e. Null);  9999 Residual
Qstar <- predict(update.fit, data = data.frame(logitQ, H), type = "response")

# in the ltmle package the fluctuation parametric model is the following:
S1 <- rep(1, nrow(df2_int))
update.fit.ltmle <- glm(df2_int$Y_death ~ -1 + S1 + offset(logitQ),
                        family = "quasibinomial",
                        weights = scale(H, center = FALSE))
# Coefficients:
#   S1
# -0.001667 #as with the ltmle package
# Degrees of Freedom: 1124 Total (i.e. Null);  1123 Residual
Qstar.tmle <- predict(update.fit.ltmle, data = data.frame(logitQ, H), type = "response")


## 4) Obtain the substition estimator of Psi_Ais1
Psi <- mean(Qstar)
# [1] 0.2874431

Psi.tmle <- mean(Qstar.tmle)
# [1] 0.2871408 # as with the ltmle package

## 5) Calculate standard errors based on the influence curve of the TMLE
IC <- H * (df2_int$Y_death - Qstar.tmle) + Qstar.tmle - Psi.tmle
head(IC)
# 0.0001003559  4.2532581791  0.0569935644 -0.0280052148  0.0569935644  0.056993564
# standard error can be estimated by :
sqrt(var(IC)/nrow(df2_int))
# [1] 0.01383821

plot(IC,Psi_Ais1$IC$tmle)
abline(a = 0, b = 1)




#### ATE With SuperLearner
library(SuperLearner)
library(xgboost)
# Below, we use the same ltmle() function than previously,
# and specify our family of algorithms to be used with the SuperLearner

## we can change the default argument of the SL.xgboost algorithm and the
## SL.step.interaction algorithm

# We can check how arguments are used in the pre-specified algorithms
SL.step.interaction
# function (Y, X, newX, family, direction = "both", trace = 0,
#     k = 2, ...)
# {
#     fit.glm <- glm(Y ~ ., data = X, family = family)
#     fit.step <- step(fit.glm, scope = Y ~ .^2, direction = direction,
#         trace = trace, k = k)
#     pred <- predict(fit.step, newdata = newX, type = "response")
#     fit <- list(object = fit.step)
#     out <- list(pred = pred, fit = fit)
#     class(out$fit) <- c("SL.step")
#     return(out)
# }
# <bytecode: 0x000001b965ed0dc0>
# <environment: namespace:SuperLearner>

# The pre-specified can be easily modified to obtain a simple backward selection
SL.interaction.back = function(...) {
  SL.step.interaction(..., direction = "backward")
}

# The same principle can be applied with xgboost
SL.xgboost
SL.xgboost.custom = function(...) {
  SL.xgboost(..., ntrees = 50)
}



# the data-adaptive algorithms can be specified separately for the Q and g functions
SL.library <- list(Q=c("SL.mean","SL.glm","SL.interaction.back", "SL.xgboost.custom"),
                   g=c("SL.mean","SL.glm","SL.interaction.back", "SL.xgboost.custom"))

set.seed(42)
Psi_ATE_tmle <- ltmle(data = data_ltmle,
                      Anodes = "A0_ace",
                      Ynodes = "Y_death",
                      Qform = Qform,
                      gform = gform,
                      gbounds = c(0.01, 1),
                      abar = list(1,0), # vector of the counterfactual treatment
                      SL.library = SL.library,
                      variance.method = "ic")
summary(Psi_ATE_tmle)
# The function give the ATE on the difference scale (as well, as RR and OR)
# Additive Treatment Effect:
# Parameter Estimate:  0.081832
#  Estimated Std Err:  0.014291
#            p-value:  1.0275e-08
#  95% Conf Interval: (0.053822, 0.10984)

## We can see how the SuperLearner used the algorithms for the g function
Psi_ATE_tmle$fit$g
# [[1]]$A0_ace
#                               Risk        Coef
# SL.mean_All             0.09976892 0.003545569 # risk is higher for the bad model
# SL.glm_All              0.09865424 0.416238369
# SL.interaction.back_All 0.09865424 0.000000000
# SL.xgboost.custom_All   0.09865550 0.580216062

# for the g function, the SuperLearner predicts the treatment mechanism
# base on a mix between the glm and the customized xgboost algorithm.

## We can see how the SuperLearner used the algorithms for the g function
Psi_ATE_tmle$fit$Q
#                              Risk       Coef
# SL.mean_All             0.1684737 0.02003166 # risk is higher for the bad model
# SL.glm_All              0.1662241 0.00000000
# SL.interaction.back_All 0.1662241 0.55956284
# SL.xgboost.custom_All   0.1662422 0.42040550

# The SuperLearner predicts both the treatment mechanism g and the Q function
# from a mix between the backward interaction glm (or the main term glm) and the
# customized xgboost algorithm.
# However, the choice between the SL.glm and the SL.interaction.back
# procedure was arbitrary: as we can see the Risk is exactly the same for both
# algorithms. The final model from the step-by-step procedure was much probably
# a main term glm.

# ---------------------------------------------------------------------------- #
## CDE using ltmle ----
# ---------------------------------------------------------------------------- #
### CDE
## for binary outcomes
library(ltmle)
# We define the formulas for the estimation of the 2 barQ functions
# Note that it is possible to specify the A*M interaction, if we really want to
# take it into account.
# Another option is to indicate prediction algorithms well adapted to the estimation
# of interaction phenomena into the SuperLearner arguments.
Qform <- c(L1="Q.kplus1 ~ L0_male + L0_parent_low_educ_lv + A0_ace",
           Y_death="Q.kplus1 ~ L0_male + L0_parent_low_educ_lv + L1 +
                    A0_ace * M_smoking")

# we define the formulas for the estimation of the 2 g function
gform <- c("A0_ace ~ L0_male + L0_parent_low_educ_lv",
           "M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1")

# the data frame should follow the time-ordering of the nodes
data_binary <- subset(df2_int, select = c(L0_male, L0_parent_low_educ_lv,
                                          A0_ace, L1,
                                          M_smoking, Y_death))



# choose a family of data-adaptive algorithms from the SuperLearner package
SL.library <- list(Q=c("SL.mean","SL.glm","SL.step.interaction","SL.xgboost"),
                   g=c("SL.mean","SL.glm","SL.step.interaction","SL.xgboost"))

set.seed(42)
## CDE, setting M=0
CDE_ltmle_M0_death <- ltmle(data = data_binary,
                            Anodes = c("A0_ace", "M_smoking"),
                            Lnodes = c("L1"), # intermediate confounders +/- baseline
                            Ynodes = c("Y_death"),
                            survivalOutcome = FALSE, # TRUE for time-to-event outcomes Y
                            Qform = Qform,
                            gform = gform,
                            abar = list(c(1,0), # counterfactual intervention do(A=1,M=0)
                                        c(0,0)), # counterfactual intervention do(A=0,M=0)
                            SL.library = SL.library,
                            estimate.time = FALSE, # estimate computation time?
                            gcomp = FALSE,
                            variance.method = "ic") # a more robust variance can
                                                    # be estimated with
                                                    # variance.method = "tmle"
summary(CDE_ltmle_M0_death)
# Parameter Estimate:  0.056766
#  Estimated Std Err:  0.018037
#            p-value:  0.0016488
#  95% Conf Interval: (0.021413, 0.092118)

## CDE, setting M=1
set.seed(42)
CDE_ltmle_M1_death <- ltmle(data = data_binary,
                            Anodes = c("A0_ace", "M_smoking"),
                            Lnodes = c("L1"), # intermediate confounders +/- baseline
                            Ynodes = c("Y_death"),
                            survivalOutcome = FALSE, # TRUE for time-to-event outcomes Y
                            Qform = Qform,
                            gform = gform,
                            abar = list(c(1,1), # counterfactual intervention do(A=1,M=1)
                                        c(0,1)), # counterfactual intervention do(A=0,M=1)
                            SL.library = SL.library,
                            estimate.time = FALSE, # estimate computation time?
                            gcomp = FALSE,
                            variance.method = "ic")
summary(CDE_ltmle_M1_death)
# Parameter Estimate:  0.094776
#  Estimated Std Err:  0.024
#            p-value:  7.8496e-05
#  95% Conf Interval: (0.047736, 0.14182)


# ---------------------------------------------------------------------------- #
## CDE for continuous outcomes ----
# ---------------------------------------------------------------------------- #
data_continuous <- subset(df2_int, select = c(L0_male, L0_parent_low_educ_lv,
                                              A0_ace, L1,
                                              M_smoking, Y_qol))

Qform <- c(L1="Q.kplus1 ~ L0_male + L0_parent_low_educ_lv + A0_ace",
           Y_qol="Q.kplus1 ~ L0_male + L0_parent_low_educ_lv + L1 +
                    A0_ace * M_smoking")

set.seed(42)
## CDE, setting M=0
CDE_ltmle_M0_qol <- ltmle(data = data_continuous,
                      Anodes = c("A0_ace", "M_smoking"),
                      Lnodes = c("L1"), # intermediate confounders +/- baseline confounders
                      Ynodes = c("Y_qol"),
                      survivalOutcome = FALSE, # TRUE for time-to-event outcomes Y
                      Qform = Qform,
                      gform = gform,
                      abar = list(c(1,0), # counterfactual intervention do(A=1,M=0)
                                  c(0,0)), # counterfactual intervention do(A=0,M=0)
                      SL.library = SL.library,
                      estimate.time = FALSE, # estimate computation time?
                      gcomp = FALSE,
                      variance.method = "ic")
summary(CDE_ltmle_M0_qol)
# Additive Treatment Effect:
#   Parameter Estimate:  -4.8023
#    Estimated Std Err:  0.43135
#              p-value:  <2e-16
#    95% Conf Interval: (-5.6477, -3.9569)

## CDE, setting M=1
set.seed(42)
CDE_ltmle_M1_qol <- ltmle(data = data_continuous,
                          Anodes = c("A0_ace", "M_smoking"),
                          Lnodes = c("L1"), # intermediate confounders +/- baseline
                          Ynodes = c("Y_qol"),
                          survivalOutcome = FALSE, # TRUE for time-to-event outcomes Y
                          Qform = Qform,
                          gform = gform,
                          abar = list(c(1,1), # counterfactual intervention do(A=1,M=1)
                                      c(0,1)), # counterfactual intervention do(A=0,M=1)
                          SL.library = SL.library,
                          estimate.time = FALSE, # estimate computation time?
                          gcomp = FALSE,
                          variance.method = "ic")
summary(CDE_ltmle_M1_qol)
# Additive Treatment Effect:
#   Parameter Estimate:  -10.219
#    Estimated Std Err:  0.544
#              p-value:  <2e-16
#    95% Conf Interval: (-11.285, -9.1523)


# ---------------------------------------------------------------------------- #
# Randomized Marginal Direct and Indirect effects ----
# ---------------------------------------------------------------------------- #
# remotes::install_github("nhejazi/medoutcon")
# remotes::install_github("nhejazi/medoutcon", INSTALL_opts=c("--no-multiarch"))
# https://arxiv.org/pdf/1912.09936
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")

library(tidyverse)
library(sl3)
library(medoutcon)
library(hal9001)

?medoutcon::medoutcon
## binary outcome ----
### one-step ----
### compute one-step estimate of the interventional direct effect
set.seed(1234)
os_de <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                   A = df2_int$A0_PM2.5, # numeric vector of the exposure
                   Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                   M = df2_int$M_diabetes, # numeric vector or matrix
                   Y = df2_int$Y_death, # numeric vector
                   effect = "direct",
                   b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                   estimator = "onestep")
os_de
# Interventional Direct Effect
# Estimator: onestep
# Estimate: 0.068
# Std. Error: 0.015
# 95% CI: [0.039, 0.096]

os_de$theta
# [1] 0.06763031

sqrt(os_de$var)
# [1] 0.01455125

set.seed(1234)
os_ie <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                   A = df2_int$A0_PM2.5, # numeric vector of the exposure
                   Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                   M = df2_int$M_diabetes, # numeric vector or matrix
                   Y = df2_int$Y_death, # numeric vector
                   effect = "indirect",
                   b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                   estimator = "onestep")
os_ie
# Interventional Indirect Effect
# Estimator: onestep
# Estimate: 0.015
# Std. Error: 0.004
# 95% CI: [0.008, 0.022]

os_ie$theta
# [1] 0.01502601

### Estimate counterfactuals E(Y_{a,G_{a*}})
set.seed(1234)
EY_1M0 <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")],
                    A = df2_int$A0_PM2.5,
                    Z = df2_int$L1,
                    M = df2_int$M_diabetes,
                    Y = df2_int$Y_death,
                    contrast = c(1,0),
                    b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                    estimator = "onestep")
EY_1M0
# Counterfactual TSM
# Contrast: A = 1, M(A = 0)
# Estimator: onestep
# Estimate: 0.273
# Std. Error: 0.014
# 95% CI: [0.246, 0.301]

EY_0M0 <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")],
                    A = df2_int$A0_PM2.5,
                    Z = df2_int$L1,
                    M = df2_int$M_diabetes,
                    Y = df2_int$Y_death,
                    contrast = c(0,0),
                    b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                    estimator = "onestep")
EY_0M0
# Counterfactual TSM
# Contrast: A = 0, M(A = 0)
# Estimator: onestep
# Estimate: 0.203
# Std. Error: 0.004
# 95% CI: [0.195, 0.212]

## randomized Natural Direct Effect =
rNDE <- EY_1M0$theta - EY_0M0$theta
# [1] 0.06964436

## se
se_rNDE <- sqrt(var(EY_1M0$eif - EY_0M0$eif) / nrow(df2_int))
# [1] 0.01459493

## 95%CI
c(rNDE - qnorm(0.975) * se_rNDE,
  rNDE + qnorm(0.975) * se_rNDE)
# [1] 0.04103882 0.09824989

### tmle ----
### compute tmle estimate of the interventional direct effect & indirect effects
set.seed(1234)
tmle_de <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                     A = df2_int$A0_PM2.5, # numeric vector of the exposure
                     Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                     M = df2_int$M_diabetes, # numeric vector or matrix
                     Y = df2_int$Y_death, # numeric vector
                     effect = "direct",
                     estimator = "tmle",
                     b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                     # estimator_args = list(cv_folds = 10L, # instead of 5L
                     #                       max_iter = 1L, # instead of 5L
                     #                       tiltmod_tol = 1) # instead of 5
                    )
tmle_de
# Interventional Direct Effect
# Estimator: tmle
# Estimate: 0.042             <- ?! biased (confirmed on simulations)
# Std. Error: 0.015
# 95% CI: [0.014, 0.071]

tmle_de$theta
# [1] 0.04222315

sqrt(tmle_de$var)
# [1] 0.01455125

set.seed(1234)
tmle_ie <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                   A = df2_int$A0_PM2.5, # numeric vector of the exposure
                   Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                   M = df2_int$M_diabetes, # numeric vector or matrix
                   Y = df2_int$Y_death, # numeric vector
                   effect = "indirect",
                   b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                   estimator = "tmle")
tmle_ie
# Interventional Indirect Effect
# Estimator: tmle
# Estimate: 0.013
# Std. Error: 0.004
# 95% CI: [0.006, 0.021]


## continuous outcome ----
### one-step ----
### compute one-step estimate of the interventional direct effect
set.seed(1234)
os_de <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                   A = df2_int$A0_PM2.5, # numeric vector of the exposure
                   Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                   M = df2_int$M_diabetes, # numeric vector or matrix
                   Y = df2_int$Y_qol, # numeric vector
                   effect = "direct",
                   b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                   estimator = "onestep")
os_de
# Interventional Direct Effect
# Estimator: onestep
# Estimate: -6.597
# Std. Error: 0.346
# 95% CI: [-7.276, -5.918]

set.seed(1234)
os_ie <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                   A = df2_int$A0_PM2.5, # numeric vector of the exposure
                   Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                   M = df2_int$M_diabetes, # numeric vector or matrix
                   Y = df2_int$Y_qol, # numeric vector
                   effect = "indirect",
                   b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                   estimator = "onestep")
os_ie
# Interventional Indirect Effect
# Estimator: onestep
# Estimate: -1.703
# Std. Error: 0.239
# 95% CI: [-2.172, -1.234]

### tmle ----
### compute tmle estimate of the interventional direct effect & indirect effects
set.seed(1234)
tmle_de <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                     A = df2_int$A0_PM2.5, # numeric vector of the exposure
                     Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                     M = df2_int$M_diabetes, # numeric vector or matrix
                     Y = df2_int$Y_qol, # numeric vector
                     effect = "direct",
                     b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                     estimator = "tmle")
tmle_de
# Interventional Direct Effect
# Estimator: tmle
# Estimate: -4.808          <- also biased ?
# Std. Error: 0.346
# 95% CI: [-5.487, -4.129]

set.seed(1234)
tmle_ie <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                     A = df2_int$A0_PM2.5, # numeric vector of the exposure
                     Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                     M = df2_int$M_diabetes, # numeric vector or matrix
                     Y = df2_int$Y_qol, # numeric vector
                     effect = "indirect",
                     b_learners = sl3::Lrnr_hal9001$new(), # outcome regression
                     estimator = "tmle")
tmle_ie
# Interventional Indirect Effect
# Estimator: tmle
# Estimate: -1.738
# Std. Error: 0.239
# 95% CI: [-2.207, -1.269]


### checking the biased estimate by tmle ----
set.seed(1234)
tmle_de <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                     A = df2_int$A0_PM2.5, # numeric vector of the exposure
                     Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                     M = df2_int$M_diabetes, # numeric vector or matrix
                     Y = df2_int$Y_death, # numeric vector
                     effect = "direct",
                     g_learners = sl3::Lrnr_hal9001$new(),
                     h_learners = sl3::Lrnr_hal9001$new(),
                     b_learners = sl3::Lrnr_hal9001$new(),
                     q_learners = sl3::Lrnr_hal9001$new(),
                     r_learners = sl3::Lrnr_hal9001$new(),
                     u_learners = sl3::Lrnr_hal9001$new(),
                     v_learners = sl3::Lrnr_hal9001$new(),
                     d_learners = sl3::Lrnr_hal9001$new(),
                     estimator = "tmle")
tmle_de
# Interventional Direct Effect
# Estimator: tmle
# Estimate: 0.043     # calculation is longer but still biased ...
# Std. Error: 0.015
# 95% CI: [0.015, 0.072]


### Example from N Hejazi github
# produces a simple data set based on ca causal model with mediation
library(data.table)
set.seed(1584)
make_example_data <- function(n_obs = 1000) {
  ## baseline covariates
  w_1 <- rbinom(n_obs, 1, prob = 0.6)
  w_2 <- rbinom(n_obs, 1, prob = 0.3)
  w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
  w <- cbind(w_1, w_2, w_3)
  w_names <- paste("W", seq_len(ncol(w)), sep = "_")

  ## exposure
  a <- as.numeric(rbinom(n_obs, 1, plogis(rowSums(w) - 2)))

  ## mediator-outcome confounder affected by treatment
  z <- rbinom(n_obs, 1, plogis(rowMeans(-log(2) + w - a) + 0.2))

  ## mediator -- could be multivariate
  m <- rbinom(n_obs, 1, plogis(rowSums(log(3) * w[, -3] + a - z)))
  m_names <- "M"

  ## outcome
  y <- rbinom(n_obs, 1, plogis(1 / (rowSums(w) - z + a + m)))

  ## construct output
  dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
  setnames(dat, c(w_names, "A", "Z", m_names, "Y"))
  return(dat)
}

# set seed and simulate example data
example_data <- make_example_data()
library(stringr)
w_names <- str_subset(colnames(example_data), "W")
m_names <- str_subset(colnames(example_data), "M")

# quick look at the data
head(example_data)
#>    W_1 W_2 W_3 A Z M Y
#> 1:   1   0   1 0 0 0 1
#> 2:   0   1   0 0 0 1 0
#> 3:   1   1   1 1 0 1 1
#> 4:   0   1   1 0 0 1 0
#> 5:   0   0   0 0 0 1 1
#> 6:   1   0   1 1 0 1 0

# compute one-step estimate of the interventional direct effect
os_de <- medoutcon(W = example_data[, ..w_names],
                   A = example_data$A,
                   Z = example_data$Z,
                   M = example_data[, ..m_names],
                   Y = example_data$Y,
                   effect = "direct",
                   estimator = "onestep")
os_de
#> Interventional Direct Effect
#> Estimator: onestep
#> Estimate: -0.065
#> Std. Error: 0.054
#> 95% CI: [-0.17, 0.041]

# compute targeted minimum loss estimate of the interventional direct effect
tmle_de <- medoutcon(W = example_data[, ..w_names],
                     A = example_data$A,
                     Z = example_data$Z,
                     M = example_data[, ..m_names],
                     Y = example_data$Y,
                     effect = "direct",
                     estimator = "tmle")
tmle_de
#> Interventional Direct Effect
#> Estimator: tmle
#> Estimate: -0.06
#> Std. Error: 0.058
#> 95% CI: [-0.173, 0.053]



### test on simulations for df2 ----
# ---------------------------------------------------------------------------- #
rm(list=ls())
library(medoutcon)
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
                                      c("MRDE.death.os","MRDE.death.os.cov","MRIE.death.os","MRIE.death.os.cov",
                                        "MRDE.death.tmle","MRDE.death.tmle.cov","MRIE.death.tmle","MRIE.death.tmle.cov",
                                        "MRDE.qol.os","MRDE.qol.os.cov","MRIE.qol.os","MRIE.qol.os.cov",
                                        "MRDE.qol.tmle","MRDE.qol.tmle.cov","MRIE.qol.tmle","MRIE.qol.tmle.cov")))
set.seed(54321)
for(i in 1:1000) {
  print(paste0("simulation ",i))
  # generate data
  df2_int <- gen.data.causal.model.2(N = 10000, A.M.inter = 1)

  ### algorithm for CRDE and CRIE
  # death
  os_de_death <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                           A = df2_int$A0_PM2.5, # numeric vector of the exposure
                           Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                           M = df2_int$M_diabetes, # numeric vector or matrix
                           Y = df2_int$Y_death, # numeric vector
                           effect = "direct",
                           estimator = "onestep")
  os_ie_death <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                           A = df2_int$A0_PM2.5, # numeric vector of the exposure
                           Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                           M = df2_int$M_diabetes, # numeric vector or matrix
                           Y = df2_int$Y_death, # numeric vector
                           effect = "indirect",
                           estimator = "onestep")

  tmle_de_death <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                             A = df2_int$A0_PM2.5, # numeric vector of the exposure
                             Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                             M = df2_int$M_diabetes, # numeric vector or matrix
                             Y = df2_int$Y_death, # numeric vector
                             effect = "direct",
                             estimator = "tmle")
  tmle_ie_death <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                             A = df2_int$A0_PM2.5, # numeric vector of the exposure
                             Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                             M = df2_int$M_diabetes, # numeric vector or matrix
                             Y = df2_int$Y_death, # numeric vector
                             effect = "indirect",
                             estimator = "tmle")

  # qol
  os_de_qol <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                         A = df2_int$A0_PM2.5, # numeric vector of the exposure
                         Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                         M = df2_int$M_diabetes, # numeric vector or matrix
                         Y = df2_int$Y_qol, # numeric vector
                         effect = "direct",
                         estimator = "onestep")
  os_ie_qol <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                         A = df2_int$A0_PM2.5, # numeric vector of the exposure
                         Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                         M = df2_int$M_diabetes, # numeric vector or matrix
                         Y = df2_int$Y_qol, # numeric vector
                         effect = "indirect",
                         estimator = "onestep")

  tmle_de_qol <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                           A = df2_int$A0_PM2.5, # numeric vector of the exposure
                           Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                           M = df2_int$M_diabetes, # numeric vector or matrix
                           Y = df2_int$Y_qol, # numeric vector
                           effect = "direct",
                           estimator = "tmle")
  tmle_ie_qol <- medoutcon(W = df2_int[,c("L0_male","L0_soc_env")], #matrix of baseline L(0)
                           A = df2_int$A0_PM2.5, # numeric vector of the exposure
                           Z = df2_int$L1, # numeric vector L(1) (only 1 variable)
                           M = df2_int$M_diabetes, # numeric vector or matrix
                           Y = df2_int$Y_qol, # numeric vector
                           effect = "indirect",
                           estimator = "tmle")

  ### save results
  results.df2[i,"MRDE.death.os"] <- os_de_death$theta - (0.073882)
  results.df2[i,"MRDE.death.os.cov"] <- (((os_de_death$theta - qnorm(0.975) * sqrt(os_de_death$var)) < 0.073882) &
                                           ((os_de_death$theta + qnorm(0.975) * sqrt(os_de_death$var)) > 0.073882))
  results.df2[i,"MRIE.death.os"] <- os_ie_death$theta - (0.0154)
  results.df2[i,"MRIE.death.os.cov"] <- (((os_ie_death$theta - qnorm(0.975) * sqrt(os_ie_death$var)) < 0.0154) &
                                           ((os_ie_death$theta + qnorm(0.975) * sqrt(os_ie_death$var)) > 0.0154))

  results.df2[i,"MRDE.death.tmle"] <- tmle_de_death$theta - (0.073882)
  results.df2[i,"MRDE.death.tmle.cov"] <- (((tmle_de_death$theta - qnorm(0.975) * sqrt(tmle_de_death$var)) < 0.073882) &
                                             ((tmle_de_death$theta + qnorm(0.975) * sqrt(tmle_de_death$var)) > 0.073882))
  results.df2[i,"MRIE.death.tmle"] <- tmle_ie_death$theta - (0.0154)
  results.df2[i,"MRIE.death.tmle.cov"] <- (((tmle_ie_death$theta - qnorm(0.975) * sqrt(tmle_ie_death$var)) < 0.0154) &
                                             ((tmle_ie_death$theta + qnorm(0.975) * sqrt(tmle_ie_death$var)) > 0.0154))

  results.df2[i,"MRDE.qol.os"] <- os_de_qol$theta - (-6.647)
  results.df2[i,"MRDE.qol.os.cov"] <- (((os_de_qol$theta - qnorm(0.975) * sqrt(os_de_qol$var)) < -6.647) &
                                         ((os_de_qol$theta + qnorm(0.975) * sqrt(os_de_qol$var)) > -6.647))
  results.df2[i,"MRIE.qol.os"] <- os_ie_qol$theta - (-1.96)
  results.df2[i,"MRIE.qol.os.cov"] <- (((os_ie_qol$theta - qnorm(0.975) * sqrt(os_ie_qol$var)) < -1.96) &
                                         ((os_ie_qol$theta + qnorm(0.975) * sqrt(os_ie_qol$var)) > -1.96))

  results.df2[i,"MRDE.qol.tmle"] <- tmle_de_qol$theta - (-6.647)
  results.df2[i,"MRDE.qol.tmle.cov"] <- (((tmle_de_qol$theta - qnorm(0.975) * sqrt(tmle_de_qol$var)) < -6.647) &
                                           ((tmle_de_qol$theta + qnorm(0.975) * sqrt(tmle_de_qol$var)) > -6.647))
  results.df2[i,"MRIE.qol.tmle"] <- tmle_ie_qol$theta - (-1.96)
  results.df2[i,"MRIE.qol.tmle.cov"] <- (((tmle_ie_qol$theta - qnorm(0.975) * sqrt(tmle_ie_qol$var)) < -1.96) &
                                           ((tmle_ie_qol$theta + qnorm(0.975) * sqrt(tmle_ie_qol$var)) > -1.96))
}

sapply(data.frame(results.df2) , mean)
#  MRDE.death.os   MRDE.death.os.cov       MRIE.death.os   MRIE.death.os.cov
#   0.0003613417        0.9520000000        0.0001405630        0.9390000000 => OK pour one step
# MRDE.death.tmle MRDE.death.tmle.cov     MRIE.death.tmle MRIE.death.tmle.cov
#   -0.0091549283        0.8480000000        0.0001887599        0.9440000000 => biased for the direct effect, OK of Indirect effect
#    MRDE.qol.os     MRDE.qol.os.cov         MRIE.qol.os     MRIE.qol.os.cov
#   0.0147858126        0.9550000000        0.0053678236        0.8610000000  => indirect effect not biased, but coverage of the one step not-optimal !
#  MRDE.qol.tmle   MRDE.qol.tmle.cov       MRIE.qol.tmle   MRIE.qol.tmle.cov
#   0.7043470959        0.5960000000        0.0053361064        0.8610000000 => biased for the direct effect and coverage not optimal for both


boxplot(subset(data.frame(results.df2), select = c("MRDE.death.os", "MRDE.death.tmle"))) # tmle biased
boxplot(subset(data.frame(results.df2), select = c("MRIE.death.os","MRIE.death.tmle"))) # both are very close
boxplot(subset(data.frame(results.df2), select = c("MRDE.qol.os", "MRDE.qol.tmle"))) # tmle biased
boxplot(subset(data.frame(results.df2), select = c("MRIE.qol.os","MRIE.qol.tmle"))) # both are very close

results.df2.rel <- results.df2
results.df2.rel[,c("MRDE.death.os","MRDE.death.tmle")] <- results.df2.rel[,c("MRDE.death.os","MRDE.death.tmle")] / (0.073882)
results.df2.rel[,c("MRIE.death.os","MRIE.death.tmle")] <- results.df2.rel[,c("MRIE.death.os","MRIE.death.tmle")] / (0.0154)
results.df2.rel[,c("MRDE.qol.os","MRDE.qol.tmle")] <- results.df2.rel[,c("MRDE.qol.os","MRDE.qol.tmle")] / (-6.647)
results.df2.rel[,c("MRIE.qol.os","MRIE.qol.tmle")] <- results.df2.rel[,c("MRIE.qol.os","MRIE.qol.tmle")] / (-1.96)
sapply(data.frame(results.df2.rel) , mean)
#    MRDE.death.os   MRDE.death.os.cov       MRIE.death.os   MRIE.death.os.cov
#      0.004890795         0.952000000         0.009127468         0.939000000 < 1% biased for one step
#  MRDE.death.tmle MRDE.death.tmle.cov     MRIE.death.tmle MRIE.death.tmle.cov
#     -0.123912839         0.848000000         0.012257135         0.944000000 > 10% biaised for the direct effect and 1.2% for indirect with tmle
#      MRDE.qol.os     MRDE.qol.os.cov         MRIE.qol.os     MRIE.qol.os.cov
#     -0.002224434         0.955000000        -0.002738686         0.861000000 < 1% biased for one step
#    MRDE.qol.tmle   MRDE.qol.tmle.cov       MRIE.qol.tmle   MRIE.qol.tmle.cov
#     -0.105964660         0.596000000        -0.002722503         0.861000000



# ---------------------------------------------------------------------------- #
# Randomized Conditional Direct and Indirect effects ----
# ---------------------------------------------------------------------------- #
# To do ...
