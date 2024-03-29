### test programs for estimations based on TMLE

rm(list=ls())

df2_int <- read.csv(file = "data/df2_int.csv")

################################################################################
######################### Estimation of the Average Total Effect (ATE)
################################################################################
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


########################
### manual calculation:
########################
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


################################################################################
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



## for continuous outcomes
#############################
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
