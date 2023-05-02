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
Qbar_Ais1 <- predict(Q.fit, newdata = data.A1, type = "response")
logitQ <- predict(Q.fit, newdata = data.A1, type = "link")

## 2) Estimate the treatment mechanism
g.L <- glm(A0_ace ~ L0_male + L0_parent_low_educ_lv,
           family = "binomial", data = df2_int)

# predict the probabilities P(A0_ace=1|L(0)) & P(A0_ace=0|L(0))
pred.g1.L <- predict(g.L, type="response")
pred.g0.L <- 1 - pred.g1.L

# the predicted probability of the observed treatment A_i=a is :
gA.L <- rep(NA, nrow(df2_int))
gA.L[df2_int$A0_ace==1] <- pred.g1.L[df2_int$A0_ace==1]
gA.L[df2_int$A0_ace==0] <- pred.g0.L[df2_int$A0_ace==0]

# its useful to check the distribution of gA.L, as values close to 0 or 1
# (< 0.01 or > 0.999)
# are indicators of near positivity violation
summary(gA.L)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.06109 0.84370 0.89011 0.80279 0.91106 0.93891
# there is no positivity issues in this example.
head(gA.L, 10)
# [1] 0.8901078 0.1562975 0.8437025 0.9110593 0.8437025 0.8437025 0.8901078 0.8901078 0.1562975 0.8437025

plot(gA.L, Psi_Ais1$cum.g)
plot(gA.L[Psi_Ais1$cum.g.used == TRUE], Psi_Ais1$cum.g[Psi_Ais1$cum.g.used == TRUE])
# en fait, on a pas besoin de calculer la probabilité d'observer sa propre observation, avec l'indicatrice, on retombe sur nos pieds !


## 3) Determine a parametric family of fluctuations of Qbar,
#     defined by parametric regression including a clever covarate chosen so the
#     loss function L(Q(epsilon)) holds with epsilon playing the role of the coefficient
#     in front of the clever covariate.

# The clever covariate H(A,L(0)) depends on g:
H <- (df2_int$A0_ace == 1) / gA.L


## 3) Update the initial fit Qbar from step 1.
# This is achieved by holding Qbar fixed (as intercept) while estimating the
# coefficient epsilon for H

# for example:
update.fit <- glm(df2_int$Y_death ~ -1 + offset(logitQ) + H,
                  family = "quasibinomial")
# Coefficients:
#   H
# -0.0001756
# Degrees of Freedom: 1124 Total (i.e. Null);  1123 Residual
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
