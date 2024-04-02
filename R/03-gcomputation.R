
### test programs for estimations based on G-computation

rm(list=ls())

df2 <- read.csv(file = "data/df2.csv")
df2_int <- read.csv(file = "data/df2_int.csv")


# ---------------------------------------------------------------------------- #
# I) Estimation of the Average Total Effect (ATE) ------------------------------
# ---------------------------------------------------------------------------- #

# 1. Estimate Qbar
Q.tot.death <- glm(Y_death ~ A0_ace + L0_male + L0_parent_low_educ_lv, family = "binomial", data = df2_int)
Q.tot.qol <- glm(Y_qol ~ A0_ace + L0_male + L0_parent_low_educ_lv, family = "gaussian", data = df2_int)

# 2. Predict an outcome for each subject, setting A=0 and A=1
# prepare data sets used to predict the outcome under the counterfactual scenarios setting A=0 and A=1
data.A1 <- data.A0 <- df2_int
data.A1$A0_ace <- 1
data.A0$A0_ace <- 0

# predict values
Y1.death.pred <- predict(Q.tot.death, newdata = data.A1, type = "response")
Y0.death.pred <- predict(Q.tot.death, newdata = data.A0, type = "response")

Y1.qol.pred <- predict(Q.tot.qol, newdata = data.A1, type = "response")
Y0.qol.pred <- predict(Q.tot.qol, newdata = data.A0, type = "response")

# 3. Plug the predicted outcome in the gformula and use the sample mean to estimate the ATE
ATE.death.gcomp <- mean(Y1.death.pred - Y0.death.pred)
ATE.death.gcomp
# [1] 0.08270821

ATE.qol.gcomp <- mean(Y1.qol.pred - Y0.qol.pred)
ATE.qol.gcomp
# [1] -8.360691


### 95% CI calculation applying a bootstrap procedure

set.seed(1234)
B <- 2000
bootstrap.estimates <- data.frame(matrix(NA, nrow = B, ncol = 2))
colnames(bootstrap.estimates) <- c("boot.death.est", "boot.qol.est")
for (b in 1:B){
  # sample the indices 1 to n with replacement
  bootIndices <- sample(1:nrow(df2_int), replace=T)
  bootData <- df2_int[bootIndices,]

  if ( round(b/100, 0) == b/100 ) print(paste0("bootstrap number ",b))

  Q.tot.death <- glm(Y_death ~ A0_ace + L0_male + L0_parent_low_educ_lv, family = "binomial", data = bootData)
  Q.tot.qol <- glm(Y_qol ~ A0_ace + L0_male + L0_parent_low_educ_lv, family = "gaussian", data = bootData)

  boot.A.1 <- boot.A.0 <- bootData
  boot.A.1$A0_ace <- 1
  boot.A.0$A0_ace <- 0

  Y1.death.boot <- predict(Q.tot.death, newdata = boot.A.1, type = "response")
  Y0.death.boot <- predict(Q.tot.death, newdata = boot.A.0, type = "response")

  Y1.qol.boot <- predict(Q.tot.qol, newdata = boot.A.1, type = "response")
  Y0.qol.boot <- predict(Q.tot.qol, newdata = boot.A.0, type = "response")

  bootstrap.estimates[b,"boot.death.est"] <- mean(Y1.death.boot - Y0.death.boot)
  bootstrap.estimates[b,"boot.qol.est"] <- mean(Y1.qol.boot - Y0.qol.boot)
}

IC95.ATE.death <- c(ATE.death.gcomp - qnorm(0.975)*sd(bootstrap.estimates[,"boot.death.est"]),
                    ATE.death.gcomp + qnorm(0.975)*sd(bootstrap.estimates[,"boot.death.est"]) )
IC95.ATE.death
# [1] 0.05571017 0.10970624

IC95.ATE.qol <- c(ATE.qol.gcomp - qnorm(0.975)*sd(bootstrap.estimates[,"boot.qol.est"]),
                  ATE.qol.gcomp + qnorm(0.975)*sd(bootstrap.estimates[,"boot.qol.est"]) )
IC95.ATE.qol
# [1] -9.156051 -7.565331



# ---------------------------------------------------------------------------- #
# II) Estimation of the Control direct Effects (CDE) ---------------------------
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
## II.1) parametric g-computation ----------------------------------------------
# ---------------------------------------------------------------------------- #
rm(list=ls())
df2_int <- read.csv(file = "df2_int.csv")


# a. fit parametric models to estimate the density of intermediate confounders
L1.model <- glm(L1 ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "binomial", data = df2_int)

# b. fit parametric models for the outcome
Y2.death.model <- glm(Y_death ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 + M_smoking + A0_ace:M_smoking, family = "binomial", data = df2_int)
Y2.qol.model <- glm(Y_qol ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 + M_smoking + A0_ace:M_smoking, family = "gaussian", data = df2_int)

# c. simulate L1 values under the counterfactual scenarios setting A0=0 or A0=1
set.seed(54321)
data.A0  <- data.A1 <- df2_int
data.A0$A0_ace <- 0
data.A1$A0_ace <- 1
p.L1.A0 <- predict(L1.model, newdata = data.A0, type="response")
p.L1.A1 <- predict(L1.model, newdata = data.A1, type="response")
sim.L1.A0 <- rbinom(n = nrow(df2_int), size = 1, prob = p.L1.A0)
sim.L1.A1 <- rbinom(n = nrow(df2_int), size = 1, prob = p.L1.A1)

# d. estimate mean outcomes under the counterfactual scenarios setting different levels of exposures for A and M:
# {A=0, M=0} or {A=1, M=0} or {A=0, M=1} or {A=1, M=1}

data.A0.M0 <- data.A0.M1 <- data.A0
data.A1.M0 <- data.A1.M1 <- data.A1

# L1 variable is replaced by the simulated values in step c)
data.A0.M0$L1 <- sim.L1.A0
data.A0.M1$L1 <- sim.L1.A0
data.A1.M0$L1 <- sim.L1.A1
data.A1.M1$L1 <- sim.L1.A1

# set M to 0 or 1
data.A0.M0$M_smoking <- 0
data.A0.M1$M_smoking <- 1
data.A1.M0$M_smoking <- 0
data.A1.M1$M_smoking <- 1

p.death.A0.M0 <- predict(Y2.death.model, newdata = data.A0.M0, type="response")
p.death.A1.M0 <- predict(Y2.death.model, newdata = data.A1.M0, type="response")
p.death.A0.M1 <- predict(Y2.death.model, newdata = data.A0.M1, type="response")
p.death.A1.M1 <- predict(Y2.death.model, newdata = data.A1.M1, type="response")

m.qol.A0.M0 <- predict(Y2.qol.model, newdata = data.A0.M0, type="response")
m.qol.A1.M0 <- predict(Y2.qol.model, newdata = data.A1.M0, type="response")
m.qol.A0.M1 <- predict(Y2.qol.model, newdata = data.A0.M1, type="response")
m.qol.A1.M1 <- predict(Y2.qol.model, newdata = data.A1.M1, type="response")

# e. Estimate CDE
# CDE setting M=0
CDE.death.m0.gcomp.param <- mean(p.death.A1.M0) - mean(p.death.A0.M0)
CDE.death.m0.gcomp.param
# [1] 0.06289087

CDE.qol.m0.gcomp.param <- mean(m.qol.A1.M0) - mean(m.qol.A0.M0)
CDE.qol.m0.gcomp.param
# [1] -4.838654

# CDE setting M=1
CDE.death.m1.gcomp.param <- mean(p.death.A1.M1) - mean(p.death.A0.M1)
CDE.death.m1.gcomp.param
# [1] 0.08751016

CDE.qol.m1.gcomp.param <- mean(m.qol.A1.M1) - mean(m.qol.A0.M1)
CDE.qol.m1.gcomp.param
# [1] -10.35059


# ---------------------------------------------------------------------------- #
## II.2) G-computation by iterative conditional expectation --------------------
# ---------------------------------------------------------------------------- #
rm(list=ls())
df2_int <- read.csv(file = "df2_int.csv")

# 1) Regress the outcome on L0, A, L1 and M (and the A*M interaction if appropriate)
Y.death.model <- glm(Y_death ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 +
                       M_smoking + A0_ace:M_smoking,
                     family = "binomial", data = df2_int)
Y.qol.model <- glm(Y_qol ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 +
                     M_smoking + A0_ace:M_smoking,
                   family = "gaussian", data = df2_int)

# 2) Generate predicted values by evaluating the regression setting the mediator
#    value to M=0 or to M=1
#    (Note: it is also possible to set A=0 or A=1 to evaluate the regression at
#     exposure history of interest: {A0=1,M=0},{A0=0,M=0},{A0=1,M=1},{A0=0,M=1})
data.Mis0 <- data.Mis1 <- df2_int
data.Mis0$M_smoking <- 0
data.Mis1$M_smoking <- 1

Q.L2.death.Mis0 <- predict(Y.death.model, newdata = data.Mis0, type="response")
Q.L2.death.Mis1 <- predict(Y.death.model, newdata = data.Mis1, type="response")

Q.L2.qol.Mis0 <- predict(Y.qol.model, newdata = data.Mis0, type="response")
Q.L2.qol.Mis1 <- predict(Y.qol.model, newdata = data.Mis1, type="response")


# 3) Regress the predicted values conditional on the exposure A
#    and baseline confounders L(0)
L1.death.Mis0.model <- glm(Q.L2.death.Mis0 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                           family = "quasibinomial", data = df2_int)
L1.death.Mis1.model <- glm(Q.L2.death.Mis1 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                           family = "quasibinomial", data = df2_int)

L1.qol.Mis0.model <- glm(Q.L2.qol.Mis0 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                         family = "gaussian", data = df2_int)
L1.qol.Mis1.model <- glm(Q.L2.qol.Mis1 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                         family = "gaussian", data = df2_int)

# 4) generate predicted values by evaluating the regression at exposure
#    of interest: {A=1} & {A=0}
data.Ais0 <- data.Ais1 <- df2_int
data.Ais0$A0_ace <- 0
data.Ais1$A0_ace <- 1

Q.L1.death.Ais0.Mis0 <- predict(L1.death.Mis0.model, newdata = data.Ais0, type="response")
Q.L1.death.Ais1.Mis0 <- predict(L1.death.Mis0.model, newdata = data.Ais1, type="response")
Q.L1.death.Ais0.Mis1 <- predict(L1.death.Mis1.model, newdata = data.Ais0, type="response")
Q.L1.death.Ais1.Mis1 <- predict(L1.death.Mis1.model, newdata = data.Ais1, type="response")

Q.L1.qol.Ais0.Mis0 <- predict(L1.qol.Mis0.model, newdata = data.Ais0, type="response")
Q.L1.qol.Ais1.Mis0 <- predict(L1.qol.Mis0.model, newdata = data.Ais1, type="response")
Q.L1.qol.Ais0.Mis1 <- predict(L1.qol.Mis1.model, newdata = data.Ais0, type="response")
Q.L1.qol.Ais1.Mis1 <- predict(L1.qol.Mis1.model, newdata = data.Ais1, type="response")

# 5) Take empirical mean of final predicted outcomes to estimate CDE
# CDE setting M=0
CDE.death.m0.gcomp.ice <- mean(Q.L1.death.Ais1.Mis0) - mean(Q.L1.death.Ais0.Mis0)
CDE.death.m0.gcomp.ice
# [1] 0.06341297

CDE.qol.m0.gcomp.ice <- mean(Q.L1.qol.Ais1.Mis0) - mean(Q.L1.qol.Ais0.Mis0)
CDE.qol.m0.gcomp.ice
# [1] -4.869509

# CDE setting M=1
CDE.death.m1.gcomp.ice <- mean(Q.L1.death.Ais1.Mis1) - mean(Q.L1.death.Ais0.Mis1)
CDE.death.m1.gcomp.ice
# [1] 0.08810508

CDE.qol.m1.gcomp.ice <- mean(Q.L1.qol.Ais1.Mis1) - mean(Q.L1.qol.Ais0.Mis1)
CDE.qol.m1.gcomp.ice
# [1] -10.38144



# ---------------------------------------------------------------------------- #
## II.3) Sequential g-estimator ------------------------------------------------
# ---------------------------------------------------------------------------- #
# Approach described for continuous outcome.
# Extension for binary outcomes using OR, in case-control studies, is described in Vansteelandt et al. Epidemiology 20(6);2009.
rm(list=ls())
df2_int <- read.csv(file = "df2_int.csv")

# 1) Regress the outcome on past
Y.qol.model <- glm(Y_qol ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 +
                     M_smoking + A0_ace:M_smoking,
                   family = "gaussian", data = df2_int)

# 2) Calculate a residual outcome Y - (coef.M * M_smoking) - (coef.A0:M * A0:M)
Y.res <- (df2_int$Y_qol -
            (Y.qol.model$coefficients["M_smoking"] * df2_int$M_smoking) -
            (Y.qol.model$coefficients["A0_ace:M_smoking"] * df2_int$A0_ace *
               df2_int$M_smoking) )

# 3) Regress the residual outcome on the exposure A and baseline confounders L(0)
Y.res.model <- glm(Y.res ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                   family = "gaussian", data = df2_int)

# 4) Use coefficients estimated from the 1st and 2nd regression to estimate CDE:
CDE.qol.m0.seq <- Y.res.model$coefficients["A0_ace"] + 0*Y.qol.model$coefficients["A0_ace:M_smoking"]
CDE.qol.m0.seq
# -4.869509

CDE.qol.m1.seq <- Y.res.model$coefficients["A0_ace"] + 1*Y.qol.model$coefficients["A0_ace:M_smoking"]
CDE.qol.m1.seq
# -10.38144


# ---------------------------------------------------------------------------- #
# III) Estimation of Marginal Randomized Direct and Indirect Effects -----------
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# III.1) parametric g-computation ----------------------------------------------
# ---------------------------------------------------------------------------- #
# described in Lin et al. Epidemiology 2017;28(2):266-74
# The approach is described as an adaptation of the parametric g-computation presented for controlled direct effects.

# The following steps are applied:
# 1) Fit parametric models for the observed data for the time-varying confounders L(1), the mediator M and the outcome Y

# 2) Estimate the joint distribution of time-varying


set.seed(54321)

# steps 1) to 3) will be repeated some fixed number k (for example k=25)
# we will save the k results in a matrix of k rows and 4 columns for the randomized
# direct and indirect effects on death (binary) and QoL (continuous outcome)
est <- matrix(NA, nrow = 25, ncol = 4)
colnames(est) <- c("rNDE.death", "rNIE.death", "rNDE.qol", "rNIE.qol")

# repeat k=25 times the following steps 1) to 3)
for (k in 1:25) {
  ## 1) Fit parametric models for the observed data for the time-varying
  ##    confounders L(1), the mediator M and the outcome Y
  ### 1a) fit parametric models of the confounders and mediators given the past
  L1.model <- glm(L1 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                  family = "binomial", data = df2_int)
  M.model <- glm(M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1,
                 family = "binomial", data = df2_int)
  ### 1b) fit parametric models of the outcomes given the past
  Y.death.model <- glm(Y_death ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 +
                         M_smoking + A0_ace:M_smoking,
                       family = "binomial", data = df2_int)
  Y.qol.model <- glm(Y_qol ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 +
                       M_smoking + A0_ace:M_smoking,
                     family = "gaussian", data = df2_int)

  ## 2) Estimate the joint distribution of time-varying confounders and of the
  ##    mediator under the counterfactual scenarios setting A0_ace = 1 or 0

  # set the exposure A0_ace to 0 or 1 in two new counterfactual data sets
  data.A0  <- data.A1 <- df2_int
  data.A0$A0_ace <- 0
  data.A1$A0_ace <- 1
  # simulate L1 values under the counterfactual exposures A0_ace=0 or A0_ace=1
  p.L1.A0 <- predict(L1.model, newdata = data.A0, type="response")
  p.L1.A1 <- predict(L1.model, newdata = data.A1, type="response")
  sim.L1.A0 <- rbinom(n = nrow(df2_int), size = 1, prob = p.L1.A0)
  sim.L1.A1 <- rbinom(n = nrow(df2_int), size = 1, prob = p.L1.A1)

  # replace L(1) by their counterfactual values in the data under A=0 or A=1
  data.A0.L <- data.A0
  data.A1.L <- data.A1
  data.A0.L$L1 <- sim.L1.A0
  data.A1.L$L1 <- sim.L1.A1

  # simulate M values under the counterfactual exposures A0_ace=0 or A0_ace=1
  p.M.A0 <- predict(M.model, newdata = data.A0.L, type="response")
  p.M.A1 <- predict(M.model, newdata = data.A1.L, type="response")
  sim.M.A0 <- rbinom(n = nrow(df2_int), size = 1, prob = p.M.A0)
  sim.M.A1 <- rbinom(n = nrow(df2_int), size = 1, prob = p.M.A1)
  # permute the n values of the joint mediator to obtain the random distributions
  # of the mediator: G_{A=0} and G_{A=1}
  marg.M.A0 <- sample(sim.M.A0, replace = FALSE)
  marg.M.A1 <- sample(sim.M.A1, replace = FALSE)

  ## 3) Simulate the outcomes Y_{A=0,G_{A=0}}
  ### 3a) use the previous permutation to replace the mediator
  ### in the counterfactual data sets for Y_{A=0,G_{A=0}}, Y_{A=1,G_{A=1}} and
  ### Y_{A=1,G_{A=0}}
  data.A0.G0 <- data.A0.G1 <- data.A0.L
  data.A1.G0 <- data.A1.G1 <- data.A1.L

  data.A0.G0$M_smoking <- marg.M.A0
  # data.A0.G1$M_smoking <- marg.M.A1 # note: this data set will not be useful

  data.A1.G0$M_smoking <- marg.M.A0
  data.A1.G1$M_smoking <- marg.M.A1

  # simulate the average outcome using the models fitted at step 1)
  p.death.A1.G1 <- predict(Y.death.model, newdata = data.A1.G1, type="response")
  p.death.A1.G0 <- predict(Y.death.model, newdata = data.A1.G0, type="response")
  p.death.A0.G0 <- predict(Y.death.model, newdata = data.A0.G0, type="response")

  m.qol.A1.G1 <- predict(Y.qol.model, newdata = data.A1.G1, type="response")
  m.qol.A1.G0 <- predict(Y.qol.model, newdata = data.A1.G0, type="response")
  m.qol.A0.G0 <- predict(Y.qol.model, newdata = data.A0.G0, type="response")

  ## save the results in row k
  # rNDE = E(Y_{A=1,G_{A=0}}) - E(Y_{A=0,G_{A=0}})
  # rNIE = E(Y_{A=1,G_{A=1}}) - E(Y_{A=1,G_{A=0}})
  est[k,"rNDE.death"]<- mean(p.death.A1.G0) - mean(p.death.A0.G0)
  est[k,"rNIE.death"] <- mean(p.death.A1.G1) - mean(p.death.A1.G0)

  est[k,"rNDE.qol"] <- mean(m.qol.A1.G0) - mean(m.qol.A0.G0)
  est[k,"rNIE.qol"] <- mean(m.qol.A1.G1) - mean(m.qol.A1.G0)
}

# take empirical mean of final predicted outcomes
rNDE.death <- mean(est[,"rNDE.death"])
rNDE.death
# [1] 0.07118987
rNIE.death <- mean(est[,"rNIE.death"])
rNIE.death
# [1] 0.0110088

rNDE.qol <- mean(est[,"rNDE.qol"])
rNDE.qol
# [1] -6.649923
rNIE.qol <-  mean(est[,"rNIE.qol"])
rNIE.qol
# [1] -1.585373

# 95% confidence intervals can be calculated by repeating this algorithm in 500
# bootstrap samples of the original data set.

# ---------------------------------------------------------------------------- #
# III.2) G-computation by iterative conditional expectation---------------------
# ---------------------------------------------------------------------------- #
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")

# We apply the g-computation algorithm used in the stremr package
# http://github.com/osofr/stremr

## The following steps are applied:

## 1) Fit a parametric model for the mediator conditional on A and L(0)
##    and generate predicted values by evaluating the regression setting the exposure
##    value to A=0 or A=1
### 1a) Fit parametric models for the mediator M, conditional on the exposure A and
###    baseline confounder Pr(M=1|A,L(0)) (but not conditional on L(1))
G.model <- glm(M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace,
               family = "binomial", data = df2_int)

### 1b) generate predicted probabilites by evaluating the regression setting the
###     exposure value to A=0 or to A=1
# create datasets corresponding to the counterfactual scenarios setting A=0 and A=1
data.Ais0  <- data.Ais1 <- df2_int
data.Ais0$A0_ace <- 0
data.Ais1$A0_ace <- 1

# estimate G_{A=0|L(0)} = Pr(M=1|A=0,L(0)) and G_{A=1|L(0)} = Pr(M=1|A=1,L(0))
G.Ais0.L0 <-predict(G.model, newdata = data.Ais0, type="response")
G.Ais1.L0 <-predict(G.model, newdata = data.Ais1, type="response")


## 2) Fit parametric models for the observed data for the outcome Y given the past
##    and generate predicted values by evaluating the regression setting the mediator
##    value to M=0 or to M=1
##    then calculate a weighted sum of the predicted Q.L2, with weights given by G
### 2a) fit parametric models of the outcomes given the past
Y.death.model <- glm(Y_death ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 +
                       M_smoking + A0_ace:M_smoking,
                     family = "binomial", data = df2_int)
Y.qol.model <- glm(Y_qol ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 +
                     M_smoking + A0_ace:M_smoking,
                   family = "gaussian", data = df2_int)

### 2b) generate predicted values by evaluating the regression setting the mediator
###     value to M=0 or to M=1
data.Mis0  <- data.Mis1 <- df2_int
data.Mis0$M_smoking <- 0
data.Mis1$M_smoking <- 1

Q.L2.death.Mis0 <- predict(Y.death.model, newdata = data.Mis0, type="response")
Q.L2.death.Mis1 <- predict(Y.death.model, newdata = data.Mis1, type="response")

Q.L2.qol.Mis0 <- predict(Y.qol.model, newdata = data.Mis0, type="response")
Q.L2.qol.Mis1 <- predict(Y.qol.model, newdata = data.Mis1, type="response")

### 2c) calculate a weighted sum of the predicted Q.L2, with weights given by the
###     predicted probabilities of the mediator G_{A=0|L(0)} or G_{A=1|L(0)}
# calculate barQ.L2_{A=0,G_{A=0|L(0)}}
Q.L2.death.A0.G0 <- Q.L2.death.Mis1 * G.Ais0.L0 + Q.L2.death.Mis0 * (1 - G.Ais0.L0)
Q.L2.qol.A0.G0 <- Q.L2.qol.Mis1 * G.Ais0.L0 + Q.L2.qol.Mis0 * (1 - G.Ais0.L0)

# calculate barQ.L2_{A=1,G_{A=0|L(0)}}
# note at this step, quantities are similar to barQ.L2_{A=0,G_{A=0|L(0)}}
Q.L2.death.A1.G0 <- Q.L2.death.Mis1 * G.Ais0.L0 + Q.L2.death.Mis0 * (1 - G.Ais0.L0)
Q.L2.qol.A1.G0 <- Q.L2.qol.Mis1 * G.Ais0.L0 + Q.L2.qol.Mis0 * (1 - G.Ais0.L0)

# calculate barQ.L2_{A=1,G_{A=1|L(0)}}
Q.L2.death.A1.G1 <- Q.L2.death.Mis1 * G.Ais1.L0 + Q.L2.death.Mis0 * (1 - G.Ais1.L0)
Q.L2.qol.A1.G1 <- Q.L2.qol.Mis1 * G.Ais1.L0 + Q.L2.qol.Mis0 * (1 - G.Ais1.L0)


## 3) Fit parametric models for the predicted values barQ.L2 conditional on the
##    exposure A and baseline confounders L(0)
##    and generate predicted values by evaluating the regression setting the exposure
##    value to A=0 or to A=1
### 3a) Fit parametric models for the predicted values barQ.L2 conditional on the
###    exposure A and baseline confounders L(0)
L1.death.A0.G0.model <- glm(Q.L2.death.A0.G0 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                            family = "quasibinomial", data = df2_int)
L1.death.A1.G0.model <- glm(Q.L2.death.A1.G0 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                            family = "quasibinomial", data = df2_int)
L1.death.A1.G1.model <- glm(Q.L2.death.A1.G1 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                            family = "quasibinomial", data = df2_int)

L1.qol.A0.G0.model <- glm(Q.L2.qol.A0.G0 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                          family = "gaussian", data = df2_int)
L1.qol.A1.G0.model <- glm(Q.L2.qol.A1.G0 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                          family = "gaussian", data = df2_int)
L1.qol.A1.G1.model <- glm(Q.L2.qol.A1.G1 ~ L0_male + L0_parent_low_educ_lv + A0_ace,
                          family = "gaussian", data = df2_int)

### 3b) generate predicted values by evaluating the regression setting the exposure
###     value to A=0 or to A=1
Q.L1.death.A0.G0 <- predict(L1.death.A0.G0.model, newdata = data.Ais0, type="response")
Q.L1.death.A1.G0 <- predict(L1.death.A1.G0.model, newdata = data.Ais1, type="response")
Q.L1.death.A1.G1 <- predict(L1.death.A1.G1.model, newdata = data.Ais1, type="response")

Q.L1.qol.A0.G0 <- predict(L1.qol.A0.G0.model, newdata = data.Ais0, type="response")
Q.L1.qol.A1.G0 <- predict(L1.qol.A1.G0.model, newdata = data.Ais1, type="response")
Q.L1.qol.A1.G1 <- predict(L1.qol.A1.G1.model, newdata = data.Ais1, type="response")

## 4) Estimate the marginal randomized natural direct and indirect effects
### MRDE = E(Y_{A=1,G_{A=0|L(0)}}) - E(Y_{A=0,G_{A=0|L(0)}})
### MRIE = E(Y_{A=1,G_{A=1|L(0)}}) - E(Y_{A=1,G_{A=0|L(0)}})

### for deaths:
MRDE.death <- mean(Q.L1.death.A1.G0) - mean(Q.L1.death.A0.G0)
MRDE.death
# [1] 0.0714693
MRIE.death <- mean(Q.L1.death.A1.G1) - mean(Q.L1.death.A1.G0)
MRIE.death
# [1] 0.01130057

### for quality of life
MRDE.qol <- mean(Q.L1.qol.A1.G0) - mean(Q.L1.qol.A0.G0)
MRDE.qol
# [1] -6.719193
MRIE.qol <-  mean(Q.L1.qol.A1.G1) - mean(Q.L1.qol.A1.G0)
MRIE.qol
# [1] -1.624645
