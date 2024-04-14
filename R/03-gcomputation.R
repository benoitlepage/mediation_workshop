
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
## III.1) parametric g-computation ---------------------------------------------
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
## III.2) G-computation by iterative conditional expectation--------------------
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



# ---------------------------------------------------------------------------- #
# IV) CMAverse -----------
# ---------------------------------------------------------------------------- #
# The DAG for this scientific setting is:
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
                 ereg = "logistic",
                 yreg = "logistic",
                 mreg = list("logistic"),
                 wmnomreg = list("logistic"),
                 wmdenomreg = list("logistic"),
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


## The g-formula Approach
res_gformula <- cmest(data = data.frame(df2_int[,c("L0_male","L0_parent_low_educ_lv","A0_ace")],
                                        L1=as.factor(df2_int$L1),
                                        df2_int[,c("M_smoking","Y_death")]),
                      model = "gformula",
                      outcome = "Y_death",
                      exposure = "A0_ace",
                      mediator = "M_smoking",
                      basec = c("L0_male", "L0_parent_low_educ_lv"),
                      postc = "L1",
                      EMint = TRUE,
                      mreg = list("logistic"),
                      yreg = "logistic",#  "linear", to get risk differences
                      postcreg = list("logistic"),
                      astar = 0,
                      a = 1,
                      mval = list(1), # mediator value at which the variable is controlled
                      estimation = "imputation",
                      inference = "bootstrap",
                      yval = 1, #the value of the outcome at which causal effects on the risk/odds ratio scale are estimated (used when the outcome is categorical)
                      nboot = 2)
summary(res_gformula)
# Causal Mediation Analysis
#
# # Outcome regression:
# Call:
#   glm(formula = Y_death ~ A0_ace + M_smoking + A0_ace * M_smoking +
#         L0_male + L0_parent_low_educ_lv + L1, family = binomial(),
#       data = getCall(x$reg.output$yreg)$data, weights = getCall(x$reg.output$yreg)$weights)
#
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -1.0859  -0.7097  -0.6033  -0.4944   2.0797
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -2.04033    0.05855 -34.849  < 2e-16 ***
#   A0_ace                 0.28922    0.10123   2.857  0.00428 **
#   M_smoking              0.44406    0.05569   7.974 1.53e-15 ***
#   L0_male                0.26913    0.04996   5.387 7.15e-08 ***
#   L0_parent_low_educ_lv  0.34603    0.05432   6.370 1.89e-10 ***
#   L1                     0.42894    0.05195   8.257  < 2e-16 ***
#   A0_ace:M_smoking       0.04387    0.14311   0.307  0.75919
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
# Null deviance: 10395  on 9999  degrees of freedom
# Residual deviance: 10079  on 9993  degrees of freedom
# AIC: 10093
#
# Number of Fisher Scoring iterations: 4
#
#
# # Mediator regressions:
#
# Call:
#   glm(formula = M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv +
#         L1, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
#
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -1.3340  -0.8581  -0.7529   1.2616   1.7835
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -1.36249    0.04783 -28.488  < 2e-16 ***
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
# # Regressions for mediator-outcome confounders affected by the exposure:
#
# Call:
#   glm(formula = L1 ~ A0_ace + L0_male + L0_parent_low_educ_lv,
#       family = gaussian(), data = getCall(x$reg.output$postcreg[[1L]])$data,
#       weights = getCall(x$reg.output$postcreg[[1L]])$weights)
#
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -0.5919  -0.3224  -0.2966   0.6342   0.7469
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
#   (Intercept)            0.296590   0.009214  32.190  < 2e-16 ***
#   A0_ace                 0.226072   0.014906  15.167  < 2e-16 ***
#   L0_male               -0.043457   0.009381  -4.632 3.66e-06 ***
#   L0_parent_low_educ_lv  0.069221   0.009814   7.053 1.86e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for gaussian family taken to be 0.2190925)
#
# Null deviance: 2259.1  on 9999  degrees of freedom
# Residual deviance: 2190.0  on 9996  degrees of freedom
# AIC: 13202
#
# Number of Fisher Scoring iterations: 2
#
#
# # Effect decomposition on the odds ratio scale via the g-formula approach
#
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
#                   Estimate Std.error  95% CIL 95% CIU  P.val
#   Rcde            1.525902  0.060467 1.398987   1.480 <2e-16 ***
#   rRpnde          1.481091  0.039678 1.478889   1.532 <2e-16 ***
#   rRtnde          1.487634  0.040141 1.465710   1.520 <2e-16 ***
#   rRpnie          1.050494  0.018867 1.037797   1.063 <2e-16 ***
#   rRtnie          1.055135  0.019266 1.028549   1.054 <2e-16 ***
#   Rte             1.562751  0.070331 1.521109   1.616 <2e-16 ***
#   ERcde           0.432152  0.054697 0.330836   0.404 <2e-16 ***
#   rERintref       0.048939  0.015019 0.127900   0.148 <2e-16 ***
#   rERintmed       0.031165  0.011786 0.004466   0.020 <2e-16 ***
#   rERpnie         0.050494  0.018867 0.037805   0.063 <2e-16 ***
#   ERcde(prop)     0.767928  0.016340 0.634669   0.657 <2e-16 ***
#   rERintref(prop) 0.086965  0.056893 0.208055   0.284 <2e-16 ***
#   rERintmed(prop) 0.055380  0.018177 0.008450   0.033 <2e-16 ***
#   rERpnie(prop)   0.089727  0.022376 0.072390   0.102 <2e-16 ***
#   rpm             0.145108  0.040553 0.080839   0.135 <2e-16 ***
#   rint            0.142345  0.038716 0.240926   0.293 <2e-16 ***
#   rpe             0.232072  0.016340 0.343378   0.365 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Rcde: controlled direct effect odds ratio;
# rRpnde: randomized analogue of pure natural direct effect odds ratio;
# rRtnde: randomized analogue of total natural direct effect odds ratio;
# rRpnie: randomized analogue of pure natural indirect effect odds ratio;
# rRtnie: randomized analogue of total natural indirect effect odds ratio;
# Rte: total effect odds ratio; ERcde: excess relative risk due to controlled direct effect;
# rERintref: randomized analogue of excess relative risk due to reference interaction;
# rERintmed: randomized analogue of excess relative risk due to mediated interaction;
# rERpnie: randomized analogue of excess relative risk due to pure natural indirect effect;
# ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref;
# rERintmed(prop): proportion rERintmed;
# rERpnie(prop): proportion rERpnie;
# rpm: randomized analogue of overall proportion mediated;
# rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)
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

