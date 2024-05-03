
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
# devtools::install_github("BS1125/CMAverse")
library(CMAverse)
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
cmdag(outcome = "Y_death", exposure = "A0_ace", mediator = "M_smoking",
      basec = c("L0_male", "L0_parent_low_educ_lv"), postc = "L1", node = TRUE, text_col = "white")
# The CMAverse package can be used to estimate CDE, rNDE and rNIE by parametric g-computation
# (using the 'gformula' option in the )

### Continuous outcome
## The g-formula Approach
set.seed(1234)
res_gformula_Qol <- cmest(data = df2_int, #data.frame(df2_int[,c("L0_male",
                                                       #"L0_parent_low_educ_lv",
                                                       #"A0_ace")],
                                           #L1=as.factor(df2_int$L1),
                                           #df2_int[,c("M_smoking","Y_qol")]),
                         model = "gformula", # for parametric g-computation
                         outcome = "Y_qol", # outcome variable
                         exposure = "A0_ace", # exposure variable
                         mediator = "M_smoking", # mediator
                         basec = c("L0_male",
                                   "L0_parent_low_educ_lv"), # confounders
                         postc = "L1", # intermediate confounder (post-exposure)
                         EMint = TRUE, # exposures*mediator interaction
                         mreg = list("logistic"), # g(M=1|L1,A,L0)
                         yreg = "linear",# Qbar.L2 = P(Y=1|M,L1,A,L0)
                         postcreg = list("logistic"), # Qbar.L1 = P(L1=1|A,L0)
                         astar = 0,
                         a = 1,
                         mval = list(0), # do(M=0) to estimate CDE_m
                         estimation = "imputation", # if model= gformula
                         inference = "bootstrap",
                         boot.ci.type = "per", # forpercentile, other option: "bca"
                         nboot = 2) # we should use a large number of bootstrap samples
summary(res_gformula_Qol)

### 1) Estimation of Qbar.Y = P(Y=1|M,L1,A,L0) with A*M interaction,
### Outcome regression:
# Call:
#   glm(formula = Y_qol ~ A0_ace + M_smoking + A0_ace * M_smoking +
#         L0_male + L0_parent_low_educ_lv + L1, family = gaussian(),
#       data = getCall(x$reg.output$yreg)$data, weights = getCall(x$reg.output$yreg)$weights)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)            74.8247     0.2133 350.823  < 2e-16 ***
#   A0_ace                 -3.7014     0.4295  -8.617  < 2e-16 ***
#   M_smoking              -8.6336     0.2331 -37.042  < 2e-16 ***
#   L0_male                -0.7280     0.2019  -3.605 0.000313 ***
#   L0_parent_low_educ_lv  -2.8828     0.2116 -13.621  < 2e-16 ***
#   L1                     -5.1668     0.2189 -23.608  < 2e-16 ***
#   A0_ace:M_smoking       -5.5119     0.6440  -8.559  < 2e-16 ***

### 2) Estimation of g(M=1|L1,A,L0), model of the mediator
### Mediator regressions:
# Call:
#   glm(formula = M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv +
#         L1, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -1.36249    0.04783 -28.488  < 2e-16 ***
#   A0_ace                 0.30994    0.06668   4.648 3.35e-06 ***
#   L0_male                0.24661    0.04369   5.644 1.66e-08 ***
#   L0_parent_low_educ_lv  0.30628    0.04650   6.587 4.50e-11 ***
#   L1                     0.86045    0.04493  19.152  < 2e-16 ***

### 3) Estimation of Qbar.L1 = P(L1=1|A,L0), model of intermediate confounder
### Regressions for mediator-outcome confounders affected by the exposure:
# Call:
#   glm(formula = L1 ~ A0_ace + L0_male + L0_parent_low_educ_lv,
#       family = binomial(), data = getCall(x$reg.output$postcreg[[1L]])$data,
#       weights = getCall(x$reg.output$postcreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -0.86983    0.04292 -20.267  < 2e-16 ***
#   A0_ace                 0.94354    0.06475  14.572  < 2e-16 ***
#   L0_male               -0.19827    0.04289  -4.622 3.80e-06 ***
#   L0_parent_low_educ_lv  0.32047    0.04556   7.034 2.01e-12 ***

### 4) Effect decomposition on the mean difference scale via the g-formula approach
#
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
# Estimate Std.error   95% CIL 95% CIU  P.val
#   cde           -5.863750  0.233488 -4.933234  -4.620 <2e-16 ***
#   rpnde         -7.565835  0.199867 -6.689581  -6.421 <2e-16 ***
#   rtnde         -8.463729  0.157383 -7.266085  -7.055 <2e-16 ***
#   rpnie         -1.406410  0.021876 -0.971101  -0.942 <2e-16 ***
#   rtnie         -2.304304  0.064359 -1.604682  -1.518 <2e-16 ***
#   te            -9.870139  0.135507 -8.207796  -8.026 <2e-16 ***
#   rintref       -1.702085  0.033622 -1.801518  -1.756 <2e-16 ***
#   rintmed       -0.897894  0.042484 -0.633581  -0.577 <2e-16 ***
#   cde(prop)      0.594090  0.018945  0.575575   0.601 <2e-16 ***
#   rintref(prop)  0.172448  0.007802  0.213991   0.224 <2e-16 ***
#   rintmed(prop)  0.090971  0.006479  0.070244   0.079 <2e-16 ***
#   rpnie(prop)    0.142491  0.004663  0.114737   0.121 <2e-16 ***
#   rpm            0.233462  0.011142  0.184981   0.200 <2e-16 ***
#   rint           0.263419  0.014282  0.284235   0.303 <2e-16 ***
#   rpe            0.405910  0.018945  0.398973   0.424 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# cde: controlled direct effect;
# rpnde: randomized analogue of pure natural direct effect;
# rtnde: randomized analogue of total natural direct effect;
# rpnie: randomized analogue of pure natural indirect effect;
# rtnie: randomized analogue of total natural indirect effect;
# te: total effect; rintref: randomized analogue of reference interaction;
# rintmed: randomized analogue of mediated interaction;
# cde(prop): proportion cde;
# rintref(prop): proportion rintref;
# rintmed(prop): proportion rintmed;
# rpnie(prop): proportion rpnie;
# rpm: randomized analogue of overall proportion mediated;
# rint: randomized analogue of overall proportion attributable to interaction;
# rpe: randomized analogue of overall proportion eliminated

## The g-formula Approach
set.seed(1234)
res_gformula_OR <- cmest(data = data.frame(df2_int[,c("L0_male","L0_parent_low_educ_lv","A0_ace")],
                                        L1=as.factor(df2_int$L1),
                                        df2_int[,c("M_smoking","Y_death")]),
                      model = "gformula",
                      outcome = "Y_death",
                      exposure = "A0_ace",
                      mediator = "M_smoking",
                      basec = c("L0_male", "L0_parent_low_educ_lv"),
                      postc = "L1",
                      EMint = TRUE,
                      mreg = list("logistic"), # g(M=1|L1,A,L0)
                      yreg = "logistic",# Qbar.L2 = P(Y=1|M,L1,A,L0)
                      postcreg = list("logistic"), # Qbar.L1 = P(L1=1|A,L0)
                      astar = 0,
                      a = 1,
                      mval = list(0), # do(M=0) to estimate CDE_m
                      estimation = "imputation", # parametric g-comp if model= gformula
                      inference = "bootstrap",
                      boot.ci.type = "per", # forpercentile, other option: "bca"
                      # yval = 1, #the value of the outcome at which causal effects on the risk/odds ratio scale are estimated (used when the outcome is categorical)
                      nboot = 2) # we should use a large number of bootstrap samples
summary(res_gformula_OR)
### Causal Mediation Analysis

## 1) Estimation of Qbar.L2 = P(Y=1|M,L1,A,L0) with A*M interaction, model of the outcome
# # Outcome regression:
#
# Call:
#   glm(formula = Y_death ~ A0_ace + M_smoking + A0_ace * M_smoking +
#         L0_male + L0_parent_low_educ_lv + L1, family = binomial(),
#       data = getCall(x$reg.output$yreg)$data, weights = getCall(x$reg.output$yreg)$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -2.04033    0.05855 -34.849  < 2e-16 ***
#   A0_ace                 0.28922    0.10123   2.857  0.00428 **
#   M_smoking              0.44406    0.05569   7.974 1.53e-15 ***
#   L0_male                0.26913    0.04996   5.387 7.15e-08 ***
#   L0_parent_low_educ_lv  0.34603    0.05432   6.370 1.89e-10 ***
#   L11                    0.42894    0.05195   8.257  < 2e-16 ***
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

## 2) Estimation of g(M=1|L1,A,L0), model of the mediator
# # Mediator regressions:
#
# Call:
#   glm(formula = M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv +
#         L1, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -1.36249    0.04783 -28.488  < 2e-16 ***
#   A0_ace                 0.30994    0.06668   4.648 3.35e-06 ***
#   L0_male                0.24661    0.04369   5.644 1.66e-08 ***
#   L0_parent_low_educ_lv  0.30628    0.04650   6.587 4.50e-11 ***
#   L11                    0.86045    0.04493  19.152  < 2e-16 ***
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

## 3) Estimation of Qbar.L1 = P(L1=1|A,L0), model of intermediate confounder
# # Regressions for mediator-outcome confounders affected by the exposure:
#
# Call:
#   glm(formula = L1 ~ A0_ace + L0_male + L0_parent_low_educ_lv,
#       family = binomial(), data = getCall(x$reg.output$postcreg[[1L]])$data,
#       weights = getCall(x$reg.output$postcreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -0.86983    0.04292 -20.267  < 2e-16 ***
#   A0_ace                 0.94354    0.06475  14.572  < 2e-16 ***
#   L0_male               -0.19827    0.04289  -4.622 3.80e-06 ***
#   L0_parent_low_educ_lv  0.32047    0.04556   7.034 2.01e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
# Null deviance: 12883  on 9999  degrees of freedom
# Residual deviance: 12587  on 9996  degrees of freedom
# AIC: 12595
#
# Number of Fisher Scoring iterations: 4
#


## 4) Effect decomposition on the odds ratio scale via the g-formula approach
#
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
# Estimate Std.error  95% CIL 95% CIU  P.val
#   Rcde            1.592115  0.121242 1.401748   1.565 <2e-16 ***
#   rRpnde          1.610141  0.064880 1.458433   1.546 <2e-16 ***
#   rRtnde          1.621371  0.046720 1.479227   1.542 <2e-16 ***
#   rRpnie          1.076471  0.001714 1.036683   1.039 <2e-16 ***
#   rRtnie          1.083979  0.014539 1.034267   1.054 <2e-16 ***
#   Rte             1.745359  0.045898 1.536897   1.599 <2e-16 ***
#   ERcde           0.392219  0.073776 0.277167   0.376 <2e-16 ***
#   rERintref       0.217921  0.008895 0.169383   0.181 <2e-16 ***
#   rERintmed       0.058747  0.017268 0.016243   0.039 <2e-16 ***
#   rERpnie         0.076471  0.001714 0.036684   0.039 <2e-16 ***
#   ERcde(prop)     0.526216  0.083694 0.515858   0.628 <2e-16 ***
#   rERintref(prop) 0.292371  0.040769 0.283120   0.338 <2e-16 ***
#   rERintmed(prop) 0.078817  0.034491 0.027264   0.074 <2e-16 ***
#   rERpnie(prop)   0.102596  0.008434 0.061314   0.073 <2e-16 ***
#   rpm             0.181413  0.042925 0.088578   0.146 <2e-16 ***
#   rint            0.371188  0.075260 0.310384   0.411 <2e-16 ***
#   rpe             0.473784  0.083694 0.371698   0.484 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Rcde: controlled direct effect odds ratio;
# rRpnde: randomized analogue of pure natural direct effect odds ratio;
# rRtnde: randomized analogue of total natural direct effect odds ratio;
# rRpnie: randomized analogue of pure natural indirect effect odds ratio;
# rRtnie: randomized analogue of total natural indirect effect odds ratio;
# Rte: total effect odds ratio;
# ERcde: excess relative risk due to controlled direct effect;
# rERintref: randomized analogue of excess relative risk due to reference interaction;
# rERintmed: randomized analogue of excess relative risk due to mediated interaction;
# rERpnie: randomized analogue of excess relative risk due to pure natural indirect effect;
# ERcde(prop): proportion ERcde;
# rERintref(prop): proportion rERintref;
# rERintmed(prop): proportion rERintmed;
# rERpnie(prop): proportion rERpnie;
# rpm: randomized analogue of overall proportion mediated;
# rint: randomized analogue of overall proportion attributable to interaction;
# rpe: randomized analogue of overall proportion eliminated)


## How to get risk differences ? apply a gaussian glm for the Qbar.L2 model
# setting do(M=0) for the controlled direct effect
set.seed(12345)
res_gformula_RD_M0 <- cmest(data = data.frame(df2_int[,c("L0_male","L0_parent_low_educ_lv","A0_ace")],
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
                            yreg = "linear",#  instead of "logistic" to get risk differences
                            postcreg = list("logistic"),
                            astar = 0,
                            a = 1,
                            mval = list(0), # mediator value at which the variable is controlled
                            estimation = "imputation",
                            inference = "bootstrap",
                            boot.ci.type = "per", # forpercentile, other option: "bca"
                            #yval = 1, # not necessary here
                            nboot = 2) # we should use a large number of bootstrap samples
summary(res_gformula_RD_M0)
## 1) Outcome regression: Qbar.L2 = P(Y=1|M,L1,A,L0) with A*M interaction
# Call:
#   glm(formula = Y_death ~ A0_ace + M_smoking + A0_ace * M_smoking +
#         L0_male + L0_parent_low_educ_lv + L1, family = gaussian(),
#       data = getCall(x$reg.output$yreg)$data, weights = getCall(x$reg.output$yreg)$weights)
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
#   (Intercept)           0.100674   0.008574  11.742  < 2e-16 ***
#   A0_ace                0.046456   0.017268   2.690  0.00715 **
#   M_smoking             0.074378   0.009370   7.938 2.27e-15 ***
#   L0_male               0.043661   0.008117   5.379 7.67e-08 ***
#   L0_parent_low_educ_lv 0.053711   0.008508   6.313 2.86e-10 ***
#   L11                   0.073753   0.008798   8.383  < 2e-16 ***
#   A0_ace:M_smoking      0.030159   0.025887   1.165  0.24404

## 2) Mediator regressions:  g(M=1|L1,A,L0)
# Call:
#   glm(formula = M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv +
#         L1, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -1.36249    0.04783 -28.488  < 2e-16 ***
#   A0_ace                 0.30994    0.06668   4.648 3.35e-06 ***
#   L0_male                0.24661    0.04369   5.644 1.66e-08 ***
#   L0_parent_low_educ_lv  0.30628    0.04650   6.587 4.50e-11 ***
#   L11                    0.86045    0.04493  19.152  < 2e-16 ***

## 3) Regressions for mediator-outcome confounders affected by the exposure: Qbar.L1 = P(L1=1|A,L0)
# Call:
#   glm(formula = L1 ~ A0_ace + L0_male + L0_parent_low_educ_lv,
#       family = binomial(), data = getCall(x$reg.output$postcreg[[1L]])$data,
#       weights = getCall(x$reg.output$postcreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -0.86983    0.04292 -20.267  < 2e-16 ***
#   A0_ace                 0.94354    0.06475  14.572  < 2e-16 ***
#   L0_male               -0.19827    0.04289  -4.622 3.80e-06 ***
#   L0_parent_low_educ_lv  0.32047    0.04556   7.034 2.01e-12 ***

## 4) Effect decomposition on the mean difference scale via the g-formula approach
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
#                  Estimate Std.error   95% CIL 95% CIU  P.val
#   cde           0.0638914 0.0158641 0.0467659   0.068 <2e-16 *** CDE=0.06289087 by previous parametric g-comp
#   rpnde         0.0734700 0.0307357 0.0598528   0.101 <2e-16 *** rPNDE = 0.07118987 by previous parametric g-comp
#   rtnde         0.0771555 0.0333264 0.0647151   0.109 <2e-16 ***
#   rpnie         0.0090890 0.0033498 0.0050769   0.010 <2e-16 ***
#   rtnie         0.0127744 0.0007591 0.0134198   0.014 <2e-16 *** rTNIE = 0.0110088 by previous parametric g-comp
#   te            0.0862445 0.0299766 0.0742924   0.115 <2e-16 ***
#   rintref       0.0095786 0.0148716 0.0130869   0.033 <2e-16 ***
#   rintmed       0.0036855 0.0025907 0.0048623   0.008 <2e-16 ***
#   cde(prop)     0.7408176 0.0263715 0.5945733   0.630 <2e-16 ***
#   rintref(prop) 0.1110634 0.0841501 0.1744983   0.288 <2e-16 ***
#   rintmed(prop) 0.0427328 0.0055167 0.0653396   0.073 <2e-16 ***
#   rpnie(prop)   0.1053862 0.0632953 0.0451212   0.130 <2e-16 ***
#   rpm           0.1481190 0.0577786 0.1178726   0.195 <2e-16 ***
#   rint          0.1537963 0.0896668 0.2398378   0.360 <2e-16 ***
#   rpe           0.2591824 0.0263715 0.3699965   0.405 <2e-16 ***
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

### How to find analogues for the 3-way and 4-way decomposition?
## TE
res_gformula_RD_M0$effect.pe["te"]
## rPNDE
res_gformula_RD_M0$effect.pe["rpnde"]
## rPNIE
res_gformula_RD_M0$effect.pe["rpnie"]
## CDE
res_gformula_RD_M0$effect.pe["cde"]


## MI= TE - rPNDE - rPNIE
res_gformula_RD_M0$effect.pe["te"] - res_gformula_RD_M0$effect.pe["rpnde"] - res_gformula_RD_M0$effect.pe["rpnie"]
# 0.003685471
res_gformula_RD_M0$effect.pe["rintmed"]
# 0.003685471 # This is how the mediated interaction is calculated

## RE = PNDE - CDE_{M=0}
res_gformula_RD_M0$effect.pe["rpnde"] - res_gformula_RD_M0$effect.pe["cde"]
# 0.009578606
res_gformula_RD_M0$effect.pe["rintref"]
# 0.009578606 # This is how the reference interaction is calculated


### setting do(M=1) for the controlled direct effect
set.seed(12345)
res_gformula_RD_M1 <- cmest(data = data.frame(df2_int[,c("L0_male","L0_parent_low_educ_lv","A0_ace")],
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
                            yreg = "linear",#  instead of "logistic" to get risk differences
                            postcreg = list("logistic"),
                            astar = 0,
                            a = 1,
                            mval = list(1), # mediator value at which the variable is controlled
                            estimation = "imputation",
                            inference = "bootstrap",
                            #yval = 1, #the value of the outcome at which causal effects on the risk/odds ratio scale are estimated (used when the outcome is categorical)
                            nboot = 2)
summary(res_gformula_RD_M1)
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
# Estimate  Std.error    95% CIL 95% CIU  P.val
#   cde            0.0940507  0.0623648  0.0863933   0.170 <2e-16 *** CDE_{M=1}=0.08751016 by previous parametric g-comp
#   rpnde          0.0734700  0.0307357  0.0598528   0.101 <2e-16 *** rPNDE = 0.07118987 by previous parametric g-comp
#   rtnde          0.0771555  0.0333264  0.0647151   0.109 <2e-16 ***
#   rpnie          0.0090890  0.0033498  0.0050769   0.010 <2e-16 ***
#   rtnie          0.0127744  0.0007591  0.0134198   0.014 <2e-16 *** rTNIE = 0.0110088 by previous parametric g-comp
#   te             0.0862445  0.0299766  0.0742924   0.115 <2e-16 ***
#   rintref       -0.0205807  0.0316291 -0.0690343  -0.027 <2e-16 ***
#   rintmed        0.0036855  0.0025907  0.0048623   0.008 <2e-16 ***
#   cde(prop)      1.0905135  0.2413267  1.1581348   1.482 <2e-16 ***
#   rintref(prop) -0.2386325  0.1835481 -0.6002305  -0.354 <2e-16 ***
#   rintmed(prop)  0.0427328  0.0055167  0.0653396   0.073 <2e-16 ***
#   rpnie(prop)    0.1053862  0.0632953  0.0451212   0.130 <2e-16 ***
#   rpm            0.1481190  0.0577786  0.1178726   0.195 <2e-16 ***
#   rint          -0.1958997  0.1780314 -0.5274792  -0.288 <2e-16 ***
#   rpe           -0.0905135  0.2413267 -0.4823579  -0.158 <2e-16 ***


### For quantitative outcomes
set.seed(12345)
res_gformula_QoL_M0 <- cmest(data = data.frame(df2_int[,c("L0_male","L0_parent_low_educ_lv","A0_ace")],
                                              L1=as.factor(df2_int$L1),
                                              df2_int[,c("M_smoking","Y_qol")]),
                            model = "gformula",
                            outcome = "Y_qol",
                            exposure = "A0_ace",
                            mediator = "M_smoking",
                            basec = c("L0_male", "L0_parent_low_educ_lv"),
                            postc = "L1",
                            EMint = TRUE,
                            mreg = list("logistic"),
                            yreg = "linear",#
                            postcreg = list("logistic"),
                            astar = 0,
                            a = 1,
                            mval = list(0), # mediator value at which the variable is controlled
                            estimation = "imputation",
                            inference = "bootstrap",
                            boot.ci.type = "per", # forpercentile, other option: "bca"
                            #yval = 1, # not necessary here
                            nboot = 2) # we should use a large number of bootstrap samples
summary(res_gformula_QoL_M0)
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
# Estimate Std.error   95% CIL 95% CIU  P.val
#   cde           -4.922873  0.284894 -5.116411  -4.734 <2e-16 *** CDE_{M=0} = -4.838654 by previous parametric g-comp
#   rpnde         -6.673463  0.394323 -7.086829  -6.557 <2e-16 *** rPNDE = -6.649923 by previous parametric g-comp
#   rtnde         -7.347022  0.601495 -7.829339  -7.021 <2e-16 ***
#   rpnie         -1.055023  0.242564 -1.054962  -0.729 <2e-16 ***
#   rtnie         -1.728582  0.449736 -1.797472  -1.193 <2e-16 *** rTNIE = -1.585373 by previous parametric g-comp
#   te            -8.402045  0.844058 -8.884301  -7.750 <2e-16 ***
#   rintref       -1.750590  0.109429 -1.970417  -1.823 <2e-16 ***
#   rintmed       -0.673558  0.207172 -0.742510  -0.464 <2e-16 ***
#   cde(prop)      0.585914  0.025973  0.576014   0.611 <2e-16 ***
#   rintref(prop)  0.208353  0.010040  0.221833   0.235 <2e-16 ***
#   rintmed(prop)  0.080166  0.017638  0.059797   0.083 <2e-16 ***
#   rpnie(prop)    0.125567  0.018375  0.093973   0.119 <2e-16 ***
#   rpm            0.205733  0.036012  0.153771   0.202 <2e-16 ***
#   rint           0.288519  0.007598  0.295119   0.305 <2e-16 ***
#   rpe            0.414086  0.025973  0.389092   0.424 <2e-16 ***

set.seed(12345)
res_gformula_QoL_M1 <- cmest(data = data.frame(df2_int[,c("L0_male","L0_parent_low_educ_lv","A0_ace")],
                                               L1=as.factor(df2_int$L1),
                                               df2_int[,c("M_smoking","Y_qol")]),
                             model = "gformula",
                             outcome = "Y_qol",
                             exposure = "A0_ace",
                             mediator = "M_smoking",
                             basec = c("L0_male", "L0_parent_low_educ_lv"),
                             postc = "L1",
                             EMint = TRUE,
                             mreg = list("logistic"),
                             yreg = "linear",#
                             postcreg = list("logistic"),
                             astar = 0,
                             a = 1,
                             mval = list(1), # mediator value at which the variable is controlled
                             estimation = "imputation",
                             inference = "bootstrap",
                             boot.ci.type = "per", # forpercentile, other option: "bca"
                             #yval = 1, # not necessary here
                             nboot = 2) # we should use a large number of bootstrap samples
summary(res_gformula_QoL_M1)
# Direct counterfactual imputation estimation with
# bootstrap standard errors, percentile confidence intervals and p-values
#
# Estimate Std.error   95% CIL 95% CIU  P.val
#   cde           -10.43481   0.53289 -11.07766 -10.362 <2e-16 *** CDE_{M=1} = -10.35059 by previous parametric g-comp
#   rpnde          -6.67346   0.39432  -7.08683  -6.557 <2e-16 *** rPNDE = -6.649923 by previous parametric g-comp
#   rtnde          -7.34702   0.60149  -7.82934  -7.021 <2e-16 ***
#   rpnie          -1.05502   0.24256  -1.05496  -0.729 <2e-16 ***
#   rtnie          -1.72858   0.44974  -1.79747  -1.193 <2e-16 *** rTNIE = -1.585373 by previous parametric g-comp
#   te             -8.40205   0.84406  -8.88430  -7.750 <2e-16 ***
#   rintref         3.76134   0.13857   3.80467   3.991 <2e-16 ***
#   rintmed        -0.67356   0.20717  -0.74251  -0.464 <2e-16 ***
#   cde(prop)       1.24194   0.06707   1.24719   1.337 <2e-16 ***
#   rintref(prop)  -0.44767   0.03106  -0.49107  -0.449 <2e-16 ***
#   rintmed(prop)   0.08017   0.01764   0.05980   0.083 <2e-16 ***
#   rpnie(prop)     0.12557   0.01837   0.09397   0.119 <2e-16 ***
#   rpm             0.20573   0.03601   0.15377   0.202 <2e-16 ***
#   rint           -0.36750   0.04870  -0.43127  -0.366 <2e-16 ***
#   rpe            -0.24194   0.06707  -0.33730  -0.247 <2e-16 ***
