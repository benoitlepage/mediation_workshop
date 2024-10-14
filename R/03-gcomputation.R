
### test programs for estimations based on G-computation

rm(list=ls())

df2 <- read.csv(file = "data/df2.csv")
df2_int <- read.csv(file = "data/df2_int.csv")


# ---------------------------------------------------------------------------- #
# I) Estimation of the Average Total Effect (ATE) ------------------------------
# ---------------------------------------------------------------------------- #

# 1. Estimate Qbar
Q.tot.death <- glm(Y_death ~ A0_PM2.5 + L0_male + L0_soc_env, family = "binomial", data = df2_int)
Q.tot.qol <- glm(Y_qol ~ A0_PM2.5 + L0_male + L0_soc_env, family = "gaussian", data = df2_int)

# 2. Predict an outcome for each subject, setting A=0 and A=1
# prepare data sets used to predict the outcome under the counterfactual scenarios setting A=0 and A=1
data.A1 <- data.A0 <- df2_int
data.A1$A0_PM2.5 <- 1
data.A0$A0_PM2.5 <- 0

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
B <- 1000
bootstrap.estimates <- data.frame(matrix(NA, nrow = B, ncol = 2))
colnames(bootstrap.estimates) <- c("boot.death.est", "boot.qol.est")
for (b in 1:B){
  # sample the indices 1 to n with replacement
  bootIndices <- sample(1:nrow(df2_int), replace=T)
  bootData <- df2_int[bootIndices,]

  if (round(b/100, 0) == b/100 ) print(paste0("bootstrap number ",b))

  Q.tot.death <- glm(Y_death ~ A0_PM2.5 + L0_male + L0_soc_env,
                     family = "binomial", data = bootData)
  Q.tot.qol <- glm(Y_qol ~ A0_PM2.5 + L0_male + L0_soc_env,
                   family = "gaussian", data = bootData)

  boot.A.1 <- boot.A.0 <- bootData
  boot.A.1$A0_PM2.5 <- 1
  boot.A.0$A0_PM2.5 <- 0

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
# [1] 0.05612907 0.10928734

IC95.ATE.qol <- c(ATE.qol.gcomp - qnorm(0.975)*sd(bootstrap.estimates[,"boot.qol.est"]),
                  ATE.qol.gcomp + qnorm(0.975)*sd(bootstrap.estimates[,"boot.qol.est"]) )
IC95.ATE.qol
# [1] -9.157856 -7.563526



# ---------------------------------------------------------------------------- #
# II) Estimation of the Control direct Effects (CDE) ---------------------------
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
## II.1) parametric g-computation ----------------------------------------------
# ---------------------------------------------------------------------------- #
rm(list=ls())
df2_int <- read.csv(file = "./data/df2_int.csv")


# a. fit parametric models to estimate the density of intermediate confounders
L1.model <- glm(L1 ~ L0_male + L0_soc_env + A0_PM2.5, family = "binomial", data = df2_int)

# b. fit parametric models for the outcome
Y.death.model <- glm(Y_death ~ L0_male + L0_soc_env + A0_PM2.5 + L1 + M_diabetes + A0_PM2.5:M_diabetes, family = "binomial", data = df2_int)
Y.qol.model <- glm(Y_qol ~ L0_male + L0_soc_env + A0_PM2.5 + L1 + M_diabetes + A0_PM2.5:M_diabetes, family = "gaussian", data = df2_int)

# c. simulate L1 values under the counterfactual scenarios setting A0=0 or A0=1
set.seed(54321)
data.A0  <- data.A1 <- df2_int
data.A0$A0_PM2.5 <- 0
data.A1$A0_PM2.5 <- 1
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
data.A0.M0$M_diabetes <- 0
data.A0.M1$M_diabetes <- 1
data.A1.M0$M_diabetes <- 0
data.A1.M1$M_diabetes <- 1

p.death.A0.M0 <- predict(Y.death.model, newdata = data.A0.M0, type = "response")
p.death.A1.M0 <- predict(Y.death.model, newdata = data.A1.M0, type = "response")
p.death.A0.M1 <- predict(Y.death.model, newdata = data.A0.M1, type = "response")
p.death.A1.M1 <- predict(Y.death.model, newdata = data.A1.M1, type = "response")

m.qol.A0.M0 <- predict(Y.qol.model, newdata = data.A0.M0, type = "response")
m.qol.A1.M0 <- predict(Y.qol.model, newdata = data.A1.M0, type = "response")
m.qol.A0.M1 <- predict(Y.qol.model, newdata = data.A0.M1, type = "response")
m.qol.A1.M1 <- predict(Y.qol.model, newdata = data.A1.M1, type = "response")

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
df2_int <- read.csv(file = "data/df2_int.csv")

# 1) Regress the outcome on L0, A, L1 and M (and the A*M interaction if appropriate)
Y.death.model <- glm(Y_death ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                       M_diabetes + A0_PM2.5:M_diabetes,
                     family = "binomial", data = df2_int)
Y.qol.model <- glm(Y_qol ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                     M_diabetes + A0_PM2.5:M_diabetes,
                   family = "gaussian", data = df2_int)

# 2) Generate predicted values by evaluating the regression setting the exposure
#    and the mediator at exposure history of interest:
#    {A=1,M=0},{A=0,M=0},{A=1,M=1},{A=0,M=1}
data.Ais0.Mis0 <- data.Ais0.Mis1 <- df2_int
data.Ais1.Mis0 <- data.Ais1.Mis1 <- df2_int

data.Ais0.Mis0$A0_PM2.5 <- 0
data.Ais0.Mis0$M_diabetes <- 0

data.Ais0.Mis1$A0_PM2.5 <- 0
data.Ais0.Mis1$M_diabetes <- 1

data.Ais1.Mis0$A0_PM2.5 <- 1
data.Ais1.Mis0$M_diabetes <- 0

data.Ais1.Mis1$A0_PM2.5 <- 1
data.Ais1.Mis1$M_diabetes <- 1


Q.L2.death.A0M0 <- predict(Y.death.model, newdata = data.Ais0.Mis0, type="response")
Q.L2.death.A0M1 <- predict(Y.death.model, newdata = data.Ais0.Mis1, type="response")
Q.L2.death.A1M0 <- predict(Y.death.model, newdata = data.Ais1.Mis0, type="response")
Q.L2.death.A1M1 <- predict(Y.death.model, newdata = data.Ais1.Mis1, type="response")

Q.L2.qol.A0M0 <- predict(Y.qol.model, newdata = data.Ais0.Mis0, type="response")
Q.L2.qol.A0M1 <- predict(Y.qol.model, newdata = data.Ais0.Mis1, type="response")
Q.L2.qol.A1M0 <- predict(Y.qol.model, newdata = data.Ais1.Mis0, type="response")
Q.L2.qol.A1M1 <- predict(Y.qol.model, newdata = data.Ais1.Mis1, type="response")


# 3) Regress the predicted values conditional on the exposure A
#    and baseline confounders L(0)
L1.death.A0M0.model <- glm(Q.L2.death.A0M0 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)
L1.death.A0M1.model <- glm(Q.L2.death.A0M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)
L1.death.A1M0.model <- glm(Q.L2.death.A1M0 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)
L1.death.A1M1.model <- glm(Q.L2.death.A1M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                           family = "quasibinomial", data = df2_int)

L1.qol.A0M0.model <- glm(Q.L2.qol.A0M0 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "gaussian", data = df2_int)
L1.qol.A0M1.model <- glm(Q.L2.qol.A0M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "gaussian", data = df2_int)
L1.qol.A1M0.model <- glm(Q.L2.qol.A1M0 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "gaussian", data = df2_int)
L1.qol.A1M1.model <- glm(Q.L2.qol.A1M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "gaussian", data = df2_int)

# 4) generate predicted values by evaluating the regression at exposure
#    of interest: {A=1} & {A=0}
Q.L1.death.A0M0 <- predict(L1.death.A0M0.model, newdata = data.Ais0.Mis0, type="response")
Q.L1.death.A0M1 <- predict(L1.death.A0M1.model, newdata = data.Ais0.Mis1, type="response")
Q.L1.death.A1M0 <- predict(L1.death.A1M0.model, newdata = data.Ais1.Mis0, type="response")
Q.L1.death.A1M1 <- predict(L1.death.A1M1.model, newdata = data.Ais1.Mis1, type="response")

Q.L1.qol.A0M0 <- predict(L1.qol.A0M0.model, newdata = data.Ais0.Mis0, type="response")
Q.L1.qol.A0M1 <- predict(L1.qol.A0M1.model, newdata = data.Ais0.Mis1, type="response")
Q.L1.qol.A1M0 <- predict(L1.qol.A1M0.model, newdata = data.Ais1.Mis0, type="response")
Q.L1.qol.A1M1 <- predict(L1.qol.A1M1.model, newdata = data.Ais1.Mis1, type="response")

# 5) Take empirical mean of final predicted outcomes to estimate CDE
# CDE setting M=0
CDE.death.m0.gcomp.ice <- mean(Q.L1.death.A1M0) - mean(Q.L1.death.A0M0)
CDE.death.m0.gcomp.ice
# [1] 0.06342833

CDE.qol.m0.gcomp.ice <- mean(Q.L1.qol.A1M0) - mean(Q.L1.qol.A0M0)
CDE.qol.m0.gcomp.ice
# [1] -4.869509

# CDE setting M=1
CDE.death.m1.gcomp.ice <- mean(Q.L1.death.A1M1) - mean(Q.L1.death.A0M1)
CDE.death.m1.gcomp.ice
# [1] 0.08812318

CDE.qol.m1.gcomp.ice <- mean(Q.L1.qol.A1M1) - mean(Q.L1.qol.A0M1)
CDE.qol.m1.gcomp.ice
# [1] -10.38144


### G-computation by iterative conditional expectation using the ltmle package
library(ltmle)
rm(list=ls())
df2_int <- read.csv(file = "./data/df2_int.csv")
df.death <- subset(df2_int, select = -Y_qol)
df.qol <- subset(df2_int, select = -Y_death)

# the data set should be composed of continuous or binary variables,
# ordered following the cause-effect sequence of each variables.
# Note that within a set of exposures or intermediate confounders measured at a
# single discrete time t, any causal sequence can be applied (for example,
# with several L1 variable, it can be {L1.1, L1.2, L1.3} or {L1.2,L1.3,L1.1},
# without any consequences on the estimation.

# 1) define Q formulas (Qbar_L1 and Qbar_Y functions)
Q_formulas.death <- c(L1 = "Q.kplus1 ~ L0_male + L0_soc_env + A0_PM2.5",
                      Y_death = "Q.kplus1 ~ L0_male + L0_soc_env + L1 +
                                 A0_PM2.5 * M_diabetes")
Q_formulas.qol <- c(L1 = "Q.kplus1 ~ L0_male + L0_soc_env + A0_PM2.5",
                    Y_qol = "Q.kplus1 ~ L0_male + L0_soc_env + L1 +
                             A0_PM2.5 * M_diabetes")
# 2) define g formulas (needed for the ltmle package) but they are not used
#    with the g-computation estimator
g_formulas <- c("A0_PM2.5 ~ L0_male + L0_soc_env",
                "M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5 + L1")


# arguments:
#  - Anodes: indicate the exposure and the mediator variables
#  - Lnodes: indicate the intermediate confounders (+/- baseline confounders)
#  - Cnodes: censoring nodes, useless in our example
#  - Ynodes: outcome variable
#  - survivalOutcome = FALSE in our example
#  - abar: list of the two values used to define counterfactual outcomes
#          for the contrast of interest. For example, setting M=0,
#          CDE(M=0) = E(Y_{A=1,M=0}) - E(Y_{A=0,M=0})
#  - rule: to define dynamic rules (useless in our example)
#  - gbounds = c(0.01, 1) by default. This parameter is not used with g-computation
#  - Yrange = NULL, can be used to define range (min,max) for continuous outcomes
#  - SL.library = "glm",  will apply main terms glm models.
#                 The argument can be used to specify SuperLearner libraries.
#                 However, simple glm models might be preferable as data.adaptive
#                 algorithms rely on cross-validation, which is difficult and long to
#                 implement with the bootstrap procedure needed for 95% confidence
#                 intervals
#  - stratify = FALSE by default. If TRUE, glm estimations are stratified for
#               each counterfactual scenario defined in abar.
#  - estimate.time = FALSE. If TRUE, print a rough estimate of computation time
#  - iptw.only = FALSE, useless with g-computation
#  - variance.method = "ic", computation is faster than with "tmle" which
#  -                     is useless with g-comp: variance estimates rely on
#                        influence curves which cannot be used with g-comp because
#                       g-computation is not a asymptotically efficient estimator.
#  - observation.weights = NULL, can be used to specify individual weights

# with binary outcome, CDE(M=1) = P(Y_{A=1,M=0} = 1) - P(Y_{A=0,M=0} = 1)
ltmle.gcomp.CDE.M0 <- ltmle(data = df.death,
                            Anodes = c("A0_PM2.5", "M_diabetes"),
                            Lnodes = c("L1"),
                            Ynodes = c("Y_death"), # binary outcome
                            survivalOutcome = FALSE,
                            Qform = Q_formulas.death, # Q formulas
                            gform = g_formulas, # g formulas
                            abar = list(c(1,0),c(0,0)), # Y_{A=1,M=0} vs Y_{A=0,M=0}
                            rule = NULL,
                            gbounds = c(0.01, 1), # by default
                            Yrange = NULL,
                            deterministic.g.function = NULL,
                            stratify = FALSE,
                            SL.library = "glm",
                            SL.cvControl = list(),
                            estimate.time = FALSE,
                            gcomp = TRUE, # should be TRUE for g-computation
                            iptw.only = FALSE,
                            deterministic.Q.function = NULL,
                            variance.method = "ic",
                            observation.weights = NULL,
                            id = NULL)
summary(ltmle.gcomp.CDE.M0)
# Additive Treatment Effect:
#   Parameter Estimate:  0.063428 # same as manual computation
#    Estimated Std Err:  0.018159
#              p-value:  0.00047789
#    95% Conf Interval: (0.027836, 0.09902) # those 95%CI should not be used
#                      => apply a bootstrap computation instead

# with continuous outcome, CDE(M=1) = E(Y_{A=1,M=1}) - E(Y_{A=0,M=1})
ltmle.gcomp.CDE.M1 <- ltmle(data = df.qol,
                            Anodes = c("A0_PM2.5", "M_diabetes"),
                            Lnodes = c("L1"),
                            Ynodes = c("Y_qol"), # continous outcome
                            survivalOutcome = FALSE,
                            Qform = Q_formulas.qol, # Q formulas
                            gform = g_formulas, # g formulas
                            abar = list(c(1,1),c(0,1)), # Y_{A=1,M=1} vs Y_{A=0,M=1}
                            rule = NULL,
                            gbounds = c(0.01, 1), # by default
                            Yrange = NULL,
                            deterministic.g.function = NULL,
                            stratify = FALSE,
                            SL.library = "glm",
                            SL.cvControl = list(),
                            estimate.time = FALSE,
                            gcomp = TRUE, # should be TRUE for g-computation
                            iptw.only = FALSE,
                            deterministic.Q.function = NULL,
                            variance.method = "ic",
                            observation.weights = NULL,
                            id = NULL)
summary(ltmle.gcomp.CDE.M1)
# Additive Treatment Effect:
#   Parameter Estimate:  -10.432
#    Estimated Std Err:  0.55975
#              p-value:  <2e-16
#    95% Conf Interval: (-11.529, -9.335) those 95%CI should not be used
#                      => apply a bootstrap computation instead

# For quantitative outcomes, the outcome is first transformed into a continuous variable
# with [0;1] range: Y' = (Y - min(Y)) / (max(Y) - min(Y)) to run a quasi-binomial
# regression, and then estimations are back-transformed on the original scale.

# ---------------------------------------------------------------------------- #
## II.3) Sequential g-estimator ------------------------------------------------
# ---------------------------------------------------------------------------- #
# Approach described for continuous outcome.
# Extension for binary outcomes using OR, in case-control studies, is described in Vansteelandt et al. Epidemiology 20(6);2009.
rm(list=ls())
df2_int <- read.csv(file = "./data/df2_int.csv")

# 1) Regress the outcome on past
Y.qol.model <- glm(Y_qol ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                     M_diabetes + A0_PM2.5:M_diabetes,
                   family = "gaussian", data = df2_int)

# 2) Calculate a residual outcome Y - (coef.M * M_diabetes) - (coef.A0:M * A0:M)
Y.res <- (df2_int$Y_qol -
            (Y.qol.model$coefficients["M_diabetes"] * df2_int$M_diabetes) -
            (Y.qol.model$coefficients["A0_PM2.5:M_diabetes"] * df2_int$A0_PM2.5 *
               df2_int$M_diabetes) )

# 3) Regress the residual outcome on the exposure A and baseline confounders L(0)
Y.res.model <- glm(Y.res ~ L0_male + L0_soc_env + A0_PM2.5,
                   family = "gaussian", data = df2_int)

# 4) Use coefficients estimated from the 1st and 2nd regression to estimate CDE:
CDE.qol.m0.seq <- (Y.res.model$coefficients["A0_PM2.5"] + 0*Y.qol.model$coefficients["A0_PM2.5:M_diabetes"])
CDE.qol.m0.seq
# -4.869509

CDE.qol.m1.seq <- (Y.res.model$coefficients["A0_PM2.5"] + 1*Y.qol.model$coefficients["A0_PM2.5:M_diabetes"])
CDE.qol.m1.seq
# -10.38144


# ---------------------------------------------------------------------------- #
# III) Estimation of Marginal Randomized Direct and Indirect Effects -----------
# ---------------------------------------------------------------------------- #
# When Natural Direct Effects and Natural Indirect Effects are identifiable (i.e. making
# the assumption that the confounder $L(1)$ of the $M-Y$ relationship is NOT affected
# by the exposure $A$ as in Causal model 1, in Figure \@ref(fig:figDAGM1)), estimation
# are based on traditional regression models as described in chapter \@ref(ChapTradRegModels).

### by hand ----
## For the example, we will use the df1_int.csv data set (with an A*M interaction
## effect on the outcome, but no intermediate confounder affected by the exposure)
rm(list=ls())
df1_int <- read.csv(file = "./data/df1_int.csv")

## 1) Regress a model of the mediator and models of the outcome (for binary and
##    continuous outcomes)
# Estimate a model of the mediator (logistic regression)
trad_m <- glm(M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env + L1, # il semble qu'il faut aussi ajuster sur L1 ici ? sur simulation, c'est mieux sans ajuster sur L1 # de plus, régression logistique => non-collapsibility cela peut jouer
              family = "binomial",
              data = df1_int)

# Estimate models of the outcome (for continuous and binary outcomes)
trad_qol_am <- lm(Y_qol ~ A0_PM2.5 + M_diabetes + A0_PM2.5:M_diabetes +
                    L0_male + L0_soc_env + L1,
                  data = df1_int)

trad_death_am <- glm(Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5:M_diabetes +
                       L0_male + L0_soc_env + L1,
                     family = "binomial",
                     data = df1_int)

## 2) Generate predicted values for every combination of A={0,1} and M={0,1}
## 2.a) Predict counterfactual probabilities of the mediator
data.Ais0 <- data.Ais1 <- df1_int
data.Ais0$A0_PM2.5 <- 0
data.Ais1$A0_PM2.5 <- 1

# Predict the counterfactual probabilities
# P(M_{A=0}|l(0),l(1)) and  P(M_{A=0}|l(0),l(1))
P.M.Ais0 <- predict(trad_m, newdata = data.Ais0, type = "response")
P.M.Ais1 <- predict(trad_m, newdata = data.Ais1, type = "response")

## 2.b) Predict counterfactual expected values of the outcomes
data.Ais0.Mis0 <- data.Ais0.Mis1 <- df1_int
data.Ais1.Mis0 <- data.Ais1.Mis1 <- df1_int

data.Ais0.Mis0$A0_PM2.5 <- 0
data.Ais0.Mis0$M_diabetes <- 0

data.Ais0.Mis1$A0_PM2.5 <- 0
data.Ais0.Mis1$M_diabetes <- 1

data.Ais1.Mis0$A0_PM2.5 <- 1
data.Ais1.Mis0$M_diabetes <- 0

data.Ais1.Mis1$A0_PM2.5 <- 1
data.Ais1.Mis1$M_diabetes <- 1

# Predict E(Y_{am} | l(0),l(1))
Q.death.A0M0 <- predict(trad_death_am, newdata = data.Ais0.Mis0, type="response")
Q.death.A0M1 <- predict(trad_death_am, newdata = data.Ais0.Mis1, type="response")
Q.death.A1M0 <- predict(trad_death_am, newdata = data.Ais1.Mis0, type="response")
Q.death.A1M1 <- predict(trad_death_am, newdata = data.Ais1.Mis1, type="response")

Q.qol.A0M0 <- predict(trad_qol_am, newdata = data.Ais0.Mis0, type="response")
Q.qol.A0M1 <- predict(trad_qol_am, newdata = data.Ais0.Mis1, type="response")
Q.qol.A1M0 <- predict(trad_qol_am, newdata = data.Ais1.Mis0, type="response")
Q.qol.A1M1 <- predict(trad_qol_am, newdata = data.Ais1.Mis1, type="response")


## 3) Plug-in the predicted values in the g-formula and estimate the population mean
PNDE.death.gcomp <- mean((Q.death.A1M0 - Q.death.A0M0) * (1 - P.M.Ais0) +
                           (Q.death.A1M1 - Q.death.A0M1) * P.M.Ais0)
# [1] 0.0638596
TNDE.death.gcomp <- mean(Q.death.A1M0 * ((1 - P.M.Ais1) - (1 - P.M.Ais0)) +
                           Q.death.A1M1 * (P.M.Ais1 - P.M.Ais0) )
# [1] 0.01044539

PNDE.qol.gcomp <- mean((Q.qol.A1M0 - Q.qol.A0M0) * (1 - P.M.Ais0) +
                           (Q.qol.A1M1 - Q.qol.A0M1) * P.M.Ais0)
# [1] -5.310172
TNDE.qol.gcomp <- mean(Q.qol.A1M0 * ((1 - P.M.Ais1) - (1 - P.M.Ais0)) +
                         Q.qol.A1M1 * (P.M.Ais1 - P.M.Ais0) )
# [1] -1.774295


### CMAverse ----
library(CMAverse)
set.seed(1234)
res_rb_param_delta <- cmest(data = df1_int,
                            model = "rb", # for "regression based" (rb) approach
                            outcome = "Y_qol",        # outcome variable
                            exposure = "A0_PM2.5",    # exposure variable
                            mediator = "M_diabetes",  # mediator
                            basec = c("L0_male",      # confounders
                                      "L0_soc_env",
                                      "L1"),
                            EMint = TRUE, # exposures*mediator interaction
                            mreg = list("logistic"), # model of the mediator
                            yreg = "linear",       # model of the outcome
                            astar = 0,
                            a = 1,
                            mval = list(0),
                            estimation = "imputation", #  closed-form parameter
                            # function estimation
                            inference = "bootstrap") # IC95% : "delta" or "bootstrap"
summary(res_rb_param_delta)
#                Estimate Std.error  95% CIL 95% CIU  P.val
#   cde          -3.71527   0.43363 -4.70609  -2.992 <2e-16 ***
#   pnde         -5.29375   0.36149 -6.13190  -4.662 <2e-16 *** vs -5.310172
#   tnde         -5.98108   0.35602 -6.76670  -5.336 <2e-16 ***
#   pnie         -1.05652   0.14171 -1.30038  -0.784 <2e-16 ***
#   tnie         -1.74384   0.23636 -2.16969  -1.305 <2e-16 *** vs -1.774295
#   te           -7.03759   0.40445 -7.83478  -6.252 <2e-16 ***

set.seed(1234)
res_rb_param_delta <- cmest(data = df1_int,
                            model = "rb", # for "regression based" (rb) approach
                            outcome = "Y_death",        # outcome variable
                            exposure = "A0_PM2.5",    # exposure variable
                            mediator = "M_diabetes",  # mediator
                            basec = c("L0_male",      # confounders
                                      "L0_soc_env",
                                      "L1"),
                            EMint = TRUE, # exposures*mediator interaction
                            mreg = list("logistic"), # model of the mediator
                            yreg = "linear",       # model of the outcome
                            astar = 0,
                            a = 1,
                            mval = list(0),
                            estimation = "imputation", #  closed-form parameter
                            # function estimation
                            inference = "bootstrap") # IC95% : "delta" or "bootstrap"
summary(res_rb_param_delta)
#                 Estimate Std.error   95% CIL 95% CIU  P.val
#   cde           0.060001  0.016971  0.029741   0.096 <2e-16 ***
#   pnde          0.064743  0.013036  0.040779   0.094 <2e-16 *** vs 0.0638596
#   tnde          0.066644  0.012808  0.041684   0.095 <2e-16 ***
#   pnie          0.006696  0.001666  0.005635   0.012 <2e-16 ***
#   tnie          0.008597  0.003576  0.003856   0.018   0.02 *    vs 0.01044539
#   te            0.073340  0.012905  0.049302   0.102 <2e-16 ***



### mediation package ----
### MORE LIKE PARAMETRIC G-COMPUTATION based on traditional regressions +++
#### Using the "mediation" package : R Package for Causal Mediation Analysis
library(mediation)
## We will use the model of the mediator (logistic regression)
trad_m
# Call:  glm(formula = M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env +
#              L1, family = "binomial", data = df1_int)
#
# Coefficients:
#   (Intercept)     A0_PM2.5      L0_male   L0_soc_env           L1
#       -1.3788       0.5626       0.2586       0.3305       0.3346

## We will use the models of the outcome (for continuous and binary outcomes)
trad_qol_am
# Call:
#   lm(formula = Y_qol ~ A0_PM2.5 + M_diabetes + A0_PM2.5:M_diabetes +
#        L0_male + L0_soc_env + L1, data = df1_int)
#
# Coefficients:
# (Intercept) A0_PM2.5  M_diabetes  L0_male L0_soc_env       L1  A0_PM2.5:M_diabetes
#    74.7669   -3.7153     -8.6317  -0.7235    -2.8899  -3.4280              -5.6154

trad_death_am
# Call:  glm(formula = Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5:M_diabetes +
#              L0_male + L0_soc_env + L1, family = "binomial", data = df1_int)
#
# Coefficients:
# (Intercept) A0_PM2.5  M_diabetes  L0_male L0_soc_env       L1  A0_PM2.5:M_diabetes
#    -2.06294  0.36668     0.40921  0.29249    0.36360  0.44716              0.01275

## Use the mediate() function to estimate the TNDE, PNIE
# estimations relies on quasi-bayesian Monte Carlo method (especially for continuous
# mediators) => set seed for reproducibility
set.seed(2024)

## For the quantitative outcome
mediation.res.qol <- mediate(trad_m,                   # model of the mediator
                             trad_qol_am,              # model of the outcome)
                             treat = "A0_PM2.5",       # exposure
                             mediator = "M_diabetes",  # mediator
                             robustSE = TRUE,          # estimate sandwich SEs
                             sims = 100)               # better to use >= 1000
summary(mediation.res.qol)
#                            Estimate 95% CI Lower 95% CI Upper p-value
#   ACME (control)             -1.074       -1.356        -0.82  <2e-16 *** PNIE
#   ACME (treated)             -1.778       -2.239        -1.38  <2e-16 *** TNIE
#   ADE (control)              -5.301       -5.950        -4.63  <2e-16 *** PNDE
#   ADE (treated)              -6.005       -6.646        -5.27  <2e-16 *** TNDE
#   Total Effect               -7.079       -7.768        -6.24  <2e-16 ***
#   Prop. Mediated (control)    0.151        0.118         0.19  <2e-16 ***
#   Prop. Mediated (treated)    0.251        0.197         0.31  <2e-16 ***
#   ACME (average)             -1.426       -1.811        -1.11  <2e-16 ***
#   ADE (average)              -5.653       -6.267        -4.95  <2e-16 ***
#   Prop. Mediated (average)    0.201        0.160         0.25  <2e-16 ***

## We can plot the estimations
plot(mediation.res.qol)

## For the binary outcome
mediation.res.death <- mediate(trad_m,                   # model of the mediator
                               trad_death_am,            # model of the outcome)
                               treat = "A0_PM2.5",       # exposure
                               mediator = "M_diabetes",  # mediator
                               robustSE = TRUE,          # estimate sandwich SEs
                               sims = 100)               # better to use >= 1000
summary(mediation.res.death)
#                          Estimate 95% CI Lower 95% CI Upper p-value
# ACME (control)            0.00860      0.00569         0.01  <2e-16 *** PNIE
# ACME (treated)            0.00979      0.00377         0.02  <2e-16 *** TNIE
# ADE (control)             0.06543      0.03658         0.09  <2e-16 *** PNDE
# ADE (treated)             0.06662      0.03946         0.09  <2e-16 *** TNDE
# Total Effect              0.07522      0.04800         0.10  <2e-16 ***
# Prop. Mediated (control)  0.11172      0.07323         0.18  <2e-16 ***
# Prop. Mediated (treated)  0.13066      0.05714         0.22  <2e-16 ***
# ACME (average)            0.00920      0.00567         0.01  <2e-16 ***
# ADE (average)             0.06603      0.03803         0.09  <2e-16 ***
# Prop. Mediated (average)  0.12119      0.07330         0.20  <2e-16 ***

## Sensitivity analysis to test sequential ignorability (assess the possible existence
## of unobserved (baseline) confounders of the M-Y relationship)
trad_m <- glm(M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env + L1,
              family = binomial("probit"), # sensitivity analysis works only for
              # probit models
              data = df1_int)
set.seed(1234)
mediation.res.qol <- mediate(trad_m,                   # model of the mediator
                             trad_qol_am,              # model of the outcome)
                             treat = "A0_PM2.5",       # exposure
                             mediator = "M_diabetes",  # mediator
                             robustSE = TRUE,          # estimate sandwich SEs
                             # long = TRUE,
                             sims = 100)               # better to use >= 1000
summary(mediation.res.qol)
#                            Estimate 95% CI Lower 95% CI Upper p-value
#   ACME (control)             -1.075       -1.430        -0.73  <2e-16 ***
#   ACME (treated)             -1.772       -2.353        -1.14  <2e-16 ***
#   ADE (control)              -5.365       -5.961        -4.76  <2e-16 ***
#   ADE (treated)              -6.063       -6.638        -5.56  <2e-16 ***
#   Total Effect               -7.138       -7.823        -6.48  <2e-16 ***
#   Prop. Mediated (control)    0.153        0.109         0.19  <2e-16 ***
#   Prop. Mediated (treated)    0.251        0.172         0.32  <2e-16 ***
#   ACME (average)             -1.424       -1.884        -0.94  <2e-16 ***
#   ADE (average)              -5.714       -6.283        -5.14  <2e-16 ***
#   Prop. Mediated (average)    0.202        0.142         0.25  <2e-16 ***

sensisitivy <- medsens(mediation.res.qol,
                       rho.by = 0.1, # sensitivity parameter = correlation between
                       # the residuals of the mediator & outcome regressions
                       # here, rho varies from -0.9 to +0.9 by 0.1 increments
                       effect.type = "both", # "direct", "indirect" or "both"
                       sims = 100)
# If there exist unobserved (baseline) confounders of the M-Y relationship, we expect
# that rho is no longer zero.
# The sensitivity analysis is conducted by varying the value of rho and examining
# how the estimated ACME and ADE change.
summary(sensisitivy)
# Mediation Sensitivity Analysis: Average Mediation Effect
# Sensitivity Region: ACME for Control Group
#         Rho ACME(control) 95% CI Lower 95% CI Upper R^2_M*R^2_Y* R^2_M~R^2_Y~
#   [1,] -0.6        0.0474      -0.0093       0.1229         0.36       0.2647
#
# Rho at which ACME for Control Group = 0: -0.6
# R^2_M*R^2_Y* at which ACME for Control Group = 0: 0.36
# R^2_M~R^2_Y~ at which ACME for Control Group = 0: 0.2647
#
# Rho at which ACME for Treatment Group = 0: -0.9
# R^2_M*R^2_Y* at which ACME for Treatment Group = 0: 0.81
# R^2_M~R^2_Y~ at which ACME for Treatment Group = 0: 0.5956
#
# Mediation Sensitivity Analysis: Average Direct Effect
# Rho at which ADE for Control Group = 0: 0.8
# R^2_M*R^2_Y* at which ADE for Control Group = 0: 0.64
# R^2_M~R^2_Y~ at which ADE for Control Group = 0: 0.4706
#
# Rho at which ADE for Treatment Group = 0: 0.8
# R^2_M*R^2_Y* at which ADE for Treatment Group = 0: 0.64
# R^2_M~R^2_Y~ at which ADE for Treatment Group = 0: 0.4706

par(mfrow = c(2,2))
plot(sensisitivy)
par(mfrow = c(2,1))

# Here, for PNIE, the confidence interval of the ACME (PNIE) contain zero when rho=-0.6
# R^2 correspond to
# When the product of the residual variance explained by the omitted confounding
# is 0.36, the point estimate of PNIE = 0
# When the product of the total variance explained by the omitted confounding
# is 0.27, the point estimate of PNIE = 0
# => The risk that PNIE is cancelled because of unmeasured M-Y confounding
#    seems low in our example (it requires rather strong correlations)

# Note: sensitivity analysis from the package cannot be applied when the mediator
#       and the outcome are both binary


# ---------------------------------------------------------------------------- #
# IV) Estimation of Marginal Randomized Direct and Indirect Effects -----------
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
## IV.1) parametric g-computation ---------------------------------------------
# ---------------------------------------------------------------------------- #
# described in Lin et al. Epidemiology 2017;28(2):266-74
# The approach is described as an adaptation of the parametric g-computation presented for controlled direct effects.

# The following steps are applied:
# 1) Fit parametric models for the observed data for the time-varying confounders L(1), the mediator M and the outcome Y

# 2) Estimate the joint distribution of time-varying
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")

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
  L1.model <- glm(L1 ~ L0_male + L0_soc_env + A0_PM2.5,
                  family = "binomial", data = df2_int)
  M.model <- glm(M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                 family = "binomial", data = df2_int)
  ### 1b) fit parametric models of the outcomes given the past
  Y.death.model <- glm(Y_death ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                         M_diabetes + A0_PM2.5:M_diabetes,
                       family = "binomial", data = df2_int)
  Y.qol.model <- glm(Y_qol ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                       M_diabetes + A0_PM2.5:M_diabetes,
                     family = "gaussian", data = df2_int)

  ## 2) Estimate the joint distribution of time-varying confounders and of the
  ##    mediator under the counterfactual scenarios setting A0_PM2.5 = 1 or 0

  # set the exposure A0_PM2.5 to 0 or 1 in two new counterfactual data sets
  data.A0  <- data.A1 <- df2_int
  data.A0$A0_PM2.5 <- 0
  data.A1$A0_PM2.5 <- 1
  # simulate L1 values under the counterfactual exposures A0_PM2.5=0 or A0_PM2.5=1
  p.L1.A0 <- predict(L1.model, newdata = data.A0, type="response")
  p.L1.A1 <- predict(L1.model, newdata = data.A1, type="response")
  sim.L1.A0 <- rbinom(n = nrow(df2_int), size = 1, prob = p.L1.A0)
  sim.L1.A1 <- rbinom(n = nrow(df2_int), size = 1, prob = p.L1.A1)

  # replace L(1) by their counterfactual values in the data under A=0 or A=1
  data.A0.L <- data.A0
  data.A1.L <- data.A1
  data.A0.L$L1 <- sim.L1.A0
  data.A1.L$L1 <- sim.L1.A1

  # simulate M values under the counterfactual exposures A0_PM2.5=0 or A0_PM2.5=1
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

  data.A0.G0$M_diabetes <- marg.M.A0
  # data.A0.G1$M_diabetes <- marg.M.A1 # note: this data set will not be useful

  data.A1.G0$M_diabetes <- marg.M.A0
  data.A1.G1$M_diabetes <- marg.M.A1

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
## IV.2) G-computation by iterative conditional expectation--------------------
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
G.model <- glm(M_diabetes ~ L0_male + L0_soc_env + A0_PM2.5,
               family = "binomial", data = df2_int)

### 1b) generate predicted probabilites by evaluating the regression setting the
###     exposure value to A=0 or to A=1
# create datasets corresponding to the counterfactual scenarios setting A=0 and A=1
data.Ais0  <- data.Ais1 <- df2_int
data.Ais0$A0_PM2.5 <- 0
data.Ais1$A0_PM2.5 <- 1

# estimate G_{A=0|L(0)} = Pr(M=1|A=0,L(0)) and G_{A=1|L(0)} = Pr(M=1|A=1,L(0))
G.Ais0.L0 <-predict(G.model, newdata = data.Ais0, type="response")
G.Ais1.L0 <-predict(G.model, newdata = data.Ais1, type="response")


## 2) Fit parametric models for the observed data for the outcome Y given the past
##    and generate predicted values by evaluating the regression setting the mediator
##    value to M=0 or to M=1
##    then calculate a weighted sum of the predicted Q.L2, with weights given by G
### 2a) fit parametric models of the outcomes given the past
Y.death.model <- glm(Y_death ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                       M_diabetes + A0_PM2.5:M_diabetes,
                     family = "binomial", data = df2_int)
Y.qol.model <- glm(Y_qol ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                     M_diabetes + A0_PM2.5:M_diabetes,
                   family = "gaussian", data = df2_int)

### 2b) generate predicted values by evaluating the regression setting the mediator
###     value to M=0 or to M=1
data.Mis0  <- data.Mis1 <- df2_int
data.Mis0$M_diabetes <- 0
data.Mis1$M_diabetes <- 1

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
L1.death.A0.G0.model <- glm(Q.L2.death.A0.G0 ~ L0_male + L0_soc_env + A0_PM2.5,
                            family = "quasibinomial", data = df2_int)
L1.death.A1.G0.model <- glm(Q.L2.death.A1.G0 ~ L0_male + L0_soc_env + A0_PM2.5,
                            family = "quasibinomial", data = df2_int)
L1.death.A1.G1.model <- glm(Q.L2.death.A1.G1 ~ L0_male + L0_soc_env + A0_PM2.5,
                            family = "quasibinomial", data = df2_int)

L1.qol.A0.G0.model <- glm(Q.L2.qol.A0.G0 ~ L0_male + L0_soc_env + A0_PM2.5,
                          family = "gaussian", data = df2_int)
L1.qol.A1.G0.model <- glm(Q.L2.qol.A1.G0 ~ L0_male + L0_soc_env + A0_PM2.5,
                          family = "gaussian", data = df2_int)
L1.qol.A1.G1.model <- glm(Q.L2.qol.A1.G1 ~ L0_male + L0_soc_env + A0_PM2.5,
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
# V) CMAverse ----
# ---------------------------------------------------------------------------- #
# The DAG for this scientific setting is:
# devtools::install_github("BS1125/CMAverse")
library(CMAverse)
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")
cmdag(outcome = "Y_death", exposure = "A0_PM2.5", mediator = "M_diabetes",
      basec = c("L0_male", "L0_soc_env"), postc = "L1", node = TRUE, text_col = "white")
# The CMAverse package can be used to estimate CDE, rNDE and rNIE by parametric g-computation
# (using the 'gformula' option in the )

### Continuous outcome
## The g-formula Approach
set.seed(1234)
res_gformula_Qol <- cmest(data = df2_int,
                         model = "gformula", # for parametric g-computation
                         outcome = "Y_qol", # outcome variable
                         exposure = "A0_PM2.5", # exposure variable
                         mediator = "M_diabetes", # mediator
                         basec = c("L0_male",
                                   "L0_soc_env"), # confounders
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
#   glm(formula = Y_qol ~ A0_PM2.5 + M_diabetes + A0_PM2.5 * M_diabetes +
#         L0_male + L0_soc_env + L1, family = gaussian(),
#       data = getCall(x$reg.output$yreg)$data, weights = getCall(x$reg.output$yreg)$weights)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)            74.8247     0.2133 350.823  < 2e-16 ***
#   A0_PM2.5                 -3.7014     0.4295  -8.617  < 2e-16 ***
#   M_diabetes              -8.6336     0.2331 -37.042  < 2e-16 ***
#   L0_male                -0.7280     0.2019  -3.605 0.000313 ***
#   L0_soc_env             -2.8828     0.2116 -13.621  < 2e-16 ***
#   L1                     -5.1668     0.2189 -23.608  < 2e-16 ***
#   A0_PM2.5:M_diabetes       -5.5119     0.6440  -8.559  < 2e-16 ***

### 2) Estimation of g(M=1|L1,A,L0), model of the mediator
### Mediator regressions:
# Call:
#   glm(formula = M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env +
#         L1, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -1.36249    0.04783 -28.488  < 2e-16 ***
#   A0_PM2.5                 0.30994    0.06668   4.648 3.35e-06 ***
#   L0_male                0.24661    0.04369   5.644 1.66e-08 ***
#   L0_soc_env  0.30628    0.04650   6.587 4.50e-11 ***
#   L1                     0.86045    0.04493  19.152  < 2e-16 ***

### 3) Estimation of Qbar.L1 = P(L1=1|A,L0), model of intermediate confounder
### Regressions for mediator-outcome confounders affected by the exposure:
# Call:
#   glm(formula = L1 ~ A0_PM2.5 + L0_male + L0_soc_env,
#       family = binomial(), data = getCall(x$reg.output$postcreg[[1L]])$data,
#       weights = getCall(x$reg.output$postcreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -0.86983    0.04292 -20.267  < 2e-16 ***
#   A0_PM2.5                 0.94354    0.06475  14.572  < 2e-16 ***
#   L0_male               -0.19827    0.04289  -4.622 3.80e-06 ***
#   L0_soc_env  0.32047    0.04556   7.034 2.01e-12 ***

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

### We can check those equalities for the analogues of the 3-way and 4-way decomposition
res_gformula_Qol$effect.pe["te"]
# TE = -9.870139
res_gformula_Qol$effect.pe["rpnde"]
# rPNDE = -7.565835
res_gformula_Qol$effect.pe["rpnie"]
# rPNIE = -1.40641
res_gformula_Qol$effect.pe["cde"]
# CDE(M=0) = -5.86375

## Check that MI = TE - rPNDE - rPNIE
(res_gformula_Qol$effect.pe["te"] - res_gformula_Qol$effect.pe["rpnde"] -
    res_gformula_Qol$effect.pe["rpnie"])
# -0.897894
## Check that MI = rTNIE - rPNIE = rTNDE - rPNDE
(res_gformula_Qol$effect.pe["rtnie"] - res_gformula_Qol$effect.pe["rpnie"])
(res_gformula_Qol$effect.pe["rtnde"] - res_gformula_Qol$effect.pe["rpnde"])
# -0.897894
res_gformula_Qol$effect.pe["rintmed"]
# -0.897894 # we have MI = TE - rPNDE - rPNIE = rTNIE - rPNIE = rTNDE - rPNDE

## Check that RE = PNDE - CDE_{M=0}
res_gformula_Qol$effect.pe["rpnde"] - res_gformula_Qol$effect.pe["cde"]
# -1.702085
res_gformula_Qol$effect.pe["rintref"]
# -1.702085 # we have RE = PNDE - CDE_{M=0}


## The g-formula Approach
set.seed(1234)
res_gformula_OR <- cmest(data = df2_int,
                         # data.frame(df2_int[,c("L0_male","L0_soc_env","A0_PM2.5")],
                         #                L1=as.factor(df2_int$L1),
                         #                df2_int[,c("M_diabetes","Y_death")]),
                      model = "gformula",
                      outcome = "Y_death",
                      exposure = "A0_PM2.5",
                      mediator = "M_diabetes",
                      basec = c("L0_male", "L0_soc_env"),
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
#   glm(formula = Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5 * M_diabetes +
#         L0_male + L0_soc_env + L1, family = binomial(),
#       data = getCall(x$reg.output$yreg)$data, weights = getCall(x$reg.output$yreg)$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -2.04033    0.05855 -34.849  < 2e-16 ***
#   A0_PM2.5                 0.28922    0.10123   2.857  0.00428 **
#   M_diabetes              0.44406    0.05569   7.974 1.53e-15 ***
#   L0_male                0.26913    0.04996   5.387 7.15e-08 ***
#   L0_soc_env  0.34603    0.05432   6.370 1.89e-10 ***
#   L11                    0.42894    0.05195   8.257  < 2e-16 ***
#   A0_PM2.5:M_diabetes       0.04387    0.14311   0.307  0.75919
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
#   glm(formula = M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env +
#         L1, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -1.36249    0.04783 -28.488  < 2e-16 ***
#   A0_PM2.5                 0.30994    0.06668   4.648 3.35e-06 ***
#   L0_male                0.24661    0.04369   5.644 1.66e-08 ***
#   L0_soc_env  0.30628    0.04650   6.587 4.50e-11 ***
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
#   glm(formula = L1 ~ A0_PM2.5 + L0_male + L0_soc_env,
#       family = binomial(), data = getCall(x$reg.output$postcreg[[1L]])$data,
#       weights = getCall(x$reg.output$postcreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -0.86983    0.04292 -20.267  < 2e-16 ***
#   A0_PM2.5                 0.94354    0.06475  14.572  < 2e-16 ***
#   L0_male               -0.19827    0.04289  -4.622 3.80e-06 ***
#   L0_soc_env  0.32047    0.04556   7.034 2.01e-12 ***
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
res_gformula_RD_M0 <- cmest(data = df2_int,
                            # data.frame(df2_int[,c("L0_male","L0_soc_env","A0_PM2.5")],
                            #                   L1=as.factor(df2_int$L1),
                            #                   df2_int[,c("M_diabetes","Y_death")]),
                            model = "gformula",
                            outcome = "Y_death",
                            exposure = "A0_PM2.5",
                            mediator = "M_diabetes",
                            basec = c("L0_male", "L0_soc_env"),
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
#   glm(formula = Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5 * M_diabetes +
#         L0_male + L0_soc_env + L1, family = gaussian(),
#       data = getCall(x$reg.output$yreg)$data, weights = getCall(x$reg.output$yreg)$weights)
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
#   (Intercept)           0.100674   0.008574  11.742  < 2e-16 ***
#   A0_PM2.5                0.046456   0.017268   2.690  0.00715 **
#   M_diabetes             0.074378   0.009370   7.938 2.27e-15 ***
#   L0_male               0.043661   0.008117   5.379 7.67e-08 ***
#   L0_soc_env 0.053711   0.008508   6.313 2.86e-10 ***
#   L11                   0.073753   0.008798   8.383  < 2e-16 ***
#   A0_PM2.5:M_diabetes      0.030159   0.025887   1.165  0.24404

## 2) Mediator regressions:  g(M=1|L1,A,L0)
# Call:
#   glm(formula = M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env +
#         L1, family = binomial(), data = getCall(x$reg.output$mreg[[1L]])$data,
#       weights = getCall(x$reg.output$mreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -1.36249    0.04783 -28.488  < 2e-16 ***
#   A0_PM2.5                 0.30994    0.06668   4.648 3.35e-06 ***
#   L0_male                0.24661    0.04369   5.644 1.66e-08 ***
#   L0_soc_env  0.30628    0.04650   6.587 4.50e-11 ***
#   L11                    0.86045    0.04493  19.152  < 2e-16 ***

## 3) Regressions for mediator-outcome confounders affected by the exposure: Qbar.L1 = P(L1=1|A,L0)
# Call:
#   glm(formula = L1 ~ A0_PM2.5 + L0_male + L0_soc_env,
#       family = binomial(), data = getCall(x$reg.output$postcreg[[1L]])$data,
#       weights = getCall(x$reg.output$postcreg[[1L]])$weights)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -0.86983    0.04292 -20.267  < 2e-16 ***
#   A0_PM2.5                 0.94354    0.06475  14.572  < 2e-16 ***
#   L0_male               -0.19827    0.04289  -4.622 3.80e-06 ***
#   L0_soc_env  0.32047    0.04556   7.034 2.01e-12 ***

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
res_gformula_RD_M1 <- cmest(data = data.frame(df2_int[,c("L0_male","L0_soc_env","A0_PM2.5")],
                                              L1=as.factor(df2_int$L1),
                                              df2_int[,c("M_diabetes","Y_death")]),
                            model = "gformula",
                            outcome = "Y_death",
                            exposure = "A0_PM2.5",
                            mediator = "M_diabetes",
                            basec = c("L0_male", "L0_soc_env"),
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
res_gformula_QoL_M0 <- cmest(data = data.frame(df2_int[,c("L0_male","L0_soc_env","A0_PM2.5")],
                                              L1=as.factor(df2_int$L1),
                                              df2_int[,c("M_diabetes","Y_qol")]),
                            model = "gformula",
                            outcome = "Y_qol",
                            exposure = "A0_PM2.5",
                            mediator = "M_diabetes",
                            basec = c("L0_male", "L0_soc_env"),
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
res_gformula_QoL_M1 <- cmest(data = data.frame(df2_int[,c("L0_male","L0_soc_env","A0_PM2.5")],
                                               L1=as.factor(df2_int$L1),
                                               df2_int[,c("M_diabetes","Y_qol")]),
                             model = "gformula",
                             outcome = "Y_qol",
                             exposure = "A0_PM2.5",
                             mediator = "M_diabetes",
                             basec = c("L0_male", "L0_soc_env"),
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



# ---------------------------------------------------------------------------- #
# VI) Estimation of Conditional Randomized Direct and Indirect Effects ----
# ---------------------------------------------------------------------------- #
rm(list=ls())
df2_int <- read.csv(file = "data/df2_int.csv")

## 1) Q^{a,a'}_R2 =  Y
Q.death_R2 <- df2_int$Y_death
Q.qol_R2 <- df2_int$Y_qol

## 2) Obtain Q^{a,a'}_L1, Q^{a,a'}_M1, Q^{a,a'}_R1
## 2.a) Obtain Q^{a,a'}_L1
# Regress Q^{a,a'}_R2 on observed values (L(0),A,L(1),M)
L1.model.death <- glm(Q.death_R2 ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                        M_diabetes + A0_PM2.5:M_diabetes,
                      family = "binomial", data = df2_int)
L1.model.qol <- glm(Q.qol_R2 ~ L0_male + L0_soc_env + A0_PM2.5 + L1 +
                      M_diabetes + A0_PM2.5:M_diabetes,
                    family = "gaussian", data = df2_int)

# Evaluate the fitted function at the observed mediator M and covariate history L(1),L(0)
# and the intervened exposure A = a
data.Ais0 <- data.Ais1 <- df2_int
data.Ais0$A0_PM2.5 <- 0
data.Ais1$A0_PM2.5 <- 1

# We will need 3 counterfactual quantities Q^{a,a'}: Q^{1,1}, Q^{1,0} and Q^{0,0}
Q11.death.L1 <- predict(L1.model.death, newdata = data.Ais1, type = "response")
Q10.death.L1 <- predict(L1.model.death, newdata = data.Ais1, type = "response")
Q00.death.L1 <- predict(L1.model.death, newdata = data.Ais0, type = "response")

Q11.qol.L1 <- predict(L1.model.qol, newdata = data.Ais1, type = "response")
Q10.qol.L1 <- predict(L1.model.qol, newdata = data.Ais1, type = "response")
Q00.qol.L1 <- predict(L1.model.qol, newdata = data.Ais0, type = "response")

## 2.b) Obtain Q^{a,a'}_M1
# Regress Q^{a,a'}_L1 on observed values (L(0),A,L(1))
M1.model.death.11 <- glm(Q11.death.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                         family = "quasibinomial", data = df2_int)
M1.model.death.10 <- glm(Q10.death.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                         family = "quasibinomial", data = df2_int)
M1.model.death.00 <- glm(Q00.death.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                         family = "quasibinomial", data = df2_int)

M1.model.qol.11 <- glm(Q11.qol.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                       family = "gaussian", data = df2_int)
M1.model.qol.10 <- glm(Q10.qol.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                       family = "gaussian", data = df2_int)
M1.model.qol.00 <- glm(Q00.qol.L1 ~ L0_male + L0_soc_env + A0_PM2.5 + L1,
                       family = "gaussian", data = df2_int)

# Evaluate the fitted function at the observed covariate history L(1),L(0)
# and the intervened exposure A = a'
Q11.death.M1 <- predict(M1.model.death.11, newdata = data.Ais1, type = "response")
Q10.death.M1 <- predict(M1.model.death.10, newdata = data.Ais0, type = "response")
Q00.death.M1 <- predict(M1.model.death.00, newdata = data.Ais0, type = "response")

Q11.qol.M1 <- predict(M1.model.qol.11, newdata = data.Ais1, type = "response")
Q10.qol.M1 <- predict(M1.model.qol.10, newdata = data.Ais0, type = "response")
Q00.qol.M1 <- predict(M1.model.qol.00, newdata = data.Ais0, type = "response")

## 2.c) Obtain Q^{a,a'}_R1
# Regress Q^{a,a'}_M1 on observed values (L(0),A)
R1.model.death.11 <- glm(Q11.death.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "quasibinomial", data = df2_int)
R1.model.death.10 <- glm(Q10.death.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "quasibinomial", data = df2_int)
R1.model.death.00 <- glm(Q00.death.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                         family = "quasibinomial", data = df2_int)

R1.model.qol.11 <- glm(Q11.qol.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                       family = "gaussian", data = df2_int)
R1.model.qol.10 <- glm(Q10.qol.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                       family = "gaussian", data = df2_int)
R1.model.qol.00 <- glm(Q00.qol.M1 ~ L0_male + L0_soc_env + A0_PM2.5,
                       family = "gaussian", data = df2_int)

# Evaluate the fitted function at the observed covariate history L(0)
# and the intervened exposure A = a
Q11.death.R1 <- predict(R1.model.death.11, newdata = data.Ais1, type = "response")
Q10.death.R1 <- predict(R1.model.death.10, newdata = data.Ais1, type = "response")
Q00.death.R1 <- predict(R1.model.death.00, newdata = data.Ais0, type = "response")

Q11.qol.R1 <- predict(R1.model.qol.11, newdata = data.Ais1, type = "response")
Q10.qol.R1 <- predict(R1.model.qol.10, newdata = data.Ais1, type = "response")
Q00.qol.R1 <- predict(R1.model.qol.00, newdata = data.Ais0, type = "response")

## 4) Estimate the Conditional Randomized Direct and Indirect Effects
## For death
CRDE.death <- mean(Q10.death.R1) - mean(Q00.death.R1)
CRDE.death
# [1] 0.07585836
CRIE.death <- mean(Q11.death.R1) - mean(Q10.death.R1)
CRIE.death
# [1] 0.006802907

## For quality of life
CRDE.qol <- mean(Q10.qol.R1) - mean(Q00.qol.R1)
CRDE.qol
# [1] -7.278295
CRIE.qol <- mean(Q11.qol.R1) - mean(Q10.qol.R1)
CRIE.qol
# [1] -1.015901
