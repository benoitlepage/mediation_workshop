###

rm(list=ls())

df1 <- read.csv(file = "data/df1.csv")
df1_int <- read.csv(file = "data/df1_int.csv")

# Baron & Kenny ----

## Model 1 to estimate the total effect:
model.tot.A.QoL <- lm(Y_qol ~ A0_PM2.5 + L0_male + L0_soc_env,
                      data = df1)
summary(model.tot.A.QoL)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
#   (Intercept)  71.8820     0.2155 333.565  < 2e-16 ***
#   A0_PM2.5     -5.0961     0.3486 -14.617  < 2e-16 *** <- Total effect
#   L0_male      -1.1486     0.2194  -5.235 1.68e-07 ***
#   L0_soc_env   -3.4441     0.2295 -15.005  < 2e-16 ***

# The total effect of being exposed to high levels of PM_2.5 on Quality of life
# is approximately equal to an average decrease of -5.1 on the QoL scale,
# given by the A0_PM2.5 coefficient:
model.tot.A.QoL$coefficients["A0_PM2.5"]
# -5.096057

## Model 2 to estimate the effect of the exposure on the mediator
## because the mediator is binary, we might want to use a logistic or probit regression
## for example
logit.model.A.M <- glm(M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env,
                       data = df1, family = "binomial")
summary(logit.model.A.M)
# effects estimated on the logit scale:
# Coefficients:
#                          Estimate Std. Error z value Pr(>|z|)
#   (Intercept) -1.27152    0.04584 -27.736  < 2e-16 ***
#   A0_PM2.5     0.56168    0.06537   8.592  < 2e-16 *** <- effect of A on M
#   L0_male      0.25455    0.04425   5.753 8.77e-09 ***
#   L0_soc_env   0.32683    0.04731   6.908 4.91e-12 ***
exp(coefficients(logit.model.A.M)["A0_PM2.5"])
# Odds ratio = 1.753609  for the effect of being exposed to high levels of PM_2.5
# on the mediator (probability of type 2 diabetes)

## Model 3 to estimate the direct effect of the exposure (conditional on the outcome) and
## the effect of M on Y, adjusted for confounders of the A-Y and M-Y relationships
model.A.M.QoL <- lm(Y_qol ~ A0_PM2.5 + M_diabetes + L1 + L0_male + L0_soc_env,
                    data = df1)
summary(model.A.M.QoL)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
#   (Intercept)  74.7858     0.2130 351.178  < 2e-16 ***
#   A0_PM2.5     -3.9650     0.3212 -12.345  < 2e-16 *** <- Direct effect
#   M_diabetes   -8.7138     0.2221 -39.237  < 2e-16 *** <- effect of M on Y
#   L1           -3.4252     0.2212 -15.483  < 2e-16 ***
#   L0_male      -0.7193     0.2017  -3.566 0.000364 ***
#   L0_soc_env   -2.8876     0.2112 -13.674  < 2e-16 ***

# The direct effect of PM_2.5 is approximately -4.0 given by the A0_PM2.5 coefficient
model.A.M.QoL$coefficients["A0_PM2.5"]
# -3.965038

## Following the Baron & Kenny Steps, we would conclude that :
# - There is a significant total effect of PM_2.5 on Quality of Life (Model 1)
# - There is a significant effect of PM_2.5 on the mediator (diabetes) (Model 2)
# - There is a significant effect of the mediator (diabetes) on Qol (model 3)
# - The direct effect is significantly non-null
# => Conclusion: Diabetes partially mediates the relationship between PM_2.5 and QoL

### Estimation of the indirect effect:
### We can apply the difference in coefficient method to estimate the indirect effect:
### substract the direct effect from the Total effect:
ind.effect.dif.meth <- (model.tot.A.QoL$coefficients["A0_PM2.5"] -
                          model.A.M.QoL$coefficients["A0_PM2.5"])

# -1.131019
# the indirect effect is approximately -1.1
# The confidence interval of the indirect effect can be computed by bootstrap.

# Because the mediator is binary and we applied a logistic regression for Model 2,
# we cannot apply the product of coefficients combining a coefficient from
# Model 2 (logit scale) and from Model 3 (difference scale)

### Model 2bis
# Surprisingly, another possibility is to run a linear model of the binary mediator
# instead of the logistic regression to apply the "product of coefficient method"
# in order to estimate the indirect effect:
linear.model.A.M <- lm(M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env,
                       data = df1)
summary(linear.model.A.M)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.215398   0.008929  24.123  < 2e-16 ***
# A0_PM2.5    0.127180   0.014446   8.804  < 2e-16 ***
# L0_male     0.052363   0.009091   5.760 8.67e-09 ***
# L0_soc_env  0.065806   0.009511   6.919 4.83e-12 ***

## product of coefficient method:
ind.effect.prod.meth <- (linear.model.A.M$coefficients["A0_PM2.5"] *
                           model.A.M.QoL$coefficients["M_diabetes"])
# -1.108213
# which also gives an indirect effect of approximately -1.1

### The Baron & Kenny approach is usually applied for continuous outcomes,
### using linear regressions. It is less adapted for binary outcomes.
### However, as for the binary mediator, using linear regression of the mediator
### and the outcome could still give some results.

### For binary outcomes, the Baron & Kenny approach
## Model 1: linear model of the probability of death to estimate the total effect:
model.tot.A.death <- lm(Y_death ~ A0_PM2.5 + L0_male + L0_soc_env,
                        data = df1)
summary(model.tot.A.death)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
#   (Intercept) 0.135282   0.007909  17.104  < 2e-16 ***
#   A0_PM2.5    0.060247   0.012796   4.708 2.53e-06 *** <- Total effect
#   L0_male     0.050285   0.008053   6.244 4.43e-10 ***
#   L0_soc_env  0.059565   0.008425   7.070 1.65e-12 ***
# On a risk difference scale the total effect of being exposed to high levels of
# PM_2.5 on the probability of death is approximately +6.0%

## Model 3: linear model to estimate the direct effect of the exposure (conditional
## on the outcome) and the effect of M on Y, adjusted for confounders of the A-Y
## and M-Y relationships
model.A.M.death <- lm(Y_death ~ A0_PM2.5 + M_diabetes + L1 + L0_male + L0_soc_env,
                      data = df1)
summary(model.A.M.death)
#   Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)
#   (Intercept) 0.098691   0.008460  11.666  < 2e-16 ***
#   A0_PM2.5    0.051509   0.012759   4.037 5.45e-05 *** <- Direct effect
#   M_diabetes  0.064751   0.008822   7.340 2.31e-13 *** <- effect of M on Y
#   L1          0.075533   0.008788   8.595  < 2e-16 ***
#   L0_male     0.047490   0.008013   5.927 3.19e-09 ***
#   L0_soc_env  0.055676   0.008389   6.637 3.36e-11 ***

# The direct effect is approximately +5.2% given by the A0_PM2.5 coefficient
model.A.M.death$coefficients["A0_PM2.5"]
# 0.05150901

# The indirect effect can be calculated by the "difference in coefficient" method
# using coefficients from models 1 and 3
model.tot.A.death$coefficients["A0_PM2.5"] - model.A.M.death$coefficients["A0_PM2.5"]
# 0.008737889, i.e. approximately 0.9%

# or the product of coefficients using the previous model 2bis and model 3:
linear.model.A.M$coefficients["A0_PM2.5"] * model.A.M.death$coefficients["M_diabetes"]
# 0.008234978, i.e. approximately 0.8%


# SEM ----
library(lavaan)
library(semPlot)

## quality of life ----
df1.qol.sem <- data.frame(L0_male = df1$L0_male,
                          L0_soc_env = df1$L0_soc_env,
                          A0_PM2.5 = ordered(df1$A0_PM2.5), #df1$A0_PM2.5, # better not to declare
                          L1 = df1$L1,
                          M_diabetes = ordered(df1$M_diabetes), #df1$M_diabetes, # better not to declare
                          Y_qol = df1$Y_qol)
sem.QoL <- "                                 # models are written between quotes
## Regression models
## we can also add the names of some path coefficients
  A0_PM2.5 ~ L0_male + L0_soc_env
  M_diabetes ~ L0_male + L0_soc_env + b.A * A0_PM2.5 + L1   # label b.A path coef
  Y_qol ~ L0_male + L0_soc_env + c.A*A0_PM2.5 + L1 + c.M*M_diabetes # label c.A and c.M

## Covariances
# Assuming a non-null covariance between confounders
# covariances are represented with a double tilde ~~ (double arrow)
# (note: in the data-generating system, the null assumptions were true)
  L0_male ~~ L0_soc_env
  L0_male ~~ L1
  L0_soc_env ~~ L1

## We can define other parameters using the := syntax
## we want the direct, indirect and total effects:
  direct := c.A
  indirect := b.A * c.M
  total := (b.A * c.M) + c.A
"

fit.qol <- lavaan::sem(model = sem.QoL,
                       data = df1)
summary(fit.qol)
# lavaan 0.6-18 ended normally after 23 iterations
#
# Estimator                                         ML
# Optimization method                           NLMINB
# Number of model parameters                        20
#
# Number of observations                         10000
#
# Model Test User Model:
#
#   Test statistic                                 0.211
#   Degrees of freedom                                 1
#   P-value (Chi-square)                           0.646
#
# Parameter Estimates:
#
#   Standard errors                             Standard
#   Information                                 Expected
#   Information saturated (h1) model          Structured
#
# Regressions:
#                  Estimate  Std.Err  z-value  P(>|z|)
# A0_PM2.5 ~
#   L0_male           0.040    0.006    6.337    0.000
#   L0_sc_nv          0.058    0.007    8.869    0.000
# M_diabetes ~
#   L0_male           0.053    0.009    5.836    0.000
#   L0_sc_nv          0.066    0.009    6.974    0.000
#   A0_PM2.5 (b.A)    0.127    0.014    8.795    0.000
#   L1                0.070    0.010    7.078    0.000
# Y_qol ~
#   L0_male          -0.719    0.202   -3.567    0.000
#   L0_sc_nv         -2.888    0.211  -13.678    0.000
#   A0_PM2.5 (c.A)   -3.965    0.321  -12.349    0.000
#   L1               -3.425    0.221  -15.488    0.000
#   M_diabts (c.M)   -8.714    0.222  -39.249    0.000
#
# Covariances:
#                  Estimate  Std.Err  z-value  P(>|z|)
# L0_male ~~
#   L0_soc_env       -0.003    0.002   -1.167    0.243
#   L1               -0.002    0.002   -0.832    0.406
# L0_soc_env ~~
#   L1               -0.001    0.002   -0.468    0.640
#
# Variances:
#                  Estimate  Std.Err  z-value  P(>|z|)
#  .A0_PM2.5          0.099    0.001   70.711    0.000
#  .M_diabetes        0.205    0.003   70.711    0.000
#  .Y_qol           100.882    1.427   70.711    0.000
#   L0_male           0.250    0.004   70.711    0.000
#   L0_soc_env        0.229    0.003   70.711    0.000
#   L1                0.207    0.003   70.711    0.000
#
# Defined Parameters:
#                Estimate  Std.Err  z-value  P(>|z|)
# direct           -3.965    0.321  -12.349    0.000
# indirect         -1.104    0.129   -8.582    0.000
# total            -5.069    0.344  -14.753    0.000

## using the Baron & Kenny approach, results were:
# - direct effect = -3.965038
# - indirect effect = -1.131019 (difference coef) or -1.108213 (product coef)
# - total effect = -5.096057
# these results are quite different from the Baron & Kenny approach !

fit.qol <- lavaan::sem(model = sem.QoL,
                       data = df1.qol.sem)
summary(fit.qol)
# lavaan 0.6-18 ended normally after 76 iterations
#
# Estimator                                       DWLS
# Optimization method                           NLMINB
# Number of model parameters                        24
#
# Number of observations                         10000
#
# Model Test User Model:
#                                             Standard      Scaled
# Test Statistic                                 0.209       0.214
# Degrees of freedom                                 1           1
# P-value (Chi-square)                           0.648       0.644
# Scaling correction factor                                  0.977
# Shift parameter                                           -0.000
# simple second-order correction
#
# Parameter Estimates:
#   Parameterization                               Delta
#   Standard errors                           Robust.sem
#   Information                                 Expected
#   Information saturated (h1) model        Unstructured
#
# Regressions:
#                  Estimate  Std.Err  z-value   P(>|z|)
# A0_PM2.5 ~
#   L0_mal  (a.01)    0.209    0.033     6.394    0.000
#   L0_sc_  (a.02)    0.324    0.035     9.129    0.000
# M_diabetes ~
#   L0_mal (b.L01)    0.129    0.027     4.846    0.000
#   L0_sc_ (b.L02)    0.156    0.029     5.463    0.000
#   A0_PM2   (b.A)    0.182    0.021     8.498    0.000
#   L1      (b.L1)    0.199    0.028     7.083    0.000
# Y_qol ~
#   L0_mal  (c.01)   -0.193    0.214    -0.899    0.369
#   L0_sc_  (c.02)   -2.130    0.232    -9.195    0.000
#   A0_PM2   (c.A)   -1.754    0.169   -10.351    0.000
#   L1      (c.L1)   -3.070    0.229   -13.411    0.000
#   M_dbts   (c.M)   -4.927    0.142   -34.802    0.000
#
# Covariances:
#                   Estimate  Std.Err  z-value   P(>|z|)
# L0_male ~~
#   L0_soc_env       -0.003    0.002    -1.167    0.243
#   L1               -0.002    0.002    -0.801    0.423
# L0_soc_env ~~
#   L1               -0.001    0.002    -0.425    0.671
#
# Intercepts:
#                  Estimate  Std.Err  z-value   P(>|z|)
#  .Y_qol            72.816    0.228   319.590    0.000
#   L0_male           0.501    0.005   100.197    0.000
#   L0_soc_env        0.644    0.003   187.422    0.000
#   L1                0.293    0.002   128.425    0.000
#
# Thresholds:                                           # Thresholds for the
#                Estimate  Std.Err  z-value   P(>|z|)   # latent response variables
# A0_PM2.5|t1       1.527    0.031    48.583    0.000   # (binary M and Y)
# M_diabetes|t1     0.810    0.026    30.789    0.000
#
# Variances:
#                 Estimate  Std.Err  z-value   P(>|z|)
# .A0_PM2.5          0.965
# .M_diabetes        0.975
# .Y_qol            88.521    1.538    57.551    0.000
#  L0_male           0.250    0.000 12499.500    0.000
#  L0_soc_env        0.229    0.002   106.346    0.000
#  L1                0.207    0.002    91.060    0.000
#
# Defined Parameters:
#                Estimate  Std.Err  z-value  P(>|z|)
# direct           -1.754    0.169   -10.351    0.000
# indirect         -0.896    0.099    -9.041    0.000
# total            -2.651    0.177   -14.942    0.000

## using the Baron & Kenny approach, results were:
# - direct effect = -3.965038
# - indirect effect = -1.131019 (difference coef) or -1.108213 (product coef)
# - total effect = -5.096057
# these results are quite different from the Baron & Kenny approach !

# If we keep the mediator as a continuous variable:
fit.qol <- lavaan::sem(model = sem.QoL,
                       data = df1, # original data.frame with continuous mediator
                       se = "bootstrap", # estimate SE and 95%CI by bootstrap
                       bootstrap = 100) # number of bootstrap sample (1000 is better)
summary(fit.qol)
# Defined Parameters:
#                Estimate  Std.Err  z-value  P(>|z|)
# direct           -3.965    0.375  -10.568    0.000
# indirect         -1.104    0.137   -8.057    0.000
# total            -5.069    0.400  -12.675    0.000
## With continuous mediator, we find the same results as with the Baron & Kenny
## approach

## Figure
semPaths(fit.qol,
         what = "est",
         layout = "tree2", # tree, tree2, spring
         rotation = 2, #exogenous on the left, endogenous on the right
         sizeMan = 10, # font size of manifest variable names
         nCharNodes = 0,
         nCharEdges = 0, # don't limit variable name lengths
         edge.label.cex = 0.6,
         curvePivot = TRUE,
         fade = FALSE)


## death (binary) ----
df1.death.sem <- data.frame(L0_male = df1$L0_male,
                            L0_soc_env = df1$L0_soc_env,
                            A0_PM2.5 = ordered(df1$A0_PM2.5), #df1$A0_PM2.5, #o
                            L1 = df1$L1,
                            M_diabetes = ordered(df1$M_diabetes), #df1$M_diabetes, #
                            Y_death = ordered(df1$Y_death)) #df1$Y_death) #

sem.death <- "
## SEM for binary outcome (QoL)
# Regression models
  A0_PM2.5 ~ a.01 * L0_male + a.02 * L0_soc_env
  M_diabetes ~ b.L01 * L0_male + b.L02 * L0_soc_env + b.A * A0_PM2.5 + b.L1 * L1
  Y_death ~ c.01 * L0_male + c.02 * L0_soc_env + c.A * A0_PM2.5 + c.L1 * L1 + c.M * M_diabetes

# Assuming the possibility of non-null covariance between confounders
# (note: in the data-generating system, the null assumptions were true)
  # L0_male ~~ L0_soc_env
  # L0_male ~~ L1
  # L0_soc_env ~~ L1

# define other parameters: direct, indirect and total effects
  direct := c.A
  indirect := b.A * c.M
  total := (b.A * c.M) + c.A
"

fit.death <- lavaan::sem(model = sem.death,
                       data = df1.death.sem)
summary(fit.death)
# lavaan 0.6-18 ended normally after 36 iterations
#
# Estimator                                       DWLS
# Optimization method                           NLMINB
# Number of model parameters                        14
#
# Number of observations                         10000
#
# Model Test User Model:
#   Standard      Scaled
# Test Statistic                                 0.268       0.268
# Degrees of freedom                                 1           1
# P-value (Chi-square)                           0.605       0.605
# Scaling correction factor                                  1.000
# Shift parameter                                           -0.000
# simple second-order correction
#
# Parameter Estimates:
#
#   Parameterization                               Delta
#   Standard errors                           Robust.sem
#   Information                                 Expected
#   Information saturated (h1) model        Unstructured
#
# Regressions:
#                  Estimate  Std.Err  z-value  P(>|z|)
# A0_PM2.5 ~
#   L0_mal  (a.01)    0.214    0.034    6.365    0.000
#   L0_sc_  (a.02)    0.331    0.037    8.916    0.000
# M_diabetes ~
#   L0_mal (b.L01)    0.131    0.027    4.833    0.000
#   L0_sc_ (b.L02)    0.158    0.029    5.415    0.000
#   A0_PM2   (b.A)    0.183    0.021    8.708    0.000
#   L1      (b.L1)    0.203    0.029    7.078    0.000
# Y_death ~
#   L0_mal  (c.01)    0.152    0.029    5.226    0.000
#   L0_sc_  (c.02)    0.181    0.032    5.687    0.000
#   A0_PM2   (c.A)    0.080    0.024    3.363    0.001
#   L1      (c.L1)    0.249    0.031    8.117    0.000
#   M_dbts   (c.M)    0.125    0.019    6.434    0.000
#
# Thresholds:
#                Estimate  Std.Err  z-value  P(>|z|)
# A0_PM2.5|t1       1.562    0.039   39.854    0.000
# M_diabetes|t1     0.822    0.029   28.586    0.000
# Y_death|t1        1.164    0.032   36.566    0.000
#
# Variances:
#                 Estimate  Std.Err  z-value  P(>|z|)
# .A0_PM2.5          1.000
# .M_diabetes        0.967
# .Y_death           0.974
#
# Defined Parameters:
#                Estimate  Std.Err  z-value  P(>|z|)
# direct            0.080    0.024    3.363    0.001
# indirect          0.023    0.004    5.419    0.000
# total             0.103    0.023    4.465    0.000

## using the Baron & Kenny approach, results were:
# - direct effect = 0.05150901
# - indirect effect = 0.008737889 (difference coef) or 0.008234978 (product coef)
# - total effect = 0.060247
# Results from sem are quite different from the Baron & Kenny approach !

# If we keep the mediator and the outcome as a continuous variable:
fit.death <- sem(sem.death,
                 data = df1, # original data.frame with continuous mediator
                 se = "bootstrap", # estimate SE and 95%CI by bootstrap
                 bootstrap = 100) # number of bootstrap sample (1000 is better)
summary(fit.death)
# Defined Parameters:
#                Estimate  Std.Err  z-value  P(>|z|)
# direct            0.052    0.012    4.140    0.000
# indirect          0.008    0.002    4.903    0.000
# total             0.060    0.012    4.912    0.000

# Same results as with Baron & Kenny

# conclusion: for binary variables => it might be better not to declare them as
#                                     ordered variables

# for multicategorical mediator: define dummy variables, and use them as continuous
# variable, and add non-null covariances between the dummy variables


## quality of life with A*M interaction----
df1.qol.int <- data.frame(L0_male = df1_int$L0_male,
                          L0_soc_env = df1_int$L0_soc_env,
                          A0_PM2.5 = df1_int$A0_PM2.5, # keep exposure as numeric
                          L1 = df1_int$L1,
                          M_diabetes = df1_int$M_diabetes, # keep mediator as numeric
                          AM_inter = df1_int$A0_PM2.5 * df1_int$M_diabetes,
                          Y_qol = df1_int$Y_qol)
sem.QoL.int <- "
## SEM for quantitative outcome (QoL)
# Regression models
  A0_PM2.5 ~ a.01 * L0_male + a.02 * L0_soc_env
  M_diabetes ~ b.L01 * L0_male + b.L02 * L0_soc_env + b.A * A0_PM2.5 + b.L1 * L1
  AM_inter ~ d.L01 * L0_male + d.L02 * L0_soc_env + d.A * A0_PM2.5 + d.L1 * L1
  Y_qol ~ c.01*L0_male+c.02*L0_soc_env+c.A*A0_PM2.5+c.L1*L1+c.M*M_diabetes+c.AM*AM_inter

# Assuming the possibility of non-null covariance between confounders
# (note: in the data-generating system, the null assumptions were true)
  L0_male ~~ L0_soc_env
  L0_male ~~ L1
  L0_soc_env ~~ L1

# define other parameters: direct, indirect and total effects
  direct := c.A
  indirect := b.A * c.M
  mie := d.A * c.AM
  total := (b.A * c.M) + c.A + (d.A * c.AM)
"

fit.qol <- lavaan::sem(model = sem.QoL.int,
                       data = df1.qol.int,
                       se = "bootstrap",
                       bootstrap = 100)
semPaths(fit.qol,
         what = "est",
         layout = "tree2", # tree, tree2, spring
         rotation = 2, #exogenous on the left, endogenous on the right
         sizeMan = 10, # font size of manifest variable names
         nCharNodes = 0,
         nCharEdges = 0, # don't limit variable name lengths
         edge.label.cex = 0.6,
         curvePivot = TRUE,
         fade = FALSE)
summary(fit.qol)
# lavaan 0.6-18 ended normally after 27 iterations
#
# Estimator                                         ML
# Optimization method                           NLMINB
# Number of model parameters                        26
#
# Number of observations                         10000
#
# Model Test User Model:
#
#   Test statistic                              1423.845
#   Degrees of freedom                                 2
#   P-value (Chi-square)                           0.000
#
# Parameter Estimates:
#
#   Standard errors                            Bootstrap
#   Number of requested bootstrap draws              100
#   Number of successful bootstrap draws             100
#
# Regressions:
#                  Estimate  Std.Err  z-value  P(>|z|)
# A0_PM2.5 ~
#   L0_mal  (a.01)    0.040    0.005    7.480    0.000
#   L0_sc_  (a.02)    0.058    0.007    8.625    0.000
#   M_diabetes ~
#   L0_mal (b.L01)    0.053    0.008    6.741    0.000
#   L0_sc_ (b.L02)    0.066    0.009    7.364    0.000
#   A0_PM2   (b.A)    0.127    0.017    7.628    0.000
#   L1      (b.L1)    0.070    0.010    6.834    0.000
# AM_inter ~
#   L0_mal (d.L01)    0.000    0.003    0.048    0.962
#   L0_sc_ (d.L02)    0.005    0.003    1.599    0.110
#   A0_PM2   (d.A)    0.423    0.017   25.407    0.000
#   L1      (d.L1)    0.005    0.004    1.325    0.185
# Y_qol ~
#   L0_mal  (c.01)   -0.724    0.210   -3.448    0.001
#   L0_sc_  (c.02)   -2.890    0.197  -14.667    0.000
#   A0_PM2   (c.A)   -3.715    0.435   -8.540    0.000
#   L1      (c.L1)   -3.428    0.274  -12.533    0.000
#   M_dbts   (c.M)   -8.632    0.239  -36.107    0.000
#   AM_ntr  (c.AM)   -5.615    0.652   -8.608    0.000
#
# Covariances:
#                  Estimate  Std.Err  z-value  P(>|z|)
# L0_male ~~
#   L0_soc_env       -0.003    0.002   -1.133    0.257
#   L1               -0.002    0.003   -0.743    0.458
# L0_soc_env ~~
#   L1               -0.001    0.002   -0.509    0.611
#
# Variances:
#                  Estimate  Std.Err  z-value  P(>|z|)
#  .A0_PM2.5          0.099    0.002   48.024    0.000
#  .M_diabetes        0.205    0.002  108.839    0.000
#  .AM_inter          0.027    0.001   35.725    0.000
#  .Y_qol           100.873    1.560   64.667    0.000
#   L0_male           0.250    0.000 7276.584    0.000
#   L0_soc_env        0.229    0.001  180.231    0.000
#   L1                0.207    0.002  101.196    0.000
#
# Defined Parameters:
#                Estimate  Std.Err  z-value  P(>|z|)
# direct           -3.715    0.437   -8.497    0.000
# indirect         -1.094    0.148   -7.369    0.000
# mie              -2.374    0.281   -8.446    0.000
# total            -7.183    0.383  -18.759    0.000


# ---------------------------------------------------------------------------- #
### test on simulations for interactions ----
# ---------------------------------------------------------------------------- #
# functions to generate data
param.causal.model.1 <- function(A.M.interaction = NULL) {
  # L0
  p_L0_male <- 0.5
  p_L0_soc_env <- 0.65

  # A: A0_PM2.5 <- rbinom( 0.05 + 0.04 * L0_male + 0.06 * L0_soc_env )
  b_A <- 0.05   # reference prevalence is 5%
  b_male_A <- 0.04  # + 0.04 for the effect of L0_male -> A0_PM2.5
  b_soc_env_A <- 0.06  # +0.06 for the effect of L0_soc_env -> A0_PM2.5

  # L1: intermediate confounder between M and Y, not influenced by A
  p_L1 <- 0.3

  # M: M_diabetes <- rbinom( 0.2 + 0.05 * L0_male + 0.06 * L0_soc_env + 0.07 * L1 +
  #                         0.1 * A0_PM2.5 )
  b_M <- 0.2 # reference prevalence is 20%
  b_male_M <- 0.05 # +0.05 for the effect of L0_male -> M_diabetes
  b_soc_env_M <- 0.06 # +0.06 for the effect of L0_soc_env -> M_diabetes
  b_L1_M <- 0.07 # +0.07 for the effect of L1 -> M_diabetes
  b_A_M <- 0.1 # +0.10 for the effect of A0_PM2.5 -> M_diabetes

  # Y binary: rbinom( 0.10 + 0.06 * L0_male + 0.04 * L0_soc_env + 0.05 * A0_PM2.5 +
  #                   0.07 * L1 + 0.08 * M_diabetes +
  #                   0.03 * A0_PM2.5 * M_diabetes * A.M.inter )
  b_Y <- 0.1 # reference prevalence is 10%
  b_male_Y <- 0.06 # +0.06 for the effect of L0_male -> Y
  b_soc_env_Y <- 0.04 # +0.04 for the effect of L0_soc_env -> Y
  b_A_Y <- 0.05 # 0.05 for the effect of A0_PM2.5 -> Y
  b_L1_Y <- 0.07 # +0.07 for the effect of L1 -> Y
  b_M_Y <- 0.08 # 0.08 for the effect of M_diabetes -> Y
  b_AM_Y <- 0.03 # 0.03 for the interaction effect A0_PM2.5 * M_diabetes -> Y

  # Y continuous: (75 - 1 * L0_male - 3 * L0_soc_env - 4 * A0_PM2.5 -3.5 * L1 -
  #                9 * M_diabetes -5 * A0_PM2.5 * M_diabetes * A.M.inter ) +
  #                rnorm(N, mean = 0, sd = 10)
  mu_Y <- 75 # reference mean for QoL
  c_male_Y <- -1 # -1 for the effect of L0_male -> Y
  c_soc_env_Y <- -3 # -3 for the effect of L0_soc_env -> Y
  c_A_Y <- -4 # -4 for the effect of A0_PM2.5 -> Y
  c_L1_Y <- -3.5 # -3.5 for the effect of L1 -> Y
  c_M_Y <- -9 # -9 for the effect of M_diabetes -> Y
  c_AM_Y <- -5  # - 5 for the interaction effect A0_PM2.5 * M_diabetes  -> Y
  sd_Y <- 10 # standard deviation of the residuals

  # A*M interaction ?
  A.M.inter <- A.M.interaction

  coef <- c( p_L0_male = p_L0_male, p_L0_soc_env = p_L0_soc_env,
             b_A = b_A, b_male_A = b_male_A, b_soc_env_A = b_soc_env_A,
             p_L1 = p_L1,
             b_M = b_M, b_male_M = b_male_M, b_soc_env_M = b_soc_env_M,
             b_L1_M = b_L1_M, b_A_M = b_A_M,
             b_Y = b_Y, b_male_Y = b_male_Y, b_soc_env_Y = b_soc_env_Y,
             b_A_Y = b_A_Y, b_L1_Y = b_L1_Y, b_M_Y = b_M_Y, b_AM_Y = b_AM_Y,
             mu_Y = mu_Y, c_male_Y = c_male_Y, c_soc_env_Y = c_soc_env_Y,
             c_A_Y = c_A_Y, c_L1_Y = c_L1_Y, c_M_Y = c_M_Y, c_AM_Y = c_AM_Y,
             sd_Y = sd_Y, A.M.inter = A.M.inter)

  return(coef)
}

gen.data.causal.model.1 <- function(N, A.M.inter) { # input parameters are the
  #   sample size N and the presence of A*M interaction with A.M.inter = 0 or 1

  b <- param.causal.model.1(A.M.interaction = A.M.inter)

  # baseline confounders: parent's educational level=L0_soc_env & sex=L0_male
  L0_male <- rbinom(N, size = 1, prob = b["p_L0_male"])
  L0_soc_env <- rbinom(N, size = 1, prob = b["p_L0_soc_env"])

  # exposure: A0_PM2.5
  A0_PM2.5 <- rbinom(N, size = 1, prob =  b["b_A"] +
                       b["b_male_A"] * L0_male +
                       b["b_soc_env_A"] * L0_soc_env )

  # intermediate confounder between M_diabetes and Y, not affected by A0 L1
  L1 <- rbinom(N, size = 1, prob = b["p_L1"])

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
results.int <- matrix(NA, nrow = 1000, ncol = 12,
                      dimnames = list(c(1:1000),
                                      c("direct.1","indirect.1","mie.1","total.1",
                                        "direct.2","indirect.2","mie.2","total.2",
                                        "direct.3","indirect.3","mie.3","total.3")))

set.seed(54321)
for(i in 1:1000) {
  print(paste0("simulation ",i))
  # generate data
  df1_int <- gen.data.causal.model.1(N=10000, A.M.inter=1)

  # create data.frame with interaction term
  df1.qol.int <- data.frame(L0_male = df1_int$L0_male,
                            L0_soc_env = df1_int$L0_soc_env,
                            A0_PM2.5 = df1_int$A0_PM2.5, # keep exposure as numeric
                            L1 = df1_int$L1,
                            M_diabetes = df1_int$M_diabetes, # keep mediator as numeric
                            AM_inter = df1_int$A0_PM2.5 * df1$M_diabetes,
                            Y_qol = df1_int$Y_qol)
  # define model 1 (arrow from A to A*M)
  sem.QoL.int.1 <- "
  ## SEM for quantitative outcome (QoL)
  # Regression models
    A0_PM2.5 ~ a.01 * L0_male + a.02 * L0_soc_env
    M_diabetes ~ b.L01 * L0_male + b.L02 * L0_soc_env + b.A * A0_PM2.5 + b.L1 * L1
    AM_inter ~ d.L01 * L0_male + d.L02 * L0_soc_env + d.A * A0_PM2.5 + d.L1 * L1
    Y_qol ~ c.01*L0_male+c.02*L0_soc_env+c.A*A0_PM2.5+c.L1*L1+c.M*M_diabetes+c.AM*AM_inter

  # Assuming the possibility of non-null covariance between confounders
  # (note: in the data-generating system, the null assumptions were true)
    L0_male ~~ L0_soc_env
    L0_male ~~ L1
    L0_soc_env ~~ L1

  # define other parameters: direct, indirect and total effects
    direct := c.A
    indirect := b.A * c.M
    mie := d.A * c.AM
    total := (b.A * c.M) + c.A + (d.A * c.AM)
  "

  # define model 2 (no arrow from A to A*M)
  sem.QoL.int.2 <- "
  ## SEM for quantitative outcome (QoL)
  # Regression models
    A0_PM2.5 ~ a.01 * L0_male + a.02 * L0_soc_env
    M_diabetes ~ b.L01 * L0_male + b.L02 * L0_soc_env + b.A * A0_PM2.5 + b.L1 * L1
    Y_qol ~ c.01*L0_male+c.02*L0_soc_env+c.A*A0_PM2.5+c.L1*L1+c.M*M_diabetes+c.AM*AM_inter

  # Assuming the possibility of non-null covariance between confounders
  # (note: in the data-generating system, the null assumptions were true)
    L0_male ~~ L0_soc_env
    L0_male ~~ L1
    L0_soc_env ~~ L1

  # define other parameters: direct, indirect and total effects
    direct := c.A
    indirect := b.A * c.M
    mie := b.A * c.AM
    total := (b.A * c.M) + c.A + (b.A * c.AM)
  "

  fit.qol.1 <- lavaan::sem(model = sem.QoL.int.1,
                           data = df1.qol.int)
  res1 <- summary(fit.qol.1)

  fit.qol.2 <- lavaan::sem(model = sem.QoL.int.2,
                           data = df1.qol.int)
  res2 <- summary(fit.qol.2)

  # traditional regressions
  trad_ATE_qol <- lm(Y_qol ~ A0_PM2.5 + L0_male + L0_soc_env,
                     data = df1_int)

  trad_qol_am <- lm(Y_qol ~ A0_PM2.5 + M_diabetes + A0_PM2.5:M_diabetes +
                      L0_male + L0_soc_env + L1,
                    data = df1_int)
  gamma.A.q <- coef(trad_qol_am)["A0_PM2.5"]
  gamma.M.q <- coef(trad_qol_am)["M_diabetes"]
  gamma.AM.q <- coef(trad_qol_am)["A0_PM2.5:M_diabetes"]

  trad_m <- glm(M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env, # + L1, # il semble qu'il faut aussi ajuster sur L1 ici ?
                family = "binomial",
                data = df1_int)
  beta.0 <- coef(trad_m)["(Intercept)"]
  beta.A <- coef(trad_m)["A0_PM2.5"]

  ## The PNDE is:
  trad_PNDE_qol <- gamma.A.q + gamma.AM.q * (exp(beta.0)) / (1 + exp(beta.0))
  ## The PNIE is:
  trad_PNIE_qol <- gamma.M.q * (exp(beta.0 + beta.A) / (1 + exp(beta.0 + beta.A)) - exp(beta.0) / (1 + exp(beta.0)))
  ## The MIE is:
  trad_MIE_qol <- gamma.AM.q * (exp(beta.0 + beta.A) / (1 + exp(beta.0 + beta.A)) - exp(beta.0) / (1 + exp(beta.0)))

  # save results
  results.int[i,"direct.1"] <- res1$pe[res1$pe$lhs == "direct","est"] - (-5.425) # direct -5.425
  results.int[i,"indirect.1"] <- res1$pe[res1$pe$lhs == "indirect","est"] - (-0.9) # indirect -0.9
  results.int[i,"mie.1"] <- res1$pe[res1$pe$lhs == "mie","est"] - (-0.5) # mie -0.5
  results.int[i,"total.1"] <- res1$pe[res1$pe$lhs == "total","est"] - (-6.825) # ATE -6.825
  results.int[i,"direct.2"] <- res2$pe[res2$pe$lhs == "direct","est"] - (-5.425)
  results.int[i,"indirect.2"] <- res2$pe[res2$pe$lhs == "indirect","est"] - (-0.9)
  results.int[i,"mie.2"] <- res2$pe[res2$pe$lhs == "mie","est"] - (-0.5)
  results.int[i,"total.2"] <- res2$pe[res2$pe$lhs == "total","est"] - (-6.825)
  results.int[i,"direct.3"] <- trad_PNDE_qol - (-5.425)
  results.int[i,"indirect.3"] <- trad_PNIE_qol - (-0.9)
  results.int[i,"mie.3"] <- trad_MIE_qol - (-0.5)
  results.int[i,"total.3"] <- coef(trad_ATE_qol)["A0_PM2.5"] - (-6.825)
}

sapply(data.frame(results.int) , mean)
## résultat en ajustant sur L1 dans le modèle du médiateur
#    direct.1  indirect.1       mie.1     total.1    direct.2  indirect.2       mie.2     total.2    direct.3  indirect.3       mie.3
# -0.49736317 -0.05784399  0.49819219 -0.05701497 -0.49736317 -0.05784399  0.49909038 -0.05611678  0.37653796  0.15473379  0.08993298
#     total.3
# -0.05628162

## résultat sans ajuster sur L1 dans le modèle du médiateur
#    direct.1  indirect.1       mie.1     total.1    direct.2  indirect.2       mie.2     total.2    direct.3  indirect.3       mie.3
# -0.49736317 -0.05784399  0.49819219 -0.05701497 -0.49736317 -0.05784399  0.49909038 -0.05611678  0.28705535  0.12002520  0.07081778
#     total.3
# -0.05628162

boxplot(data.frame(results.int))

boxplot(subset(data.frame(results.int),
               select = c("direct.1","direct.2","direct.3"))) # modèles 1 et 2 similaires ; modèle 3 un peu moins biaisé
boxplot(subset(data.frame(results.int),
               select = c("indirect.1","indirect.2","indirect.3"))) # modèles 1 et 2 similaires ; modèle 3 un peu plus biaisé
boxplot(subset(data.frame(results.int),
               select = c("mie.1","mie.2","mie.3"))) # modèles 2 plus précis que le 3 ; modèle 3 moins biaisé
boxplot(subset(data.frame(results.int),
               select = c("total.1","total.2","total.3"))) # 3 modèles très proches, le modèle 3 plus précis

results.int.rel <- results.int
results.int.rel[,c(1,5,9)] <- results.int[,c(1,5,9)] / (-5.425)
results.int.rel[,c(2,6,10)] <- results.int[,c(2,6,10)] / (-0.9)
results.int.rel[,c(3,7,11)] <- results.int[,c(3,7,11)] / (-0.5)
results.int.rel[,c(4,8,12)] <- results.int[,c(4,8,12)] / (-6.825)

sapply(data.frame(results.int.rel) , mean)
## résultats en ajustant sur L1
#    direct.1   indirect.1        mie.1      total.1     direct.2   indirect.2        mie.2      total.2     direct.3   indirect.3
# 0.091679847  0.064271096 -0.996384375  0.008353842  0.091679847  0.064271096 -0.998180759  0.008222239 -0.069407919 -0.171926433
#        mie.3      total.3
# -0.179865966  0.008246392

## résultats sans ajuster sur L1
#    direct.1   indirect.1        mie.1      total.1     direct.2   indirect.2        mie.2      total.2     direct.3   indirect.3
# 0.091679847  0.064271096 -0.996384375  0.008353842  0.091679847  0.064271096 -0.998180759  0.008222239 -0.052913428 -0.133361335
#        mie.3      total.3
# -0.141635557  0.008246392
# le modèles 1 et 2 font de mauvaises estimations de l'interaction médiée
# le modèle 3 semble globalement meilleur, mais biais > 10% pour l'effet indirect et l'interaction médiée





## quality of life confusion intermédiaire ----
df2 <- read.csv(file = "data/df2.csv")
names(df2)
# [1] "L0_male"    "L0_soc_env" "A0_PM2.5"   "L1"         "M_diabetes" "Y_death"    "Y_qol"

## quality of life ----
sem.QoL.df2 <- "
## SEM for quantitative outcome (QoL)
# Regression models
  A0_PM2.5 ~ a.L01 * L0_male + a.L02 * L0_soc_env
  L1 ~ b.L01 * L0_male + b.L02 * L0_soc_env + b.A * A0_PM2.5
  M_diabetes ~ c.L01 * L0_male + c.L02 * L0_soc_env + c.A * A0_PM2.5 + c.L1 * L1
  Y_qol ~ d.01 * L0_male + d.02 * L0_soc_env + d.A * A0_PM2.5 + d.L1 * L1 + d.M * M_diabetes

# Assuming the possibility of non-null covariance between confounders
# (note: in the data-generating system, the null assumptions were true)
  L0_male ~~ L0_soc_env

# define other parameters: specific paths, direct, indirect and total effects
  path.A_Y := d.A
  path.A_L1_Y := b.A * d.L1
  path.A_M_Y := c.A * d.M
  path.A_L1_M_Y := b.A * c.L1 * d.M
  marg.direct := d.A + (b.A * d.L1)
  marg.indirect := (c.A * d.M) + (b.A * c.L1 * d.M)
  total := d.A + (b.A * d.L1) + (c.A * d.M) + (b.A * c.L1 * d.M)
  cond.direct := d.A + (b.A * d.L1) + (b.A * c.L1 * d.M)
  cond.indirect := (c.A * d.M)
"

fit.qol.df2 <- lavaan::sem(model = sem.QoL.df2,
                           data = df2)
summary(fit.qol.df2)
# lavaan 0.6-18 ended normally after 18 iterations
#
# Estimator                                         ML
# Optimization method                           NLMINB
# Number of model parameters                        21
#
# Number of observations                         10000
#
# Model Test User Model:
#
#   Test statistic                                 0.000
#   Degrees of freedom                                 0
#
# Parameter Estimates:
#
#   Standard errors                             Standard
#   Information                                 Expected
#   Information saturated (h1) model          Structured
#
# Regressions:
#                  Estimate  Std.Err  z-value  P(>|z|)
# A0_PM2.5 ~
#   L0_mal (a.L01)    0.040    0.006    6.337    0.000
#   L0_sc_ (a.L02)    0.058    0.007    8.869    0.000
# L1 ~
#   L0_mal (b.L01)   -0.043    0.009   -4.633    0.000
#   L0_sc_ (b.L02)    0.069    0.010    7.055    0.000
#   A0_PM2   (b.A)    0.226    0.015   15.170    0.000
#   M_diabetes ~
#   L0_mal (c.L01)    0.052    0.009    5.658    0.000
#   L0_sc_ (c.L02)    0.064    0.010    6.575    0.000
#   A0_PM2   (c.A)    0.072    0.015    4.846    0.000
#   L1      (c.L1)    0.194    0.010   19.779    0.000
# Y_qol ~
#   L0_mal  (d.01)   -0.725    0.202   -3.594    0.000
#   L0_sc_  (d.02)   -2.881    0.212  -13.617    0.000
#   A0_PM2   (d.A)   -3.926    0.324  -12.120    0.000
#   L1      (d.L1)   -5.165    0.219  -23.611    0.000
#   M_dbts   (d.M)   -8.698    0.218  -39.849    0.000
#
# Covariances:
#                  Estimate  Std.Err  z-value  P(>|z|)
# L0_male ~~
#   L0_soc_env       -0.003    0.002   -1.167    0.243
#
# Variances:
#                  Estimate  Std.Err  z-value  P(>|z|)
#  .A0_PM2.5          0.099    0.001   70.711    0.000
#  .L1                0.219    0.003   70.711    0.000
#  .M_diabetes        0.212    0.003   70.711    0.000
#  .Y_qol           100.878    1.427   70.711    0.000
#   L0_male           0.250    0.004   70.711    0.000
#   L0_soc_env        0.229    0.003   70.711    0.000
#
# Defined Parameters:
#                Estimate  Std.Err  z-value  P(>|z|)
# path.A_Y         -3.926    0.324  -12.120    0.000
# path.A_L1_Y      -1.168    0.092  -12.762    0.000
# path.A_M_Y       -0.625    0.130   -4.810    0.000    true : -0.9
# path.A_L1_M_Y    -0.382    0.033  -11.523    0.000
# marg.direct      -5.093    0.330  -15.437    0.000    true : -5
# marg.indirect    -1.007    0.132   -7.609    0.000    true : -1.26
# total            -6.101    0.359  -16.986    0.000    true : -6.26
# cond.direct      -5.476    0.337  -16.263    0.000    true : -5.36
# cond.indirect    -0.625    0.130   -4.810    0.000    true : -0.9

## comparing these estimations with results from the CMAverse package
library(CMAverse)
cmdag(outcome = "Y_qol", exposure = "A0_PM2.5", mediator = "M_diabetes",
      basec = c("L0_male", "L0_soc_env"), postc = "L1", node = TRUE, text_col = "white")

set.seed(1234)
res_gformula_Qol <- cmest(data = df2,
                          model = "gformula", # for parametric g-computation
                          outcome = "Y_qol", # outcome variable
                          exposure = "A0_PM2.5", # exposure variable
                          mediator = "M_diabetes", # mediator
                          basec = c("L0_male",
                                    "L0_soc_env"), # confounders
                          postc = "L1", # intermediate confounder (post-exposure)
                          EMint = FALSE, # exposures*mediator interaction
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
res.cmavers <- summary(res_gformula_Qol)
res.cmavers$effect.pe


# ---------------------------------------------------------------------------- #
### test on simulations for df2 ----
# ---------------------------------------------------------------------------- #
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
results.df2 <- matrix(NA, nrow = 1000, ncol = 6,
                      dimnames = list(c(1:1000),
                                      c("direct.1","indirect.1","total.1",
                                        "direct.2","indirect.2","total.2")))

set.seed(54321)
for(i in 1:1000) {
  print(paste0("simulation ",i))
  # generate data
  df2 <- gen.data.causal.model.2(N=10000, A.M.inter=0)

  # create data.frame with interaction term
  df1.qol.int <- data.frame(L0_male = df1_int$L0_male,
                            L0_soc_env = df1_int$L0_soc_env,
                            A0_PM2.5 = df1_int$A0_PM2.5, # keep exposure as numeric
                            L1 = df1_int$L1,
                            M_diabetes = df1_int$M_diabetes, # keep mediator as numeric
                            AM_inter = df1_int$A0_PM2.5 * df1$M_diabetes,
                            Y_qol = df1_int$Y_qol)
  ## quality of life ----
  sem.QoL.df2 <- "
  ## SEM for quantitative outcome (QoL)
  # Regression models
    A0_PM2.5 ~ a.L01 * L0_male + a.L02 * L0_soc_env
    L1 ~ b.L01 * L0_male + b.L02 * L0_soc_env + b.A * A0_PM2.5
    M_diabetes ~ c.L01 * L0_male + c.L02 * L0_soc_env + c.A * A0_PM2.5 + c.L1 * L1
    Y_qol ~ d.01 * L0_male + d.02 * L0_soc_env + d.A * A0_PM2.5 + d.L1 * L1 + d.M * M_diabetes

  # Assuming the possibility of non-null covariance between confounders
  # (note: in the data-generating system, the null assumptions were true)
    L0_male ~~ L0_soc_env

  # define other parameters: specific paths, direct, indirect and total effects
    path.A_Y := d.A
    path.A_L1_Y := b.A * d.L1
    path.A_M_Y := c.A * d.M
    path.A_L1_M_Y := b.A * c.L1 * d.M
    marg.direct := d.A + (b.A * d.L1)
    marg.indirect := (c.A * d.M) + (b.A * c.L1 * d.M)
    total := d.A + (b.A * d.L1) + (c.A * d.M) + (b.A * c.L1 * d.M)
    cond.direct := d.A + (b.A * d.L1) + (b.A * c.L1 * d.M)
    cond.indirect := (c.A * d.M)
  "

  fit.qol.df2 <- lavaan::sem(model = sem.QoL.df2,
                             data = df2)
  res.lavaan <- summary(fit.qol.df2)


  # CMAverse
  res_gformula_Qol <- cmest(data = df2,
                            model = "gformula", # for parametric g-computation
                            outcome = "Y_qol", # outcome variable
                            exposure = "A0_PM2.5", # exposure variable
                            mediator = "M_diabetes", # mediator
                            basec = c("L0_male",
                                      "L0_soc_env"), # confounders
                            postc = "L1", # intermediate confounder (post-exposure)
                            EMint = FALSE, # exposures*mediator interaction
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
  res.cmavers <- summary(res_gformula_Qol)
  res.cmavers$effect.pe

  # save results
  results.df2[i,"direct.1"] <- res.lavaan$pe[res.lavaan$pe$lhs == "marg.direct","est"] - (-5) # direct -5
  results.df2[i,"indirect.1"] <- res.lavaan$pe[res.lavaan$pe$lhs == "marg.indirect","est"] - (-1.26) # indirect -1.26
  results.df2[i,"total.1"] <- res.lavaan$pe[res.lavaan$pe$lhs == "total","est"] - (-6.26) # ATE -6.26
  results.df2[i,"direct.2"] <- res.cmavers$effect.pe["rpnde"] - (-5) # direct -5
  results.df2[i,"indirect.2"] <- res.cmavers$effect.pe["rpnie"] - (-1.26) # indirect -1.26
  results.df2[i,"total.2"] <- res.cmavers$effect.pe["te"] - (-6.26) # ATE -6.26
}

sapply(data.frame(results.df2) , mean)
#    direct.1  indirect.1     total.1    direct.2  indirect.2     total.2
# 0.003239160 0.008998731 0.012237891 0.007879356 0.026796054 0.034675410

boxplot(data.frame(results.df2))
# globalement les performances ont l'air très proches !!
boxplot(subset(data.frame(results.df2), select = c("direct.1", "direct.2")))
boxplot(subset(data.frame(results.df2), select = c("indirect.1", "indirect.2")))
boxplot(subset(data.frame(results.df2), select = c("total.1", "total.2")))

results.df2.rel <- results.df2
results.df2.rel[,c("direct.1","direct.2")] <- results.df2.rel[,c("direct.1","direct.2")] / (-5)
results.df2.rel[,c("indirect.1","indirect.2")] <- results.df2.rel[,c("indirect.1","indirect.2")] / (-1.26)
results.df2.rel[,c("total.1","total.2")] <- results.df2.rel[,c("total.1","total.2")] / (-6.26)
sapply(data.frame(results.df2.rel) , mean)
#     direct.1   indirect.1      total.1     direct.2   indirect.2      total.2
# -0.000647832 -0.007141850 -0.001954935 -0.001575871 -0.021266710 -0.005539203


