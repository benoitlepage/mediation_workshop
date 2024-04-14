###

rm(list=ls())

df1 <- read.csv(file = "data/df1.csv")


### Baron & Kenny

## Model 1 to estimate the total effect:
model.tot.A.QoL <- lm(Y_qol ~ A0_ace + L0_male + L0_parent_low_educ_lv,
                      data = df1)
summary(model.tot.A.QoL)
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)
#   (Intercept)            71.8820     0.2155 333.565  < 2e-16 ***
#   A0_ace                 -5.0961     0.3486 -14.617  < 2e-16 *** <- Total effect
#   L0_male                -1.1486     0.2194  -5.235 1.68e-07 ***
#   L0_parent_low_educ_lv  -3.4441     0.2295 -15.005  < 2e-16 ***

# The total effect of ACE on Quality of life is approximately equal to -5.1,
# given by the A0_ace coefficient:
model.tot.A.QoL$coefficients["A0_ace"]
# -5.096057

## Model 2 to estimate the effect of the exposure on the mediator
## because the mediator is binary, we might want to use a logistic or probit regression
## for example
logit.model.A.M <- glm(M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv,
                       data = df1, family = "binomial")
summary(logit.model.A.M)
# effects estimated on the logit scale:
# Coefficients:
#                         Estimate Std. Error z value Pr(>|z|)
#   (Intercept)           -1.27152    0.04584 -27.736  < 2e-16 ***
#   A0_ace                 0.56168    0.06537   8.592  < 2e-16 *** <- effect of A on M
#   L0_male                0.25455    0.04425   5.753 8.77e-09 ***
#   L0_parent_low_educ_lv  0.32683    0.04731   6.908 4.91e-12 ***
exp(coefficients(logit.model.A.M)["A0_ace"])
# Odds ratio = 1.753609  for the effect of ACE on the mediator (probability of smoking)

## Model 3 to estimate the direct effect of the exposure (conditional on the outcome) and
## the effect of M on Y, adjusted for confounders of the A-Y and M-Y relationships
model.A.M.QoL <- lm(Y_qol ~ A0_ace + M_smoking + L1 + L0_male + L0_parent_low_educ_lv,
                    data = df1)
summary(model.A.M.QoL)
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)
#   (Intercept)            74.7858     0.2130 351.178  < 2e-16 ***
#   A0_ace                 -3.9650     0.3212 -12.345  < 2e-16 *** <- Direct effect
#   M_smoking              -8.7138     0.2221 -39.237  < 2e-16 *** <- effect of M on Y
#   L1                     -3.4252     0.2212 -15.483  < 2e-16 ***
#   L0_male                -0.7193     0.2017  -3.566 0.000364 ***
#   L0_parent_low_educ_lv  -2.8876     0.2112 -13.674  < 2e-16 ***

# The direct effect is approximately -4.0 given by the A0_ace coefficient
model.A.M.QoL$coefficients["A0_ace"]
# -3.965038

## Following the Baron & Kenny Steps, we would conclude that :
# - There is a significant total effect of ACEs on Quality of Life (Model 1)
# - There is a significant effect of ACEs on the mediator (smoking) (Model 2)
# - There is a significant effect of the mediator (smoking) on Qol (model 3)
# - The direct effect is significantly non-null
# => Conclusion: Smoking partially mediates the relationship between ACEs and QoL

### Estimation of the indirect effect:
### We can apply the difference in coefficient method to estimate the indirect effect:
### substract the direct effect from the Total effect:
ind.effect.dif.meth <- (model.tot.A.QoL$coefficients["A0_ace"] -
                          model.A.M.QoL$coefficients["A0_ace"])

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
linear.model.A.M <- lm(M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv,
                       data = df1)
summary(linear.model.A.M)
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)
#   (Intercept)           0.215398   0.008929  24.123  < 2e-16 ***
#   A0_ace                0.127180   0.014446   8.804  < 2e-16 ***
#   L0_male               0.052363   0.009091   5.760 8.67e-09 ***
#   L0_parent_low_educ_lv 0.065806   0.009511   6.919 4.83e-12 ***

## product of coefficient method:
ind.effect.prod.meth <- (linear.model.A.M$coefficients["A0_ace"] *
                           model.A.M.QoL$coefficients["M_smoking"])
# -1.108213
# which also gives an indirect effect of approximately -1.1

### The Baron & Kenny approach is usually applied for continuous outcomes,
### using linear regressions. It is less adapted for binary outcomes.
### However, as for the binary mediator, using linear regression of the mediator
### and the outcome could still give some results.

### For binary outcomes, the Baron & Kenny approach
## Model 1: linear model of the probability of death to estimate the total effect:
model.tot.A.death <- lm(Y_death ~ A0_ace + L0_male + L0_parent_low_educ_lv,
                        data = df1)
summary(model.tot.A.death)
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)
#   (Intercept)           0.135282   0.007909  17.104  < 2e-16 ***
#   A0_ace                0.060247   0.012796   4.708 2.53e-06 *** <- Total effect
#   L0_male               0.050285   0.008053   6.244 4.43e-10 ***
#   L0_parent_low_educ_lv 0.059565   0.008425   7.070 1.65e-12 ***
# On a risk difference scale the total effect of ACE on the probability of death is
# approximately +6.0%

## Model 3: linear model to estimate the direct effect of the exposure (conditional on
## the outcome) and the effect of M on Y, adjusted for confounders of the A-Y
## and M-Y relationships
model.A.M.death <- lm(Y_death ~ A0_ace + M_smoking + L1 + L0_male + L0_parent_low_educ_lv,
                      data = df1)
summary(model.A.M.death)
#   Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)
#   (Intercept)           0.098691   0.008460  11.666  < 2e-16 ***
#   A0_ace                0.051509   0.012759   4.037 5.45e-05 *** <- Direct effect
#   M_smoking             0.064751   0.008822   7.340 2.31e-13 *** <- effect of M on Y
#   L1                    0.075533   0.008788   8.595  < 2e-16 ***
#   L0_male               0.047490   0.008013   5.927 3.19e-09 ***
#   L0_parent_low_educ_lv 0.055676   0.008389   6.637 3.36e-11 ***


# The direct effect is approximately +5.2% given by the A0_ace coefficient
model.A.M.death$coefficients["A0_ace"]
# 0.05150901

# The indirect effect can be calculated by the "difference in coefficient" method
# using coefficients from models 1 and 3
model.tot.A.death$coefficients["A0_ace"] - model.A.M.death$coefficients["A0_ace"]
# 0.008737889, i.e. approximately 0.9%

# or the product of coefficients using the previous model 2bis and model 3:
linear.model.A.M$coefficients["A0_ace"] * model.A.M.death$coefficients["M_smoking"]
# 0.008234978, i.e. approximately 0.8%


#### SEM
library(lavaan)
library(semPlot)

df1.qol.sem <- data.frame(L0_male = df1$L0_male,
                          L0_parent_low_educ_lv = df1$L0_parent_low_educ_lv,
                          A0_ace = df1$A0_ace, # better not to declare ordered(df1$A0_ace), #
                          L1 = df1$L1,
                          M_smoking = df1$M_smoking, # better not to declare ordered(df1$M_smoking),
                          Y_qol = df1$Y_qol)
sem.QoL <- '#SEM for quantitative outcome (QoL)
A0_ace ~ a.01 * L0_male
A0_ace ~ a.02 * L0_parent_low_educ_lv
M_smoking ~ b.L01 * L0_male
M_smoking ~ b.L02 * L0_parent_low_educ_lv
M_smoking ~ b.A * A0_ace
M_smoking ~ b.L1 * L1
Y_qol ~ g.01 * L0_male
Y_qol ~ g.02 * L0_parent_low_educ_lv
Y_qol ~ g.A * A0_ace
Y_qol ~ g.L1 * L1
Y_qol ~ g.M * M_smoking
# direct effect
dir := g.A
#indirect effect
ind := b.A * g.M
#total effect
total := (b.A * g.M) + g.A'

fit.qol <- sem(sem.QoL,
               data = df1.qol.sem)
summary(fit.qol)
# Defined Parameters:
#                Estimate  Std.Err  z-value  P(>|z|)
# dir              -3.965    0.321  -12.349    0.000
# ind              -1.104    0.129   -8.582    0.000
# total            -5.069    0.344  -14.753    0.000
# results are closer to the Baron & Kenny approach without declaring categorical variables... !

# semPaths(fit.qol,
#          what = "paths",
#          whatLabels = "par",
#          rotation = 2)


### binary
df1.death.sem <- data.frame(L0_male = df1$L0_male,
                            L0_parent_low_educ_lv = df1$L0_parent_low_educ_lv,
                            A0_ace = df1$A0_ace, #ordered(df1$A0_ace), #
                            L1 = df1$L1,
                            M_smoking = df1$M_smoking, #ordered(df1$M_smoking), #
                            Y_death = df1$Y_death) # ordered(df1$Y_death))

sem.death <- '#SEM for quantitative outcome (QoL)
A0_ace ~ a.01 * L0_male
A0_ace ~ a.02 * L0_parent_low_educ_lv
M_smoking ~ b.L01 * L0_male
M_smoking ~ b.L02 * L0_parent_low_educ_lv
M_smoking ~ b.A * A0_ace
M_smoking ~ b.L1 * L1
Y_death ~ g.01 * L0_male
Y_death ~ g.02 * L0_parent_low_educ_lv
Y_death ~ g.A * A0_ace
Y_death ~ g.L1 * L1
Y_death ~ g.M * M_smoking
# direct effect
dir := g.A
#indirect effect
ind := b.A * g.M
#total effect
total := (b.A * g.M) + g.A'

fit.death <- sem(sem.death,
                 data = df1.death.sem)
summary(fit.death)
# Defined Parameters:
#                Estimate  Std.Err  z-value  P(>|z|)
# dir               0.080    0.024    3.363    0.001
# ind               0.023    0.004    5.419    0.000
# total             0.103    0.023    4.465    0.000

# for binary variables => better to not declared them as ordered variables
# ordered variables are probably good for multicatogorical (>2) ordered variables,
# but not for binary variables
