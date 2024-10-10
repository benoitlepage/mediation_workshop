### test programs for estimations based on traditional regressions

rm(list=ls())

df1 <- read.csv(file = "df1.csv")
df1_int <- read.csv(file = "df1_int.csv")


################################################################################
######################### Estimation of the Average Total Effect (ATE)
################################################################################

### ATE
# For quantitative outcomes, apply a linear regression of Y on A (A0_PM2.5),
# adjusted for the baseline confounders L(0):
trad_ATE_qol <- lm(Y_qol ~ A0_PM2.5 + L0_male + L0_soc_env,
                   data = df1_int)

# For binary outcomes, apply a GLM of Y on A with a Gaussian distribution and
# identity link, adjusted for the baseline confounders:
trad_ATE_death <- glm(Y_death ~ A0_PM2.5 + L0_male + L0_soc_env,
                      family = gaussian("identity"),
                      data = df1_int)

# Use the regression coefficient of the exposure (A0_PM2.5) to estimate the ATE
ATE_trad_qol <- coefficients(trad_ATE_qol)["A0_PM2.5"]
ATE_trad_death <- coefficients(trad_ATE_death)["A0_PM2.5"]

ATE_trad_qol
# -7.210089

ATE_trad_death
# 0.07720726

## IC95%
library(sandwich)
# for the Quality of Life outcome
ATE_trad_qol <- list(ATE = coef(trad_ATE_qol)["A0_PM2.5"],
                     lo = coef(trad_ATE_qol)["A0_PM2.5"] - qnorm(0.975) *
                       sqrt(sandwich(trad_ATE_qol)["A0_PM2.5","A0_PM2.5"]),
                     hi = coef(trad_ATE_qol)["A0_PM2.5"] + qnorm(0.975) *
                       sqrt(sandwich(trad_ATE_qol)["A0_PM2.5","A0_PM2.5"]))

ATE_trad_qol
# $ATE
# A0_PM2.5
# -7.210089
# $lo
# [1] -7.978234
# $hi
# [1] -6.441944

# for death outcome
ATE_trad_death <- list(ATE = coef(trad_ATE_death)["A0_PM2.5"],
                       lo = coef(trad_ATE_death)["A0_PM2.5"] - qnorm(0.975) *
                         sqrt(sandwich(trad_ATE_death)["A0_PM2.5","A0_PM2.5"]),
                       hi = coef(trad_ATE_death)["A0_PM2.5"] + qnorm(0.975) *
                         sqrt(sandwich(trad_ATE_death)["A0_PM2.5","A0_PM2.5"]))
ATE_trad_death
# $ATE
# A0_PM2.5
# 0.07720726
# $lo
# A0_PM2.5
# 0.04945859
# $hi
# A0_PM2.5
# 0.1049559

# 95% CI calculation applying a bootstrap procedure
library(boot)
bootfunc <- function(data,index){
  boot_dat <- data[index,]
  mod.qol <- lm(Y_qol ~ A0_PM2.5 + L0_male + L0_soc_env,
                data = boot_dat)
  mod.death <- glm(Y_death ~ A0_PM2.5 + L0_male + L0_soc_env,
                   family = gaussian("identity"),
                   data = boot_dat)
  est <- c(coef(mod.qol)["A0_PM2.5"],
           coef(mod.death)["A0_PM2.5"])
  return(est)
}

set.seed(1234)
start.time <- Sys.time()
boot_est <- boot(df1_int,bootfunc,R=2000) # , parallel = "multicore" ça ne change rien sur le temps de calcul
end.time <- Sys.time()
end.time - start.time

# the 95% CI for the estimation of the ATE of ACE on QoL is:
boot.ci(boot_est, index = 1, type = "norm")
# Intervals :
# Level      Normal
# 95%   (-7.978, -6.444 )

# the 95% CI for the estimation of the ATE of ACE on death is:
boot.ci(boot_est, index = 2, type = "norm")
# Intervals :
# Level      Normal
# 95%   ( 0.0502,  0.1040 )


### For risk of death, expressed using Odds Ratio conditional on L(0):
TE_death_model <- glm(Y_death ~ A0_PM2.5 + L0_male + L0_soc_env,
                      family = "binomial",
                      data = df1_int)
res_TE_death <- summary(TE_death_model)
tot.effect.death.OR <- list(OR = exp(coef(res_TE_death)["A0_PM2.5","Estimate"]),
                            lo = exp(coef(res_TE_death)["A0_PM2.5","Estimate"] -
                                       qnorm(0.975) *
                                       coef(res_TE_death)["A0_PM2.5","Std. Error"]),
                            hi = exp(coef(res_TE_death)["A0_PM2.5","Estimate"] +
                                       qnorm(0.975) *
                                       coef(res_TE_death)["A0_PM2.5","Std. Error"]))
tot.effect.death.OR
# $OR
# A0_PM2.5
# 1.523254
# $lo
# [1] 1.323317
# $hi
# [1] 1.753398


################################################################################
######################### Two-way decomposition
################################################################################

################################################################################
######################### TNDE & PNIE, or PNDE & TNIE
################################################################################

### Quantitative outcome
trad_qol_am <- lm(Y_qol ~ A0_PM2.5 + M_diabetes + A0_PM2.5:M_diabetes +
                    L0_male + L0_soc_env + L1,
                  data = df1_int)
gamma.A.q <- coef(trad_qol_am)["A0_PM2.5"]
gamma.M.q <- coef(trad_qol_am)["M_diabetes"]
gamma.AM.q <- coef(trad_qol_am)["A0_PM2.5:M_diabetes"]


trad_m <- glm(M_diabetes ~ A0_PM2.5 + L0_male + L0_soc_env + L1, # il semble qu'il faut aussi ajuster sur L1 ici ? sur simulation, c'est mieux sans ajuster sur L1 # de plus, régression logistique => non-collapsibility cela peut jouer
              family = "binomial",
              data = df1_int)
beta.0 <- coef(trad_m)["(Intercept)"]
beta.A <- coef(trad_m)["A0_PM2.5"]


### binary outcome
trad_death_am <- glm(Y_death ~ A0_PM2.5 + M_diabetes + A0_PM2.5:M_diabetes +
                       L0_male + L0_soc_env + L1,
                     family = "binomial",
                     data = df1_int)
gamma.A.d <- coef(trad_death_am)["A0_PM2.5"]
gamma.M.d <- coef(trad_death_am)["M_diabetes"]
gamma.AM.d <- coef(trad_death_am)["A0_PM2.5:M_diabetes"]

### CDE
### For a continuous outcome
# setting the mediator to M=0
trad_CDE_qol_m0 <- gamma.A.q + gamma.AM.q * 0
trad_CDE_qol_m0
# -3.715265
# setting the mediator to M=1
trad_CDE_qol_m1 <- gamma.A.q + gamma.AM.q * 1
trad_CDE_qol_m1
# -9.330657

### For a binary outcome
## setting the mediator to M=0
trad_OD_CDE_death_m0 <- exp(gamma.A.d + gamma.AM.d * 0)
trad_OD_CDE_death_m0
# 1.442942

## setting the mediator to M=1
trad_OD_CDE_death_m1 <- exp(gamma.A.d + gamma.AM.d * 1)
trad_OD_CDE_death_m1
# 1.461464


### Natural direct and indirect effects
### For a continuous outcome
## The PNDE and TNIE are:
trad_PNDE_qol <- gamma.A.q + gamma.AM.q * (exp(beta.0)) / (1 + exp(beta.0))
trad_PNDE_qol
# -4.845089

trad_TNIE_qol <- (gamma.M.q + gamma.AM.q) *
  (exp(beta.0 + beta.A) / (1 + exp(beta.0 + beta.A)) -
     exp(beta.0) / (1 + exp(beta.0)))
trad_TNIE_qol
# -1.50119

## The TNDE and PNIE are:
trad_TNDE_qol <- gamma.A.q +
  gamma.AM.q * exp(beta.0 + beta.A) / (1 + exp(beta.0 + beta.A))
trad_TNDE_qol
# -5.436773

trad_PNIE_qol <- gamma.M.q *
  (exp(beta.0 + beta.A) /
     (1 + exp(beta.0 + beta.A)) - exp(beta.0) / (1 + exp(beta.0)))
trad_PNIE_qol
# -0.9095061



### For a binary outcome
## The PNDE and TNIE are:
trad_OR_PNDE_death <- exp(gamma.A.d) *
  (1 + exp(gamma.M.d + gamma.AM.d + beta.0 )) /
  (1 + exp(gamma.M.d + beta.0))
trad_OR_PNDE_death
# 1.448035

trad_OR_TNIE_death <- (1 + exp(beta.0)) *
  (1 + exp(gamma.M.d + gamma.AM.d + beta.0 + beta.A)) /
  ((1 + exp(beta.0 + beta.A)) * (1 + exp(gamma.M.d + gamma.AM.d + beta.0)))
trad_OR_TNIE_death
# 1.050029


## The TNDE and PNIE are:
trad_OR_TNDE_death <- exp(gamma.A.d) *
  (1 + exp(gamma.M.d + gamma.AM.d + beta.0 + beta.A)) /
  (1 + exp(gamma.M.d + beta.0 + beta.A))
trad_OR_TNDE_death
# 1.450344

trad_OR_PNIE_death <- (1 + exp(beta.0)) *
  (1 + exp(gamma.M.d + beta.0 + beta.A)) /
  ((1 + exp(beta.0 + beta.A)) * (1 + exp(gamma.M.d + beta.0)))
trad_OR_PNIE_death
# 1.048358


#### Using the regmedint package = Regression-Based Causal Mediation Analysis
#### with Interaction and Effect Modification Terms
library(regmedint)
## For continuous outcomes
regmedint_cont <- regmedint(data = df1_int,
                            ## Variables
                            yvar = "Y_qol",                   # outcome variable
                            avar = "A0_PM2.5",                  # exposure
                            mvar = "M_diabetes",               # mediator
                            cvar = c("L0_male",               # confounders
                                     "L0_soc_env",
                                     "L1"),
                            #eventvar = "event",     # only for survival outcome
                            ## Values at which effects are evaluated
                            a0 = 0,
                            a1 = 1,
                            m_cde = 0,
                            c_cond = c(0,0,0),                 # covariate level
                            ## Model types
                            mreg = "logistic",
                            yreg = "linear",
                            ## Additional specification
                            interaction = TRUE, # presence of A:M interaction term
                                                # in the outcome model
                            casecontrol = FALSE)
summary(regmedint_cont)
# ### Mediation analysis
#             est         se          Z            p      lower      upper
# cde  -3.7152652 0.41600219  -8.930879 0.000000e+00 -4.5306145 -2.8999159
# pnde -4.8450888 0.35052810 -13.822255 0.000000e+00 -5.5321113 -4.1580663
# tnie -1.5011902 0.20821830  -7.209694 5.608847e-13 -1.9092905 -1.0930898
# tnde -5.4367728 0.34049175 -15.967414 0.000000e+00 -6.1041244 -4.7694213
# pnie -0.9095061 0.12266064  -7.414817 1.219025e-13 -1.1499166 -0.6690957
# te   -6.3462790 0.38788368 -16.361294 0.000000e+00 -7.1065170 -5.5860409
# pm    0.2365465 0.02947624   8.024991 1.110223e-15  0.1787742  0.2943189

## For binary outcomes
regmedint_bin <- regmedint(data = df1_int,
                            ## Variables
                            yvar = "Y_death",                 # outcome variable
                            avar = "A0_PM2.5",                  # exposure
                            mvar = "M_diabetes",               # mediator
                            cvar = c("L0_male",               # confounders
                                     "L0_soc_env",
                                     "L1"),
                            #eventvar = "event",     # only for survival outcome
                            ## Values at which effects are evaluated
                            a0 = 0,
                            a1 = 1,
                            m_cde = 0,
                            c_cond = c(0,0,0),                 # covariate level
                            ## Model types
                            mreg = "logistic",
                            yreg = "logistic",
                            ## Additional specification
                            interaction = TRUE,
                            casecontrol = FALSE)
results.binary <- summary(regmedint_bin)
exp(results.binary$summary_myreg[,c("est","lower","upper")])
# ### Mediation analysis
#           est    lower    upper
# cde  1.442942 1.191195 1.747893
# pnde 1.448035 1.245470 1.683545
# tnie 1.050029 1.013842 1.087509
# tnde 1.450344 1.257042 1.673371
# pnie 1.048358 1.029285 1.067783
# te   1.520479 1.316954 1.755457
# pm   1.149340 1.031922 1.280118

################################################################################
######################### Three-way decomposition
################################################################################

### For a continuous outcome
## The PNDE is:
trad_PNDE_qol <- gamma.A.q + gamma.AM.q * (exp(beta.0)) / (1 + exp(beta.0))
trad_PNDE_qol
# -4.845089


## The PNIE is:
trad_PNIE_qol <- gamma.M.q *
  (exp(beta.0 + beta.A) /
     (1 + exp(beta.0 + beta.A)) - exp(beta.0) / (1 + exp(beta.0)))
trad_PNIE_qol
# -0.9095061

## The MIE is:
trad_MIE_qol <- gamma.AM.q *
  (exp(beta.0 + beta.A) / (1 + exp(beta.0 + beta.A)) -
     exp(beta.0) / (1 + exp(beta.0)))
trad_MIE_qol
# -0.591684


### For a binary outcome
## The excess relative risk is 52.3% (calculated from the OR of the total effect)
1.523254 - 1
# 0.523254

## The excess relative risk is decomposed into 3 components:
## The component of the excess relative risk due to PNDE is:
comp_PNDE_death <- exp(gamma.A.d) * (1 + exp(beta.0 + gamma.M.d + gamma.AM.d)) /
  (1 + exp(beta.0 + gamma.M.d)) - 1
comp_PNDE_death
# 0.4480347

## The component of the excess relative risk due to PNIE is:
comp_PNIE_death <- (1 + exp(beta.0)) * (1 + exp(beta.0 + beta.A + gamma.M.d)) /
  ((1 + exp(beta.0 + beta.A)) * (1 + exp(beta.0 + gamma.M.d))) - 1
comp_PNIE_death
# 0.04835753

## The component of the excess relative risk due to the mediated interactive
## effect is:
comp_MIE_qol <- exp(gamma.A.d) *
  (1 + exp(beta.0 + beta.A + gamma.M.d + gamma.AM.d)) * (1 + exp(beta.0)) /
  ((1 + exp(beta.0 + gamma.M.d)) * (1 + exp(beta.0 + beta.A))) -
  (1 + exp(beta.0 + beta.A + gamma.M.d)) * (1 + exp(beta.0)) /
  ((1 + exp(beta.0 + gamma.M.d)) * (1 + exp(beta.0 + beta.A))) -
  exp(gamma.A.d) * (1 + exp(beta.0 + gamma.M.d + gamma.AM.d)) /
  (1 + exp(beta.0 + gamma.M.d)) + 1
comp_MIE_qol
# 0.02408674





################################################################################
######################### Four-way decomposition
################################################################################

### For a continuous outcome
## The CDE_(M=0) is:
trad_CDE_qol_m0 <- gamma.A.q + gamma.AM.q * 0
trad_CDE_qol_m0
# -3.715265

## The PNIE is:
trad_PNIE_qol <- gamma.M.q *
  (exp(beta.0 + beta.A) /
     (1 + exp(beta.0 + beta.A)) - exp(beta.0) / (1 + exp(beta.0)))
trad_PNIE_qol
# -0.9095061

## The MIE is:
trad_MIE_qol <- gamma.AM.q *
  (exp(beta.0 + beta.A) / (1 + exp(beta.0 + beta.A)) -
     exp(beta.0) / (1 + exp(beta.0)))
trad_MIE_qol
# -0.591684

## The RIE is:
trad_RIE_qol <- gamma.AM.q * (exp(beta.0)) / (1 + exp(beta.0))
trad_RIE_qol
# -1.129824



### For a binary outcome
## The excess relative risk is 52.3% (calculated from the OR of the total effect)
1.523254 - 1
# 0.523254

## The excess relative risk is decomposed into 4 components:
## The component of the excess relative risk due to the CDE(M=0) is:
comp_CDE_death_m0 <- exp(gamma.A.d) * (1 + exp(beta.0)) /
  (1 + exp(beta.0 + gamma.M.d)) -
  (1 + exp(beta.0)) / (1 + exp(beta.0 + gamma.M.d))
comp_CDE_death_m0
# 0.402041

## The component of the excess relative risk due to the PNIE is:
comp_PNIE_death <- (1 + exp(beta.0)) * (1 + exp(beta.0 + beta.A + gamma.M.d)) /
  ((1 + exp(beta.0 + beta.A)) * (1 + exp(beta.0 + gamma.M.d))) - 1
comp_PNIE_death
# 0.04835753

## The component of the excess relative risk due to the MIE is:
comp_MIE_death <- exp(gamma.A.d) * (1 + exp(beta.0 + beta.A + gamma.M.d +
                                              gamma.AM.d)) * (1 + exp(beta.0)) /
  ((1 + exp(beta.0 + gamma.M.d)) * (1 + exp(beta.0 + beta.A))) -
  (1 + exp(beta.0 + beta.A + gamma.M.d)) * (1 + exp(beta.0)) /
  ((1 + exp(beta.0 + gamma.M.d)) * (1 + exp(beta.0 + beta.A))) -
  exp(gamma.A.d) * (1 + exp(beta.0 + gamma.M.d + gamma.AM.d)) /
  (1 + exp(beta.0 + gamma.M.d)) + 1
comp_MIE_death
# 0.02408674

## The component of the excess relative risk due to the RIE is:
comp_RIE_death <- exp(gamma.A.d) * (1 + exp(beta.0 + gamma.M.d + gamma.AM.d)) /
  (1 + exp(beta.0 + gamma.M.d)) - 1 -
  exp(gamma.A.d) * (1 + exp(beta.0)) / (1 + exp(beta.0 + gamma.M.d)) +
  (1 + exp(beta.0)) / (1 + exp(beta.0 + gamma.M.d))
comp_RIE_death
# 0.04599376



#### Using the CMAverse package = Regression-Based Causal Mediation Analysis
#### with Interaction and Effect Modification Terms
#### https://bs1125.github.io/CMAverse/
library(CMAverse)

### For the continuous outcome
## Closed-form parameter function estimation and delta method inferece
res_rb_param_delta <- cmest(data = df1_int,
                            model = "rb", # for "regression based" (rb) approach
                            outcome = "Y_qol",        # outcome variable
                            exposure = "A0_PM2.5",      # exposure variable
                            mediator = "M_diabetes",   # mediator
                            basec = c("L0_male",      # confounders
                                      "L0_soc_env",
                                      "L1"),
                            EMint = TRUE, # exposures*mediator interaction
                            mreg = list("logistic"), # model of the mediator
                            yreg = "linear",       # model of the outcome
                            astar = 0,
                            a = 1,
                            mval = list(0),
                            basecval = list(0,0,0),      # covariate level
                            estimation = "paramfunc", #  closed-form parameter
                                                      # function estimation
                            inference = "delta") # IC95% : delta method
summary(res_rb_param_delta)
# Closed-form parameter function estimation with
# delta method standard errors, confidence intervals and p-values
#
#                Estimate Std.error  95% CIL 95% CIU    P.val
#   cde          -3.71527   0.41600 -4.53061  -2.900  < 2e-16 ***   CDE(M=0)
#   pnde         -4.84509   0.35053 -5.53211  -4.158  < 2e-16 ***
#   tnde         -5.43677   0.34049 -6.10412  -4.769  < 2e-16 ***
#   pnie         -0.90951   0.12266 -1.14992  -0.669 1.22e-13 ***   PNIE
#   tnie         -1.50119   0.20822 -1.90929  -1.093 5.61e-13 ***
#   te           -6.34628   0.38788 -7.10652  -5.586  < 2e-16 ***
#   intref       -1.12982   0.13824 -1.40076  -0.859 2.22e-16 ***   INTref
#   intmed       -0.59168   0.10398 -0.79547  -0.388 1.27e-08 ***   MIE
#   cde(prop)     0.58542   0.04505  0.49712   0.674  < 2e-16 ***
#   intref(prop)  0.17803   0.02560  0.12786   0.228 3.51e-12 ***
#   intmed(prop)  0.09323   0.01586  0.06216   0.124 4.11e-09 ***
#   pnie(prop)    0.14331   0.01655  0.11087   0.176  < 2e-16 ***
#   pm            0.23655   0.02948  0.17877   0.294 1.11e-15 ***
#   int           0.27126   0.03780  0.19717   0.345 7.20e-13 ***
#   pe            0.41458   0.04505  0.32627   0.503  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# for the 3-way decomposition, the PNDE, PNIE, and MIE are given by:
data.frame("Estimate" = res_rb_param_delta$effect.pe,
           "lower95CI" = res_rb_param_delta$effect.ci.low,
           "upper95CI" = res_rb_param_delta$effect.ci.high,
           "P.value" = res_rb_param_delta$effect.pval)[c("pnde","pnie", "intmed"),]
#          Estimate  lower95CI  upper95CI      P.value
# pnde   -4.8450888 -5.5321113 -4.1580663 0.000000e+00
# pnie   -0.9095061 -1.1499166 -0.6690957 1.219025e-13
# intmed -0.5916840 -0.7954739 -0.3878941 1.266206e-08

# for the 4-way decomposition, the CDE(M=0), Intref, MIE and PNIE are given by:
data.frame("Estimate" = res_rb_param_delta$effect.pe,
           "lower95CI" = res_rb_param_delta$effect.ci.low,
           "upper95CI" = res_rb_param_delta$effect.ci.high,
           "P.value" = res_rb_param_delta$effect.pval)[c("cde","intref","intmed","pnie"),]
#          Estimate  lower95CI  upper95CI      P.value
# cde    -3.7152652 -4.5306145 -2.8999159 0.000000e+00
# intref -1.1298236 -1.4007600 -0.8588871 2.220446e-16
# intmed -0.5916840 -0.7954739 -0.3878941 1.266206e-08
# pnie   -0.9095061 -1.1499166 -0.6690957 1.219025e-13


### For the binary outcome
## Closed-form parameter function estimation and delta method inference
res_rb_param_delta <- cmest(data = df1_int,
                            model = "rb", # for "regression based" (rb) approach
                            outcome = "Y_death",        # outcome variable
                            exposure = "A0_PM2.5",      # exposure variable
                            mediator = "M_diabetes",   # mediator
                            basec = c("L0_male",      # confounders
                                      "L0_soc_env",
                                      "L1"),
                            EMint = TRUE, # exposures*mediator interaction
                            mreg = list("logistic"), # model of the mediator
                            yreg = "logistic",       # model of the outcome
                            astar = 0,
                            a = 1,
                            mval = list(0),
                            basecval = list(0,0,0),      # covariate level
                            estimation = "paramfunc", #  closed-form parameter
                            # function estimation
                            inference = "delta") # IC95% : delta method
summary(res_rb_param_delta)
# Closed-form parameter function estimation with
# delta method standard errors, confidence intervals and p-values
#
# Estimate Std.error  95% CIL 95% CIU    P.val
#   Rcde            1.44294   0.14115  1.19120   1.748 0.000178 ***
#   Rpnde           1.44803   0.11133  1.24547   1.684 1.47e-06 ***
#   Rtnde           1.45034   0.10585  1.25704   1.673 3.50e-07 ***
#   Rpnie           1.04836   0.00982  1.02929   1.068 4.62e-07 ***
#   Rtnie           1.05003   0.01879  1.01384   1.088 0.006368 **
#   Rte             1.52048   0.11148  1.31695   1.755 1.10e-08 ***
#   ERcde           0.40204   0.12699  0.15315   0.651 0.001546 **  CDE(M=0)
#   ERintref        0.04599   0.04821 -0.04850   0.140 0.340093     INTref
#   ERintmed        0.02409   0.02543 -0.02576   0.074 0.343607     MIE
#   ERpnie          0.04836   0.00982  0.02911   0.068 8.47e-07 *** PNIE
#   ERcde(prop)     0.77244   0.14282  0.49252   1.052 6.36e-08 ***
#   ERintref(prop)  0.08837   0.09311 -0.09413   0.271 0.342600
#   ERintmed(prop)  0.04628   0.04901 -0.04978   0.142 0.345058
#   ERpnie(prop)    0.09291   0.02602  0.04192   0.144 0.000356 ***
#   pm              0.13919   0.05498  0.03142   0.247 0.011359 *
#   int             0.13465   0.14186 -0.14340   0.413 0.342562
#   pe              0.22756   0.14282 -0.05237   0.507 0.111092
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# for the 3-way decomposition, the excess relative risk due to
# the PNDE, PNIE and MIE are given by:
data.frame("Estimate" = res_rb_param_delta$effect.pe - 1,
           "lower95CI" = res_rb_param_delta$effect.ci.low - 1,
           "upper95CI" = res_rb_param_delta$effect.ci.high - 1)[c("Rpnde"),]
#        Estimate lower95CI upper95CI
# Rpnde 0.4480347 0.2454699 0.6835449
# and
data.frame("Estimate" = res_rb_param_delta$effect.pe,
           "lower95CI" = res_rb_param_delta$effect.ci.low,
           "upper95CI" = res_rb_param_delta$effect.ci.high,
           "P.value" = res_rb_param_delta$effect.pval)[c("ERpnie", "ERintmed"),]
#            Estimate   lower95CI  upper95CI      P.value
# ERpnie   0.04835753  0.02910970 0.06760535 8.473169e-07
# ERintmed 0.02408674 -0.02576124 0.07393472 3.436070e-01




# for the 4-way decomposition, the CDE(M=0), the excess relative risk due to
# CDE(M=0), Intref, MIE and PNIE are given by:
data.frame("Estimate" = res_rb_param_delta$effect.pe,
           "lower95CI" = res_rb_param_delta$effect.ci.low,
           "upper95CI" = res_rb_param_delta$effect.ci.high,
           "P.value" = res_rb_param_delta$effect.pval)[c("ERcde","ERintref","ERintmed","ERpnie"),]
#            Estimate   lower95CI  upper95CI      P.value
# ERcde    0.40204098  0.15315072 0.65093124 1.545523e-03
# ERintref 0.04599376 -0.04850084 0.14048835 3.400930e-01
# ERintmed 0.02408674 -0.02576124 0.07393472 3.436070e-01
# ERpnie   0.04835753  0.02910970 0.06760535 8.473169e-07


