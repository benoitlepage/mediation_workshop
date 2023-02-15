### test programs for estimations based on IPTW

rm(list=ls())

df2 <- read.csv(file = "df2.csv")
df2_int <- read.csv(file = "df2_int.csv")


################################################################################
######################### Estimation of the Average Total Effect (ATE)
################################################################################

### IPTW for the ATE

# 1. Estimate g
g.L <- glm(A0_ace ~ L0_male + L0_parent_low_educ_lv,
           family = "binomial", data = df2_int)

# 2. Predict each individual's probability of being exposed to her own exposure
# predict the probabilities P(A0_ace=1|L(0)) & P(A0_ace=0|L(0))
pred.g1.L <- predict(g.L, type="response")
pred.g0.L <- 1 - pred.g1.L
# the predicted probability of the observed treatment A_i=a is :
gA.L <- rep(NA, nrow(df2_int))
gA.L[df2_int$A0_ace==1] <- pred.g1.L[df2_int$A0_ace==1]
gA.L[df2_int$A0_ace==0] <- pred.g0.L[df2_int$A0_ace==0]

# 3. Apply weights corresponding to the inverse of the predicted probability
wt <- 1 / gA.L

# 4. Use the empirical mean of the weighted outcome
# point estimates:
IPTW.death <- mean(wt * as.numeric(df2_int$A0_ace==1) * df2_int$Y_death) -
  mean(wt * as.numeric(df2_int$A0_ace==0) * df2_int$Y_death)
IPTW.death
# [1] 0.08224947

IPTW.qol <- mean(wt * as.numeric(df2_int$A0_ace==1) * df2_int$Y_qol) -
  mean(wt * as.numeric(df2_int$A0_ace==0) * df2_int$Y_qol)
IPTW.qol
# [1] -8.436797



### Stabilized IPTW for the ATE

# 3. For example, applying g^*(A) = 1
# 4. Applying the stabilized estimator
# point estimates:
s.IPTW.death <- (mean(wt * as.numeric(df2_int$A0_ace==1) * df2_int$Y_death) /
                   mean(wt * as.numeric(df2_int$A0_ace==1))) -
  (mean(wt * as.numeric(df2_int$A0_ace==0) * df2_int$Y_death) /
     mean(wt * as.numeric(df2_int$A0_ace==0)))
s.IPTW.death
# [1] 0.08294185

s.IPTW.qol <- (mean(wt * as.numeric(df2_int$A0_ace==1) * df2_int$Y_qol) /
                 mean(wt * as.numeric(df2_int$A0_ace==1))) -
  (mean(wt * as.numeric(df2_int$A0_ace==0) * df2_int$Y_qol) /
     mean(wt * as.numeric(df2_int$A0_ace==0)))
s.IPTW.qol
# [1] -8.291992
