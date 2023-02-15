param.causal.model.1 <- function(A.M.interaction = NULL) {
  # L0
  p_L0_male <- 0.5
  p_L0_parent_low_educ_lv <- 0.65

  # A: A0_ace <- rbinom( 0.05 + 0.04 * L0_male + 0.06 * L0_parent_low_educ_lv )
  b_A <- 0.05   # reference prevalence is 5%
  b_male_A <- 0.04  # + 0.04 for the effect of L0_male -> A0_ace
  b_parent_educ_A <- 0.06  # +0.06 for the effect of L0_parent_low_educ_lv -> A0_ace

  # L1: intermediate confounder between M and Y, not influenced by A
  p_L1 <- 0.3

  # M: M_smoking <- rbinom( 0.2 + 0.05 * L0_male + 0.06 * L0_parent_low_educ_lv + 0.07 * L1 +
  #                         0.1 * A0_ace )
  b_M <- 0.2 # reference prevalence is 20%
  b_male_M <- 0.05 # +0.05 for the effect of L0_male -> M_smoking
  b_parent_educ_M <- 0.06 # +0.06 for the effect of L0_parent_low_educ_lv -> M_smoking
  b_L1_M <- 0.07 # +0.07 for the effect of L1 -> M_smoking
  b_A_M <- 0.1 # +0.10 for the effect of A0_ace -> M_smoking

  # Y binary: rbinom( 0.10 + 0.06 * L0_male + 0.04 * L0_parent_low_educ_lv + 0.05 * A0_ace +
  #                   0.07 * L1 + 0.08 * M_smoking +
  #                   0.03 * A0_ace * M_smoking * A.M.inter )
  b_Y <- 0.1 # reference prevalence is 10%
  b_male_Y <- 0.06 # +0.06 for the effect of L0_male -> Y
  b_parent_educ_Y <- 0.04 # +0.04 for the effect of L0_parent_low_educ_lv -> Y
  b_A_Y <- 0.05 # 0.05 for the effect of A0_ace -> Y
  b_L1_Y <- 0.07 # +0.07 for the effect of L1 -> Y
  b_M_Y <- 0.08 # 0.08 for the effect of M_smoking -> Y
  b_AM_Y <- 0.03 # 0.03 for the interaction effect A0_ace * M_smoking -> Y

  # Y continuous: (75 - 1 * L0_male - 3 * L0_parent_low_educ_lv - 4 * A0_ace -3.5 * L1 -
  #                9 * M_smoking -5 * A0_ace * M_smoking * A.M.inter ) +
  #                rnorm(N, mean = 0, sd = 10)
  mu_Y <- 75 # reference mean for QoL
  c_male_Y <- -1 # -1 for the effect of L0_male -> Y
  c_parent_educ_Y <- -3 # -3 for the effect of L0_parent_low_educ_lv -> Y
  c_A_Y <- -4 # -4 for the effect of A0_ace -> Y
  c_L1_Y <- -3.5 # -3.5 for the effect of L1 -> Y
  c_M_Y <- -9 # -9 for the effect of M_smoking -> Y
  c_AM_Y <- -5  # - 5 for the interaction effect A0_ace * M_smoking  -> Y
  sd_Y <- 10 # standard deviation of the residuals

  # A*M interaction ?
  A.M.inter <- A.M.interaction

  coef <- c( p_L0_male = p_L0_male, p_L0_parent_low_educ_lv = p_L0_parent_low_educ_lv,
             b_A = b_A, b_male_A = b_male_A, b_parent_educ_A = b_parent_educ_A,
             p_L1 = p_L1,
             b_M = b_M, b_male_M = b_male_M, b_parent_educ_M = b_parent_educ_M,
             b_L1_M = b_L1_M, b_A_M = b_A_M,
             b_Y = b_Y, b_male_Y = b_male_Y, b_parent_educ_Y = b_parent_educ_Y,
             b_A_Y = b_A_Y, b_L1_Y = b_L1_Y, b_M_Y = b_M_Y, b_AM_Y = b_AM_Y,
             mu_Y = mu_Y, c_male_Y = c_male_Y, c_parent_educ_Y = c_parent_educ_Y,
             c_A_Y = c_A_Y, c_L1_Y = c_L1_Y, c_M_Y = c_M_Y, c_AM_Y = c_AM_Y,
             sd_Y = sd_Y, A.M.inter = A.M.inter)

  return(coef)
}


gen.data.causal.model.1 <- function(N, A.M.inter) { # input parameters are the
  #   sample size N and the presence of A*M interaction with A.M.inter = 0 or 1

  b <- param.causal.model.1(A.M.interaction = A.M.inter)

  # baseline confounders: parent's educational level=L0_parent_low_educ_lv & sex=L0_male
  L0_male <- rbinom(N, size = 1, prob = b["p_L0_male"])
  L0_parent_low_educ_lv <- rbinom(N, size = 1, prob = b["p_L0_parent_low_educ_lv"])

  # exposure: A0_ace
  A0_ace <- rbinom(N, size = 1, prob =  b["b_A"] +
                     b["b_male_A"] * L0_male +
                     b["b_parent_educ_A"] * L0_parent_low_educ_lv )

  # intermediate confounder between M_smoking and Y, not affected by A0 L1
  L1 <- rbinom(N, size = 1, prob = b["p_L1"])

  # mediator: M_smoking
  M_smoking <- rbinom(N, size = 1, prob = b["b_M"] +
                        b["b_male_M"] * L0_male +
                        b["b_parent_educ_M"] * L0_parent_low_educ_lv +
                        b["b_A_M"] * A0_ace +
                        b["b_L1_M"] * L1)

  # Y_death
  Y_death <- rbinom(N, size = 1, prob = b["b_Y"] +
                      b["b_male_Y"] * L0_male +
                      b["b_parent_educ_Y"] * L0_parent_low_educ_lv +
                      b["b_A_Y"] * A0_ace +
                      b["b_L1_Y"] * L1 +
                      b["b_M_Y"] * M_smoking +
                      b["b_AM_Y"] * A0_ace * M_smoking * A.M.inter )

  # Y_qol
  Y_qol <- ( b["mu_Y"] +
               b["c_male_Y"] * L0_male +
               b["c_parent_educ_Y"] * L0_parent_low_educ_lv +
               b["c_A_Y"] * A0_ace +
               b["c_L1_Y"] * L1 +
               b["c_M_Y"] * M_smoking +
               b["c_AM_Y"] * A0_ace * M_smoking * A.M.inter ) +
    rnorm(N, mean = 0, sd = b["sd_Y"])

  # data.frame
  data.sim <- data.frame(L0_male, L0_parent_low_educ_lv, A0_ace, L1, M_smoking,
                         Y_death, Y_qol)

  return( data.sim )
}


######## Effet total
## note : l'effet théorique P(Y_a | l0) dans le sous-groupe où l0 = 0  est =
## sum_{m,l1} P(Y=1|a,m,l0,l1)P(m|a,l0,l1)P(l1)
true.OR_ATE1 <- function(interaction = NULL) {
  b <- param.causal.model.1(A.M.interaction = interaction)

  S <- cbind(expand.grid(c(0,1),c(0,1)), rep(NA,n=2^2))
  colnames(S) <- list("L1","M","sum")
  for (n in 1:4) {
    S[n,"sum"] <- ( b["b_Y"] +
                      b["b_male_Y"] * 0 +
                      b["b_parent_educ_Y"] * 0 +
                      b["b_A_Y"] * 1 +
                      b["b_L1_Y"] * S[n,"L1"] +
                      b["b_M_Y"] * S[n,"M"] +
                      b["b_AM_Y"] * 1 * S[n,"M"] * b["A.M.inter"] ) *
      (( b["b_M"] +
           b["b_male_M"] * 0 +
           b["b_parent_educ_M"] * 0 +
           b["b_L1_M"] * S[n,"L1"] +
           b["b_A_M"] * 1 )^( S[n,"M"] )) *
      (( 1 - (b["b_M"] +
                b["b_male_M"] * 0 +
                b["b_parent_educ_M"] * 0 +
                b["b_L1_M"] * S[n,"L1"] +
                b["b_A_M"] * 1) )^( 1 - S[n,"M"] )) *
      ((b["p_L1"])^(S[n,"L1"])) *
      ((1 - b["p_L1"])^(1 - S[n,"L1"]))
  }
  EY_A1 <- sum(S[,"sum"])

  S <- cbind(expand.grid(c(0,1),c(0,1)), rep(NA,n=2^2))
  colnames(S) <- list("L1","M","sum")
  for (n in 1:4) {
    S[n,"sum"] <- ( b["b_Y"] +
                    b["b_male_Y"] * 0 +
                    b["b_parent_educ_Y"] * 0 +
                    b["b_A_Y"] * 0 +
                    b["b_L1_Y"] * S[n,"L1"] +
                    b["b_M_Y"] * S[n,"M"] +
                    b["b_AM_Y"] * 0 * S[n,"M"] * b["A.M.inter"] ) *
      (( b["b_M"] +
         b["b_male_M"] * 0 +
         b["b_parent_educ_M"] * 0 +
         b["b_L1_M"] * S[n,"L1"] +
         b["b_A_M"] * 0 )^( S[n,"M"] )) *
      (( 1 - (b["b_M"] +
              b["b_male_M"] * 0 +
              b["b_parent_educ_M"] * 0 +
              b["b_L1_M"] * S[n,"L1"] +
              b["b_A_M"] * 0) )^( 1 - S[n,"M"] )) *
      ((b["p_L1"])^(S[n,"L1"])) *
      ((1 - b["p_L1"])^(1 - S[n,"L1"]))
    }
  EY_A0 <- sum(S[,"sum"])

  OR_ATE = (EY_A1/(1-EY_A1)) / (EY_A0/(1-EY_A0))
  return = list(EY_A1 = EY_A1, EY_A0 = EY_A0, OR_ATE = OR_ATE)
}

true.OR_TE1.with.inter <- true.OR_ATE1(interaction = 1)
true.OR_TE1.with.inter
# $EY_A1
# [1] 0.20631
#
# $EY_A0
# [1] 0.13868
#
# $OR_ATE
# [1] 1.614433

df1_int <- gen.data.causal.model.1(N = 1e7, A.M.inter =1)
TE_death_model <- glm(Y_death ~ A0_ace * L0_male +
                        A0_ace * L0_parent_low_educ_lv +
                        L0_male * L0_parent_low_educ_lv +
                        A0_ace:L0_male:L0_parent_low_educ_lv,
                      family = "binomial",
                      data = df1_int)
exp(coef(TE_death_model)["A0_ace"])
# 1.615188 OK là ça marche ...
base.test <- df1_int[1:2,]
base.test$L0_male <- rep(0, 2)
base.test$L0_parent_low_educ_lv <- rep(0, 2)
base.test$A0_ace[1] <- 1
base.test$A0_ace[2] <- 0
predict(TE_death_model, newdata = base.test, type = "response")
# 1         2
# 0.2076379 0.1384742
# le problème est que c'est un modèle difficile à prédire
# il y a besoin de interactions de degré 3 car il n'est pas simulé sur une échelle multiplicative...
# avec la triple interaction, je retrouve bien les bonnes valeurs ...

# Donc le modèle simple de VanderWeele est biaisé car les données sont simulées à partir d'un modèle additif !



########## Estimation des vrais OR_CDE_conditionnels à L0=0 et L1=0
true.OR_CDEm <- function(interaction = NULL) {
  b <- param.causal.model.1(A.M.interaction = interaction)

  S <- cbind(expand.grid(c(0,1),c(0,1)), rep(NA,n=2^2))
  colnames(S) <- list("a","m","E_Yam")
  for (n in 1:4) {
    S[n,"E_Yam"] <- ( b["b_Y"] +
                      b["b_male_Y"] * 0 +
                      b["b_parent_educ_Y"] * 0 +
                      b["b_A_Y"] * S[n,"a"] +
                      b["b_L1_Y"] * 0 +
                      b["b_M_Y"] * S[n,"m"] +
                      b["b_AM_Y"] * S[n,"a"] * S[n,"m"] * b["A.M.inter"] )
  }
  EY_00 <- S[S["a"] == 0 & S["m"] == 0,"E_Yam"]
  EY_01 <- S[S["a"] == 0 & S["m"] == 1,"E_Yam"]
  EY_10 <- S[S["a"] == 1 & S["m"] == 0,"E_Yam"]
  EY_11 <- S[S["a"] == 1 & S["m"] == 1,"E_Yam"]

  OR_CDE_m0 = (EY_10/(1-EY_10)) / (EY_00/(1-EY_00))
  OR_CDE_m1 = (EY_11/(1-EY_11)) / (EY_01/(1-EY_01))
  return = list(OR_CDE_m0 = OR_CDE_m0, OR_CDE_m1 = OR_CDE_m1)
}

true.OR_CDE.with.inter <- true.OR_CDEm(interaction = 1)
true.OR_CDE.with.inter
# $OR_CDE_m0
# [1] 1.588235
#
# $OR_CDE_m1
# [1] 1.600601


########## Estimation des vrais OR_PNDE_TNIE_conditionnels à L0=0 et L1=0
true.OR_NDE.NIE1 <- function(interaction = NULL) {
  b <- param.causal.model.1(A.M.interaction = interaction)

  # binary outcome (death)
  S <- cbind(expand.grid(c(0,1)), rep(NA,n=2^4), rep(NA,n=2^4), rep(NA,n=2^4),
             rep(NA,n=2^4))
  colnames(S) <- list("M","EY_1.1", "EY_1.0", "EY_0.1", "EY_0.0")
  for (n in 1:2) {
    # PNDE
    S[n,"EY_1.1"] <- ( b["b_Y"] +
                         b["b_male_Y"] * 0 +
                         b["b_parent_educ_Y"] * 0 +
                         b["b_A_Y"] * 1 +
                         b["b_L1_Y"] * 0 +
                         b["b_M_Y"] * S[n,"M"] +
                         b["b_AM_Y"] * 1 * S[n,"M"] * b["A.M.inter"] ) *
      (( b["b_M"] +
           b["b_male_M"] * 0 +
           b["b_parent_educ_M"] * 0 +
           b["b_L1_M"] * 0 +
           b["b_A_M"] * 1 )^( S[n,"M"] )) *
      (( 1 - (b["b_M"] +
                b["b_male_M"] * 0 +
                b["b_parent_educ_M"] * 0 +
                b["b_L1_M"] * 0 +
                b["b_A_M"] * 1) )^( 1 - S[n,"M"] ))

    S[n,"EY_1.0"] <- ( b["b_Y"] +
                         b["b_male_Y"] * 0 +
                         b["b_parent_educ_Y"] * 0 +
                         b["b_A_Y"] * 1 +
                         b["b_L1_Y"] * 0 +
                         b["b_M_Y"] * S[n,"M"] +
                         b["b_AM_Y"] * 1 * S[n,"M"] * b["A.M.inter"] ) *
      (( b["b_M"] +
           b["b_male_M"] * 0 +
           b["b_parent_educ_M"] * 0 +
           b["b_L1_M"] * 0 +
           b["b_A_M"] * 0 )^( S[n,"M"] )) *
      (( 1 - (b["b_M"] +
                b["b_male_M"] * 0 +
                b["b_parent_educ_M"] * 0 +
                b["b_L1_M"] * 0 +
                b["b_A_M"] * 0) )^( 1 - S[n,"M"] ))

    S[n,"EY_0.1"] <- ( b["b_Y"] +
                         b["b_male_Y"] * 0 +
                         b["b_parent_educ_Y"] * 0 +
                         b["b_A_Y"] * 0 +
                         b["b_L1_Y"] * 0 +
                         b["b_M_Y"] * S[n,"M"] +
                         b["b_AM_Y"] * 0 * S[n,"M"] * b["A.M.inter"] ) *
      (( b["b_M"] +
           b["b_male_M"] * 0 +
           b["b_parent_educ_M"] * 0 +
           b["b_L1_M"] * 0 +
           b["b_A_M"] * 1 )^( S[n,"M"] )) *
      (( 1 - (b["b_M"] +
                b["b_male_M"] * 0 +
                b["b_parent_educ_M"] * 0 +
                b["b_L1_M"] * 0 +
                b["b_A_M"] * 1) )^( 1 - S[n,"M"] ))

    S[n,"EY_0.0"] <- ( b["b_Y"] +
                         b["b_male_Y"] * 0 +
                         b["b_parent_educ_Y"] * 0 +
                         b["b_A_Y"] * 0 +
                         b["b_L1_Y"] * 0 +
                         b["b_M_Y"] * S[n,"M"] +
                         b["b_AM_Y"] * 0 * S[n,"M"] * b["A.M.inter"] ) *
      (( b["b_M"] +
           b["b_male_M"] * 0 +
           b["b_parent_educ_M"] * 0 +
           b["b_L1_M"] * 0 +
           b["b_A_M"] * 0 )^( S[n,"M"] )) *
      (( 1 - (b["b_M"] +
                b["b_male_M"] * 0 +
                b["b_parent_educ_M"] * 0 +
                b["b_L1_M"] * 0 +
                b["b_A_M"] * 0) )^( 1 - S[n,"M"] ))

  }

  EY_11 <- sum(S[,"EY_1.1"])
  EY_10 <- sum(S[,"EY_1.0"])
  EY_01 <- sum(S[,"EY_0.1"])
  EY_00 <- sum(S[,"EY_0.0"])

  OR_PNDE = (EY_10/(1-EY_10)) / (EY_00/(1-EY_00))
  OR_TNIE = (EY_11/(1-EY_11)) / (EY_10/(1-EY_10))
  OR_TNDE = (EY_11/(1-EY_11)) / (EY_01/(1-EY_01))
  OR_PNIE = (EY_01/(1-EY_01)) / (EY_00/(1-EY_00))

  return = list(OR_PNDE = OR_PNDE, OR_TNIE = OR_TNIE,
                OR_TNDE = OR_TNDE, OR_PNIE = OR_PNIE)
}

true.OR_NDE.NIE1.with.inter <- true.OR_NDE.NIE1(interaction = 1)
true.OR_NDE.NIE1.with.inter
# $OR_PNDE
# [1] 1.583042
#
# $OR_TNIE
# [1] 1.078278
#
# $OR_TNDE
# [1] 1.582382
#
# $OR_PNIE
# [1] 1.078728
# Les résultats sont très proches entre pure et total

################################################################################
### simulation
ATE_biasedQ <- rep(NA, 1000)
ATE_correct <- rep(NA, 1000)

biased_CDE_M0_ATE <- rep(NA, 1000)
biased_CDE_M1_ATE <- rep(NA, 1000)
flexible_CDE_M0_ATE <- rep(NA, 1000)
flexible_CDE_M1_ATE <- rep(NA, 1000)
satu_CDE_M0_ATE <- rep(NA, 1000)
satu_CDE_M1_ATE <- rep(NA, 1000)

biased_PNDE <- rep(NA, 1000)
biased_TNIE <- rep(NA, 1000)
flexible_PNDE <- rep(NA, 1000)
flexible_TNIE <- rep(NA, 1000)
satu_PNDE <- rep(NA, 1000)
satu_TNIE <- rep(NA, 1000)

biased_TNDE <- rep(NA, 1000)
biased_PNIE <- rep(NA, 1000)
flexible_TNDE <- rep(NA, 1000)
flexible_PNIE <- rep(NA, 1000)
satu_TNDE <- rep(NA, 1000)
satu_PNIE <- rep(NA, 1000)

set.seed(1234)
for (j in 1:1000){
  print(paste("je suis à la simulation", j))
  df1_int <- gen.data.causal.model.1(N = 10000, A.M.inter =1)

  ### effet total
  TE_death_biased_model <- glm(Y_death ~ A0_ace + L0_male + L0_parent_low_educ_lv,
                                family = "binomial",
                                data = df1_int)
  ATE_biasedQ[j] <- exp(coef(TE_death_biased_model)["A0_ace"])

  TE_death_correct_model <- glm(Y_death ~ A0_ace * L0_male * L0_parent_low_educ_lv,
                        family = "binomial",
                        data = df1_int)
  ATE_correct[j] <- exp(coef(TE_death_correct_model)["A0_ace"])

  ### CDE
  trad_death_biased_CDEmodl <- glm(Y_death ~ A0_ace + M_smoking + A0_ace:M_smoking +
                                     L0_male + L0_parent_low_educ_lv + L1,
                                   family = "binomial",
                                   data = df1_int)

  trad_death_moreflex_CDEmodl <- glm(Y_death ~ (A0_ace + M_smoking + L0_male + L0_parent_low_educ_lv + L1)^2,
                                     family = "binomial",
                                     data = df1_int)

  trad_death_satu_CDEmodl <- glm(Y_death ~ A0_ace * M_smoking * L0_male * L0_parent_low_educ_lv * L1,
                                     family = "binomial",
                                     data = df1_int)

  biased_CDE_M0_ATE[j] <- exp(coef(trad_death_biased_CDEmodl)["A0_ace"] +
                      coef(trad_death_biased_CDEmodl)["A0_ace:M_smoking"] * 0)

  biased_CDE_M1_ATE[j] <- exp(coef(trad_death_biased_CDEmodl)["A0_ace"] +
                                coef(trad_death_biased_CDEmodl)["A0_ace:M_smoking"] * 1)

  flexible_CDE_M0_ATE[j] <- exp(coef(trad_death_moreflex_CDEmodl)["A0_ace"] +
                                coef(trad_death_moreflex_CDEmodl)["A0_ace:M_smoking"] * 0)

  flexible_CDE_M1_ATE[j] <- exp(coef(trad_death_moreflex_CDEmodl)["A0_ace"] +
                                  coef(trad_death_moreflex_CDEmodl)["A0_ace:M_smoking"] * 1)

  satu_CDE_M0_ATE[j] <- exp(coef(trad_death_satu_CDEmodl)["A0_ace"] +
                              coef(trad_death_satu_CDEmodl)["A0_ace:M_smoking"] * 0)

  satu_CDE_M1_ATE[j] <- exp(coef(trad_death_satu_CDEmodl)["A0_ace"] +
                              coef(trad_death_satu_CDEmodl)["A0_ace:M_smoking"] * 1)

  # PNDE and TNIE
  trad_mediator_biased <- glm(M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv + L1,
                              family = "binomial",
                              data = df1_int)

  trad_mediator_moreflex <- glm(M_smoking ~ (A0_ace + L0_male + L0_parent_low_educ_lv + L1)^2,
                              family = "binomial",
                              data = df1_int)

  trad_mediator_satu <- glm(M_smoking ~ A0_ace * L0_male * L0_parent_low_educ_lv * L1,
                              family = "binomial",
                              data = df1_int)

  biased_PNDE[j] <- exp(coef(trad_death_biased_CDEmodl)["A0_ace"]) *
    (1 + exp(coef(trad_death_biased_CDEmodl)["M_smoking"] +
               coef(trad_death_biased_CDEmodl)["A0_ace:M_smoking"] +
               coef(trad_mediator_biased)["(Intercept)"])) /
    (1 + exp(coef(trad_death_biased_CDEmodl)["M_smoking"] +
               coef(trad_mediator_biased)["(Intercept)"]))

  biased_TNIE[j] <- (1 + exp(coef(trad_mediator_biased)["(Intercept)"])) *
    (1 + exp(coef(trad_death_biased_CDEmodl)["M_smoking"] +
               coef(trad_death_biased_CDEmodl)["A0_ace:M_smoking"] +
               coef(trad_mediator_biased)["(Intercept)"] +
               coef(trad_mediator_biased)["A0_ace"])) /
    ((1 + exp(coef(trad_mediator_biased)["(Intercept)"] + coef(trad_mediator_biased)["A0_ace"])) *
       (1 + exp(coef(trad_death_biased_CDEmodl)["M_smoking"] +
                  coef(trad_death_biased_CDEmodl)["A0_ace:M_smoking"] +
                  coef(trad_mediator_biased)["(Intercept)"])))

  flexible_PNDE[j] <- exp(coef(trad_death_moreflex_CDEmodl)["A0_ace"]) *
    (1 + exp(coef(trad_death_moreflex_CDEmodl)["M_smoking"] +
               coef(trad_death_moreflex_CDEmodl)["A0_ace:M_smoking"] +
               coef(trad_mediator_moreflex)["(Intercept)"])) /
    (1 + exp(coef(trad_death_moreflex_CDEmodl)["M_smoking"] +
               coef(trad_mediator_moreflex)["(Intercept)"]))

  flexible_TNIE[j] <- (1 + exp(coef(trad_mediator_moreflex)["(Intercept)"])) *
    (1 + exp(coef(trad_death_moreflex_CDEmodl)["M_smoking"] +
               coef(trad_death_moreflex_CDEmodl)["A0_ace:M_smoking"] +
               coef(trad_mediator_moreflex)["(Intercept)"] +
               coef(trad_mediator_moreflex)["A0_ace"])) /
    ((1 + exp(coef(trad_mediator_moreflex)["(Intercept)"] + coef(trad_mediator_moreflex)["A0_ace"])) *
       (1 + exp(coef(trad_death_moreflex_CDEmodl)["M_smoking"] +
                  coef(trad_death_moreflex_CDEmodl)["A0_ace:M_smoking"] +
                  coef(trad_mediator_moreflex)["(Intercept)"])))

  satu_PNDE[j] <- exp(coef(trad_death_satu_CDEmodl)["A0_ace"]) *
    (1 + exp(coef(trad_death_satu_CDEmodl)["M_smoking"] +
               coef(trad_death_satu_CDEmodl)["A0_ace:M_smoking"] +
               coef(trad_mediator_satu)["(Intercept)"])) /
    (1 + exp(coef(trad_death_satu_CDEmodl)["M_smoking"] +
               coef(trad_mediator_satu)["(Intercept)"]))

  satu_TNIE[j] <- (1 + exp(coef(trad_mediator_satu)["(Intercept)"])) *
    (1 + exp(coef(trad_death_satu_CDEmodl)["M_smoking"] +
               coef(trad_death_satu_CDEmodl)["A0_ace:M_smoking"] +
               coef(trad_mediator_satu)["(Intercept)"] +
               coef(trad_mediator_satu)["A0_ace"])) /
    ((1 + exp(coef(trad_mediator_satu)["(Intercept)"] + coef(trad_mediator_satu)["A0_ace"])) *
       (1 + exp(coef(trad_death_satu_CDEmodl)["M_smoking"] +
                  coef(trad_death_satu_CDEmodl)["A0_ace:M_smoking"] +
                  coef(trad_mediator_satu)["(Intercept)"])))

  # TNDE and PNIE
    biased_TNDE[j] <- exp(coef(trad_death_biased_CDEmodl)["A0_ace"]) *
      (1 + exp(coef(trad_death_biased_CDEmodl)["M_smoking"] +
                 coef(trad_death_biased_CDEmodl)["A0_ace:M_smoking"] +
                 coef(trad_mediator_biased)["(Intercept)"] +
                 coef(trad_mediator_biased)["A0_ace"])) /
      (1 + exp(coef(trad_death_biased_CDEmodl)["M_smoking"] +
                 coef(trad_mediator_biased)["(Intercept)"] + coef(trad_mediator_biased)["A0_ace"]))

    biased_PNIE[j] <- (1 + exp(coef(trad_mediator_biased)["(Intercept)"])) *
      (1 + exp(coef(trad_death_biased_CDEmodl)["M_smoking"] +
                 coef(trad_mediator_biased)["(Intercept)"] +
                 coef(trad_mediator_biased)["A0_ace"])) /
      ((1 + exp(coef(trad_mediator_biased)["(Intercept)"] + coef(trad_mediator_biased)["A0_ace"])) *
         (1 + exp(coef(trad_death_biased_CDEmodl)["M_smoking"] +
                    coef(trad_mediator_biased)["(Intercept)"])))

    flexible_TNDE[j] <- exp(coef(trad_death_moreflex_CDEmodl)["A0_ace"]) *
      (1 + exp(coef(trad_death_moreflex_CDEmodl)["M_smoking"] +
                 coef(trad_death_moreflex_CDEmodl)["A0_ace:M_smoking"] +
                 coef(trad_mediator_moreflex)["(Intercept)"] +
                 coef(trad_mediator_moreflex)["A0_ace"])) /
      (1 + exp(coef(trad_death_moreflex_CDEmodl)["M_smoking"] +
                 coef(trad_mediator_moreflex)["(Intercept)"] + coef(trad_mediator_moreflex)["A0_ace"]))

    flexible_PNIE[j] <- (1 + exp(coef(trad_mediator_moreflex)["(Intercept)"])) *
      (1 + exp(coef(trad_death_moreflex_CDEmodl)["M_smoking"] +
                 coef(trad_mediator_moreflex)["(Intercept)"] +
                 coef(trad_mediator_moreflex)["A0_ace"])) /
      ((1 + exp(coef(trad_mediator_moreflex)["(Intercept)"] + coef(trad_mediator_moreflex)["A0_ace"])) *
         (1 + exp(coef(trad_death_moreflex_CDEmodl)["M_smoking"] +
                    coef(trad_mediator_moreflex)["(Intercept)"])))

    satu_TNDE[j] <- exp(coef(trad_death_satu_CDEmodl)["A0_ace"]) *
      (1 + exp(coef(trad_death_satu_CDEmodl)["M_smoking"] +
                 coef(trad_death_satu_CDEmodl)["A0_ace:M_smoking"] +
                 coef(trad_mediator_satu)["(Intercept)"] +
                 coef(trad_mediator_satu)["A0_ace"])) /
      (1 + exp(coef(trad_death_satu_CDEmodl)["M_smoking"] +
                 coef(trad_mediator_satu)["(Intercept)"] + coef(trad_mediator_satu)["A0_ace"]))

    satu_PNIE[j] <- (1 + exp(coef(trad_mediator_satu)["(Intercept)"])) *
      (1 + exp(coef(trad_death_satu_CDEmodl)["M_smoking"] +
                 coef(trad_mediator_satu)["(Intercept)"] +
                 coef(trad_mediator_satu)["A0_ace"])) /
      ((1 + exp(coef(trad_mediator_satu)["(Intercept)"] + coef(trad_mediator_satu)["A0_ace"])) *
         (1 + exp(coef(trad_death_satu_CDEmodl)["M_smoking"] +
                    coef(trad_mediator_satu)["(Intercept)"])))
}

### Effet total
boxplot(data.frame(ATE_biasedQ, ATE_correct))
abline(h=1.614433)
# on retrouve la limite habituelle entre un modèle biaisé mais précis (peu flexible)
# versus un modèle moins biaisé car très flexible, mais avec une grande variance

### Effet direct contrôlé
boxplot(data.frame(biased_CDE_M0_ATE, flexible_CDE_M0_ATE, satu_CDE_M0_ATE))
abline(h=1.588235)

boxplot(data.frame(biased_CDE_M1_ATE, flexible_CDE_M1_ATE, satu_CDE_M1_ATE))
abline(h=1.600601)

### Effet direct et indirects naturels
# PNDE
boxplot(data.frame(biased_PNDE, flexible_PNDE, satu_PNDE))
abline(h=1.583042)

# TNIE
boxplot(data.frame(biased_TNIE, flexible_TNIE, satu_TNIE))
abline(h=1.078278)

# TNDE
boxplot(data.frame(biased_TNDE, flexible_TNDE, satu_TNDE))
abline(h=1.582382)

# PNIE
boxplot(data.frame(biased_PNIE, flexible_PNIE, satu_PNIE))
abline(h=1.078728)
