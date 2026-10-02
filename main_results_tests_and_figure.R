#---------------------------------------
# Load libraries and prepare data
#---------------------------------------
source("data_cleaning.R")
# workspace hygiene
KEEP <- c("df_long", "KEEP")
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Predicion 1
# onsesided exact binomial
#---------------------------------------
df_tmp <- df_long %>%
  filter(
    parent_status == "Spectators" &
	alternatives == "A1_A7" &
	meritocracy == "no"
  )
  
binom.test(
  x = sum(df_tmp$choice == "A7", na.rm = TRUE),
  n = sum(!is.na(df_tmp$choice)),
  p = 0.5,
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Predicion 2 
# one-sided exact mcnemar, Pos vs Neg
#---------------------------------------
df_tmp <- df_long %>%
  filter(
    parent_status == "Spectators",
    alternatives == "A1_A7",
    meritocracy %in% c("positive", "negative")
  ) %>%
    pivot_wider(
    id_cols = participant_code,
    names_from = meritocracy,
    values_from = choice
  )

table(df_tmp$positive, df_tmp$negative)
positive_switch <- table(df_tmp$positive, df_tmp$negative)[2]
negative_switch <- table(df_tmp$positive, df_tmp$negative)[3]
binom.test(
  x = positive_switch,
  n = positive_switch + negative_switch,
  p = 0.5,
  alternative = "greater"
)

#Below double check with exact2x2
#library(exact2x2)
#exact2x2(table(df_tmp$negative, df_tmp$positive), paired = TRUE, alternative = "greater")
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Predicion 2 
# Cochran-Armitage trend statistic,
# exact null by within-participant permutation
#---------------------------------------
#preparethe data
trend_mat <- df_long %>%
  filter(parent_status == "Spectators",
         alternatives   == "A1_A7") %>%
  mutate(fm = as.integer(choice == "A7")) %>%
  pivot_wider(id_cols     = participant_code,
              names_from  = meritocracy,
              values_from = fm) %>%
  drop_na(negative, no, positive) %>%
  select(negative, no, positive) %>%     # <- order defines the trend direction
  as.matrix()

#Linear scores (equivalent to -1,0,1)
scores <- c(1, 2, 3)                        # equally spaced, as in prop.trend.test
 
n_j  <- rep(nrow(trend_mat), 3)             # observations per condition
x_j  <- colSums(trend_mat)                  # A7 choices per condition
N    <- sum(n_j)
pbar <- sum(x_j) / N
 
U <- sum(scores * (x_j - n_j * pbar))       # CA numerator
U <- round(U)                               # exact integer; pbar division leaves fp dust
 
# equally spaced scores with unit spacing => the middle condition cancels and
# U reduces to sum_i (x_i,Pos - x_i,Neg); this is the object we permute
stopifnot(
  all.equal(diff(scores), c(1, 1)),
  U == sum(trend_mat %*% c(-1, 0, 1))
)
 
m <- sum(rowSums(trend_mat) %in% c(1, 2))   # participants informative under permutation
stopifnot(m > 0)
 
support <- -m:m
pmf <- 1
for (i in seq_len(m)) {
  new <- numeric(length(pmf) + 2)
  for (s in 0:2) new[s + seq_along(pmf)] <- new[s + seq_along(pmf)] + pmf / 3
  pmf <- new
}
 
p_exact <- sum(pmf[support >= U])           # H1: trend in predicted direction
 
c(N = nrow(trend_mat), m = m, U = U, p = p_exact)
colMeans(trend_mat)
 
# classical (independence-based) p-value, for the comparison sentence
#Var_U <- pbar * (1 - pbar) * (sum(n_j * scores^2) - sum(n_j * scores)^2 / N)
#1 - pnorm(U / sqrt(Var_U))
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 2
# one-sided exact mcnemar, Pos vs No
#---------------------------------------
df_tmp <- df_long %>%
  filter(
    parent_status == "Spectators",
    alternatives == "A1_A7",
    meritocracy %in% c("positive", "no")
  ) %>%
  pivot_wider(
    id_cols = participant_code,
    names_from = meritocracy,
    values_from = choice
  )

table(df_tmp$positive, df_tmp$no)
positive_switch <- table(df_tmp$positive, df_tmp$no)[2]
no_switch       <- table(df_tmp$positive, df_tmp$no)[3]
binom.test(
  x = positive_switch,
  n = positive_switch + no_switch,
  p = 0.5,
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 2
# one-sided exact mcnemar, No vs Neg
#---------------------------------------
df_tmp <- df_long %>%
  filter(
    parent_status == "Spectators",
    alternatives == "A1_A7",
    meritocracy %in% c("no", "negative")
  ) %>%
  pivot_wider(
    id_cols = participant_code,
    names_from = meritocracy,
    values_from = choice
  )

table(df_tmp$no, df_tmp$negative)
no_switch       <- table(df_tmp$no, df_tmp$negative)[2]
negative_switch <- table(df_tmp$no, df_tmp$negative)[3]
binom.test(
  x = no_switch,
  n = no_switch + negative_switch,
  p = 0.5,
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 3
# Fisher exact, Low vs High status, No meritocracy
#---------------------------------------
df_tmp <- df_long %>%
  filter(
    parent_status %in% c("Low", "High"),
    alternatives == "A1_A7",
    meritocracy == "no",
    !is.na(choice)
  ) %>%
  mutate(
    parent_status = factor(parent_status, levels = c("Low", "High")),
    choice        = factor(choice,        levels = c("A7", "A1"))
  )

table(df_tmp$parent_status, df_tmp$choice)
fisher.test(
  table(df_tmp$parent_status, df_tmp$choice),
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 3
# Fisher exact, Low vs High status, Pos meritocracy
#---------------------------------------
df_tmp <- df_long %>%
  filter(
    parent_status %in% c("Low", "High"),
    alternatives == "A1_A7",
    meritocracy == "positive",
    !is.na(choice)
  ) %>%
  mutate(
    parent_status = factor(parent_status, levels = c("Low", "High")),
    choice        = factor(choice,        levels = c("A7", "A1"))
  )

table(df_tmp$parent_status, df_tmp$choice)
fisher.test(
  table(df_tmp$parent_status, df_tmp$choice),
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 3
# Fisher exact, Low vs High status, Neg meritocracy
#---------------------------------------
df_tmp <- df_long %>%
  filter(
    parent_status %in% c("Low", "High"),
    alternatives == "A1_A7",
    meritocracy == "negative",
    !is.na(choice)
  ) %>%
  mutate(
    parent_status = factor(parent_status, levels = c("Low", "High")),
    choice        = factor(choice,        levels = c("A7", "A1"))
  )

table(df_tmp$parent_status, df_tmp$choice)
fisher.test(
  table(df_tmp$parent_status, df_tmp$choice),
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 3
# Cochran-Armitage trend across status, No meritocracy
# exact null by between-subject permutation of the status labels
#
# Status is coded High < Medium < Low so that the predicted direction is
# INCREASING, matching the "greater" convention used in the blocks above.
#---------------------------------------
df_tmp <- df_long %>%
  filter(
    parent_status %in% c("Low", "Medium", "High"),
    alternatives == "A1_A7",
    meritocracy == "no",
    !is.na(choice)
  ) %>%
  mutate(
    parent_status = factor(parent_status, levels = c("High", "Medium", "Low")),
    fm            = as.integer(choice == "A7")
  )

stopifnot(!any(duplicated(df_tmp$participant_code)))

table(df_tmp$parent_status, df_tmp$fm)

scores <- c(1, 2, 3)                                # High, Medium, Low
n_j    <- as.vector(table(df_tmp$parent_status))    # group sizes
x_j    <- as.vector(tapply(df_tmp$fm, df_tmp$parent_status, sum))
N      <- sum(n_j)
x_tot  <- sum(x_j)
pbar   <- x_tot / N
sbar   <- sum(n_j * scores) / N

U <- sum(scores * (x_j - n_j * pbar))               # CA numerator

# equal group sizes => sbar = 2, centred scores are (-1, 0, 1), and the Medium
# group cancels: U reduces to x_Low - x_High
stopifnot(
  all.equal(diff(scores), c(1, 1)),
  all(n_j == n_j[1]),
  isTRUE(all.equal(U, x_j[3] - x_j[1]))
)
U <- round(U)                                       # exact integer once balanced

# exact null: conditional on the margins, (x_high, x_medium, x_low) is multivariate
# It as an urn. There are 144 participants partitioned into three groups of 48. Draw 85 of them without replacement to be the A7-choosers. Under the null the count landing in each group has hypergeometric distribution. Enumerate every configuration compute U(x_h,x_m,x_l) and accumulate the upper tail when it is higher than the observed one.
p_exact <- 0
p_total <- 0
for (a in 0:min(n_j[1], x_tot)) {                        # x_High: can't exceed group size or total
  for (b in 0:min(n_j[2], x_tot - a)) {                  # x_Med: can't exceed what's left
    cc <- x_tot - a - b                                  # x_Low is determined
    if (cc < 0 || cc > n_j[3]) next                      # skip when x_Low is unfeasible
    prob <- exp(lchoose(n_j[1], a) + lchoose(n_j[2], b) +
                lchoose(n_j[3], cc) - lchoose(N, x_tot)) # P(of this configuration under the null)
    U_ab <- sum(scores * (c(a, b, cc) - n_j * pbar))     # Compute U of the configuartion
    p_total <- p_total + prob
    if (U_ab >= U - 1e-9) p_exact <- p_exact + prob      # accumulate the upper tail
  }
}
stopifnot(abs(p_total - 1) < 1e-9)

c(N = N, U = U, p = p_exact)
x_j / n_j                                           # observed A7 share by status

# classical (independence-based) p-value, for the comparison sentence
#Var_U <- pbar * (1 - pbar) * sum(n_j * (scores - sbar)^2)
#1 - pnorm(U / sqrt(Var_U))
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 4
# one-sided exact mcnemar, Pos vs Neg (optimistic, non-spectators)
#
# Fiest classify in types by first-order stochastic dominance against uniform beliefs
#
# b >_FOSD b_u  <=>  F(k) <= F_u(k) for all k, strict somewhere  -> optimistic
# b <_FOSD b_u  <=>  F(k) >= F_u(k) for all k, strict somewhere  -> pessimistic
# otherwise the CDFs cross (or coincide)                         -> ambiguous
#---------------------------------------
df_tmp <- df_long %>%
  mutate(
    belief_sum = belief_bottom + belief_middle + belief_top,
    F1         = belief_bottom / belief_sum,
    F2         = (belief_bottom + belief_middle) / belief_sum,
    belief_type = case_when(
      is.na(F1) | is.na(F2)                                        ~ NA_character_,
      F1 <= 1/3 + 1e-9 & F2 <= 2/3 + 1e-9 &
        (F1 < 1/3 - 1e-9 | F2 < 2/3 - 1e-9)                        ~ "optimistic",
      F1 >= 1/3 - 1e-9 & F2 >= 2/3 - 1e-9 &
        (F1 > 1/3 + 1e-9 | F2 > 2/3 + 1e-9)                        ~ "pessimistic",
      TRUE                                                         ~ "ambiguous"
    ),
    belief_type = factor(belief_type,
                         levels = c("pessimistic", "ambiguous", "optimistic"))
  ) %>%
  filter(
    parent_status %in% c("Low", "Medium", "High"),
	belief_type == "optimistic",
    alternatives == "A1_A7",
    meritocracy %in% c("positive", "negative")
  ) %>%
  pivot_wider(
    id_cols = participant_code,
    names_from = meritocracy,
    values_from = choice
  )

table(df_tmp$positive, df_tmp$negative)
positive_switch <- table(df_tmp$positive, df_tmp$negative)[2]
negative_switch <- table(df_tmp$positive, df_tmp$negative)[3]
binom.test(
  x = positive_switch,
  n = positive_switch + negative_switch,
  p = 0.5,
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 4
# one-sided exact mcnemar, Pos vs No (optimistic, non-spectators)
#
# First classify in types by first-order stochastic dominance against uniform beliefs
#
# b >_FOSD b_u  <=>  F(k) <= F_u(k) for all k, strict somewhere  -> optimistic
# b <_FOSD b_u  <=>  F(k) >= F_u(k) for all k, strict somewhere  -> pessimistic
# otherwise the CDFs cross (or coincide)                         -> ambiguous
#---------------------------------------
df_tmp <- df_long %>%
  mutate(
    belief_sum = belief_bottom + belief_middle + belief_top,
    F1         = belief_bottom / belief_sum,
    F2         = (belief_bottom + belief_middle) / belief_sum,
    belief_type = case_when(
      is.na(F1) | is.na(F2)                                        ~ NA_character_,
      F1 <= 1/3 + 1e-9 & F2 <= 2/3 + 1e-9 &
        (F1 < 1/3 - 1e-9 | F2 < 2/3 - 1e-9)                        ~ "optimistic",
      F1 >= 1/3 - 1e-9 & F2 >= 2/3 - 1e-9 &
        (F1 > 1/3 + 1e-9 | F2 > 2/3 + 1e-9)                        ~ "pessimistic",
      TRUE                                                         ~ "ambiguous"
    ),
    belief_type = factor(belief_type,
                         levels = c("pessimistic", "ambiguous", "optimistic"))
  ) %>%
  filter(
    parent_status %in% c("Low", "Medium", "High"),
	belief_type == "optimistic",
    alternatives == "A1_A7",
    meritocracy %in% c("positive", "no")
  ) %>%
  pivot_wider(
    id_cols = participant_code,
    names_from = meritocracy,
    values_from = choice
  )

table(df_tmp$positive, df_tmp$no)
positive_switch <- table(df_tmp$positive, df_tmp$no)[2]
no_switch       <- table(df_tmp$positive, df_tmp$no)[3]
binom.test(
  x = positive_switch,
  n = positive_switch + no_switch,
  p = 0.5,
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 4
# one-sided exact mcnemar, No vs Neg (optimistic, non-spectators)
#
# Fiest classify in types by first-order stochastic dominance against uniform beliefs
#
# b >_FOSD b_u  <=>  F(k) <= F_u(k) for all k, strict somewhere  -> optimistic
# b <_FOSD b_u  <=>  F(k) >= F_u(k) for all k, strict somewhere  -> pessimistic
# otherwise the CDFs cross (or coincide)                         -> ambiguous
#---------------------------------------
df_tmp <- df_long %>%
  mutate(
    belief_sum = belief_bottom + belief_middle + belief_top,
    F1         = belief_bottom / belief_sum,
    F2         = (belief_bottom + belief_middle) / belief_sum,
    belief_type = case_when(
      is.na(F1) | is.na(F2)                                        ~ NA_character_,
      F1 <= 1/3 + 1e-9 & F2 <= 2/3 + 1e-9 &
        (F1 < 1/3 - 1e-9 | F2 < 2/3 - 1e-9)                        ~ "optimistic",
      F1 >= 1/3 - 1e-9 & F2 >= 2/3 - 1e-9 &
        (F1 > 1/3 + 1e-9 | F2 > 2/3 + 1e-9)                        ~ "pessimistic",
      TRUE                                                         ~ "ambiguous"
    ),
    belief_type = factor(belief_type,
                         levels = c("pessimistic", "ambiguous", "optimistic"))
  ) %>%
  filter(
    parent_status %in% c("Low", "Medium", "High"),
	belief_type == "optimistic",
    alternatives == "A1_A7",
    meritocracy %in% c("no", "negative")
  ) %>%
  pivot_wider(
    id_cols = participant_code,
    names_from = meritocracy,
    values_from = choice
  )

table(df_tmp$no, df_tmp$negative)
no_switch       <- table(df_tmp$no, df_tmp$negative)[2]
negative_switch <- table(df_tmp$no, df_tmp$negative)[3]
binom.test(
  x = no_switch,
  n = no_switch + negative_switch,
  p = 0.5,
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 4
# Cochran-Armitage trend statistic (optimistic, non-spectators),
# exact null by within-participant permutation
#---------------------------------------
#prepare the data
trend_mat <- df_long %>%
  mutate(
    belief_sum = belief_bottom + belief_middle + belief_top,
    F1         = belief_bottom / belief_sum,
    F2         = (belief_bottom + belief_middle) / belief_sum,
    belief_type = case_when(
      is.na(F1) | is.na(F2)                                        ~ NA_character_,
      F1 <= 1/3 + 1e-9 & F2 <= 2/3 + 1e-9 &
        (F1 < 1/3 - 1e-9 | F2 < 2/3 - 1e-9)                        ~ "optimistic",
      F1 >= 1/3 - 1e-9 & F2 >= 2/3 - 1e-9 &
        (F1 > 1/3 + 1e-9 | F2 > 2/3 + 1e-9)                        ~ "pessimistic",
      TRUE                                                         ~ "ambiguous"
    ),
    belief_type = factor(belief_type,
                         levels = c("pessimistic", "ambiguous", "optimistic"))
  ) %>%
  filter(parent_status %in% c("Low", "Medium", "High"),
         belief_type    == "optimistic",
         alternatives   == "A1_A7") %>%
  mutate(fm = as.integer(choice == "A7")) %>%
  pivot_wider(id_cols     = participant_code,
              names_from  = meritocracy,
              values_from = fm) %>%
  drop_na(negative, no, positive) %>%
  select(negative, no, positive) %>%     # <- order defines the trend direction
  as.matrix()

#Linear scores (equivalent to -1,0,1)
scores <- c(1, 2, 3)                        # equally spaced, as in prop.trend.test

n_j  <- rep(nrow(trend_mat), 3)             # observations per condition
x_j  <- colSums(trend_mat)                  # A7 choices per condition
N    <- sum(n_j)
pbar <- sum(x_j) / N

U <- sum(scores * (x_j - n_j * pbar))       # CA numerator
U <- round(U)                               # exact integer; pbar division leaves fp dust

# equally spaced scores with unit spacing => the middle condition cancels and
# U reduces to sum_i (x_i,Pos - x_i,Neg); this is the object we permute
stopifnot(
  all.equal(diff(scores), c(1, 1)),
  U == sum(trend_mat %*% c(-1, 0, 1))
)

m <- sum(rowSums(trend_mat) %in% c(1, 2))   # participants informative under permutation
stopifnot(m > 0)

support <- -m:m
pmf <- 1
for (i in seq_len(m)) {
  new <- numeric(length(pmf) + 2)
  for (s in 0:2) new[s + seq_along(pmf)] <- new[s + seq_along(pmf)] + pmf / 3
  pmf <- new
}

p_exact <- sum(pmf[support >= U])           # H1: trend in predicted direction

c(N = nrow(trend_mat), m = m, U = U, p = p_exact)
colMeans(trend_mat)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 4 (robustness)
# one-sided exact mcnemar, Pos vs Neg (optimistic, non-spectators)
#
# Alternative classification: optimistic if the expected performance level is
# strictly positive, scoring bottom = -1, middle = 0, top = +1. With symmetric
# scores this reduces to belief_top > belief_bottom.
#---------------------------------------
df_tmp <- df_long %>%
  mutate(
    belief_sum  = belief_bottom + belief_middle + belief_top,
    belief_ev   = (-1 * belief_bottom + 0 * belief_middle + 1 * belief_top) / belief_sum,
    belief_type = case_when(
      is.na(belief_ev)   ~ NA_character_,
      belief_ev >  1e-9  ~ "optimistic",
      belief_ev < -1e-9  ~ "pessimistic",
      TRUE               ~ "neutral"
    ),
    belief_type = factor(belief_type,
                         levels = c("pessimistic", "neutral", "optimistic"))
  ) %>%
  filter(
    parent_status %in% c("Low", "Medium", "High"),
    belief_type == "optimistic",
    alternatives == "A1_A7",
    meritocracy %in% c("positive", "negative")
  ) %>%
  pivot_wider(
    id_cols = participant_code,
    names_from = meritocracy,
    values_from = choice
  )

table(df_tmp$positive, df_tmp$negative)
positive_switch <- table(df_tmp$positive, df_tmp$negative)[2]
negative_switch <- table(df_tmp$positive, df_tmp$negative)[3]
binom.test(
  x = positive_switch,
  n = positive_switch + negative_switch,
  p = 0.5,
  alternative = "greater"
)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Test Prediction 4 (robustness)
# Cochran-Armitage trend statistic (optimistic, non-spectators),
# exact null by within-participant permutation
#
# Alternative classification: optimistic if the expected performance level is
# strictly positive, scoring bottom = -1, middle = 0, top = +1.
#---------------------------------------
#prepare the data
trend_mat <- df_long %>%
  mutate(
    belief_sum  = belief_bottom + belief_middle + belief_top,
    belief_ev   = (-1 * belief_bottom + 0 * belief_middle + 1 * belief_top) / belief_sum,
    belief_type = case_when(
      is.na(belief_ev)   ~ NA_character_,
      belief_ev >  1e-9  ~ "optimistic",
      belief_ev < -1e-9  ~ "pessimistic",
      TRUE               ~ "neutral"
    ),
    belief_type = factor(belief_type,
                         levels = c("pessimistic", "neutral", "optimistic"))
  ) %>%
  filter(parent_status %in% c("Low", "Medium", "High"),
         belief_type    == "optimistic",
         alternatives   == "A1_A7") %>%
  mutate(fm = as.integer(choice == "A7")) %>%
  pivot_wider(id_cols     = participant_code,
              names_from  = meritocracy,
              values_from = fm) %>%
  drop_na(negative, no, positive) %>%
  select(negative, no, positive) %>%     # <- order defines the trend direction
  as.matrix()

#Linear scores (equivalent to -1,0,1)
scores <- c(1, 2, 3)                        # equally spaced, as in prop.trend.test

n_j  <- rep(nrow(trend_mat), 3)             # observations per condition
x_j  <- colSums(trend_mat)                  # A7 choices per condition
N    <- sum(n_j)
pbar <- sum(x_j) / N

U <- sum(scores * (x_j - n_j * pbar))       # CA numerator
U <- round(U)                               # exact integer; pbar division leaves fp dust

# equally spaced scores with unit spacing => the middle condition cancels and
# U reduces to sum_i (x_i,Pos - x_i,Neg); this is the object we permute
stopifnot(
  all.equal(diff(scores), c(1, 1)),
  U == sum(trend_mat %*% c(-1, 0, 1))
)

m <- sum(rowSums(trend_mat) %in% c(1, 2))   # participants informative under permutation
stopifnot(m > 0)

support <- -m:m
pmf <- 1
for (i in seq_len(m)) {
  new <- numeric(length(pmf) + 2)
  for (s in 0:2) new[s + seq_along(pmf)] <- new[s + seq_along(pmf)] + pmf / 3
  pmf <- new
}

p_exact <- sum(pmf[support >= U])           # H1: trend in predicted direction

c(N = nrow(trend_mat), m = m, U = U, p = p_exact)
colMeans(trend_mat)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Figure: main results (Predictions 1-4)
#
# All shares and confidence intervals are computed from df_long.
# Intervals are standard symmetric (Wald) 95% intervals for a proportion:
#
#   se    = sqrt( share * (1 - share) / n )
#   lower = share - qnorm(0.975) * se
#   upper = share + qnorm(0.975) * se
#
# with qnorm(0.975) = 1.959964. These are the usual normal-approximation intervals
#---------------------------------------
library(ggplot2)
update_geom_defaults("text", list(family = "serif")) 
 
#--- panel 1: spectators, by meritocracy (Predictions 1 & 2) -----------------
p1 <- df_long %>%
  filter(parent_status == "Spectators",
         alternatives  == "A1_A7",
         !is.na(choice)) %>%
  group_by(meritocracy) %>%
  summarise(x = sum(choice == "A7"), n = n(), .groups = "drop") %>%
  mutate(panel = "Results 1 & 2",
         xkey  = paste0("1_", meritocracy),
         xlab  = "Spectator",
         merit = as.character(meritocracy))
 
#--- panel 2: parental status, No meritocracy (Prediction 3) -----------------
p2 <- df_long %>%
  filter(parent_status %in% c("Low", "Medium", "High"),
         alternatives  == "A1_A7",
         meritocracy   == "no",
         !is.na(choice)) %>%
  group_by(parent_status) %>%
  summarise(x = sum(choice == "A7"), n = n(), .groups = "drop") %>%
  mutate(panel = "Result 3",
         xkey  = paste0("2_", parent_status),
         xlab  = as.character(parent_status),
         merit = "no")
 
#--- panel 3: optimistic non-spectators, by meritocracy (Prediction 4) -------
p3 <- df_long %>%
  mutate(
    belief_sum = belief_bottom + belief_middle + belief_top,
    F1         = belief_bottom / belief_sum,
    F2         = (belief_bottom + belief_middle) / belief_sum,
    belief_type = case_when(
      is.na(F1) | is.na(F2)                                        ~ NA_character_,
      F1 <= 1/3 + 1e-9 & F2 <= 2/3 + 1e-9 &
        (F1 < 1/3 - 1e-9 | F2 < 2/3 - 1e-9)                        ~ "optimistic",
      F1 >= 1/3 - 1e-9 & F2 >= 2/3 - 1e-9 &
        (F1 > 1/3 + 1e-9 | F2 > 2/3 + 1e-9)                        ~ "pessimistic",
      TRUE                                                         ~ "ambiguous"
    )
  ) %>%
  filter(parent_status %in% c("Low", "Medium", "High"),
         belief_type   == "optimistic",
         alternatives  == "A1_A7",
         !is.na(choice)) %>%
  group_by(meritocracy) %>%
  summarise(x = sum(choice == "A7"), n = n(), .groups = "drop") %>%
  mutate(panel = "Result 4",
         xkey  = paste0("3_", meritocracy),
         xlab  = "Member\n(Optimists)",
         merit = as.character(meritocracy))
 
#--- combine and add exact Clopper-Pearson intervals -------------------------
df_fig <- bind_rows(p1, p2, p3) %>%
  mutate(
    share = x / n,
    se    = sqrt(share * (1 - share) / n),
    lo    = pmax(0, share - qnorm(0.975) * se),
    hi    = pmin(1, share + qnorm(0.975) * se),
    # exact Clopper-Pearson alternative:
    # lo  = ifelse(x == 0, 0, qbeta(0.025, x,     n - x + 1)),
    # hi  = ifelse(x == n, 1, qbeta(0.975, x + 1, n - x)),
    panel = factor(panel, levels = c("Results 1 & 2", "Result 3", "Result 4")),
    xkey  = factor(xkey, levels = c("1_negative", "1_no", "1_positive",
                                    "2_Low", "2_Medium", "2_High",
                                    "3_negative", "3_no", "3_positive")),
    merit = factor(merit, levels = c("negative", "no", "positive"),
                   labels = c("Negative", "No", "Positive"))
  )
 
#df_fig    # inspect before plotting
 
#--- plot --------------------------------------------------------------------
# Sizing: the figure is saved with these sizes
FIG_W <- 8.5    # inches
FIG_H <- 4      # inches
 
pt <- function(x) x / .pt   # ggplot2 text `size` is in mm; .pt converts from pt
 
cols <- c("Negative" = "#D55E00",   # vermillion
          "No"       = "#0072B2",   # blue
          "Positive" = "#009E73")   # green
 
xlabs <- setNames(df_fig$xlab, as.character(df_fig$xkey))   # tick = parental status
 
ggplot(df_fig, aes(x = xkey, y = share, colour = merit)) +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             colour = "grey45", linewidth = 0.5) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 1, alpha = 0.85) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * share)),
            hjust = -0.20, size = pt(9), fontface = "bold",
            show.legend = FALSE) +
  geom_text(aes(y = 0.02, label = sprintf("(n = %d)", n)),
            colour = "grey25", size = pt(8), vjust = 0, show.legend = FALSE) +
  facet_wrap(~ panel, scales = "free_x") +
  scale_colour_manual(values = cols, name = "Meritocracy treatment") +
  scale_y_continuous(limits = c(0, 1.03),
                     breaks = seq(0, 1, 0.1),
                     labels = function(v) paste0(round(100 * v), "%"),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_x_discrete(labels = xlabs, expand = expansion(add = 0.8)) +
  labs(x = "Parental status",
       y = "Preference for mobility") +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    legend.margin      = margin(t = -2),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    strip.background   = element_rect(fill = "grey94", colour = NA),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9, lineheight = 0.95),
    axis.title.y       = element_text(size = 11, margin = margin(r = 6)),
    axis.title.x       = element_text(size = 11, margin = margin(t = 6)),
    plot.margin        = margin(4, 6, 2, 2)
  )
 
ggsave("FigureMainResults_R.pdf", width = FIG_W, height = FIG_H)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Figure v2: distribution of motives (no tremble vs with tremble)
# Two panels: the four model-based types, and the two aggregate shares.
#
# Estimates and standard errors are read from the two results files produced by
# 02_stima_no_tremble_se_bestmodel.R and 02_stima_tremble_se_bestmodel.R.
#
# The `se` column in those files is a bootstrap STANDARD ERROR, not a CI: it is
# apply(boot_results, 2, sd, na.rm = TRUE) over 1000 bootstrap replications in
# which participants are resampled with replacement within each role.
# Intervals below are therefore the usual symmetric normal-approximation ones,
#
#   lower = estimate - qnorm(0.975) * se
#   upper = estimate + qnorm(0.975) * se
#
# The aggregate model-based share M = pi_Self + pi_Mob + pi_Mer + pi_MobMer
# equals 1 - pi_Rand exactly, so its standard error is the one estimated for
# pi_Rand. No extra computation is needed.
#---------------------------------------
library(ggplot2)
library(dplyr)
library(readr)

FOLD <- "Clean version - no overlap - 5 types"

nt <- read_csv(file.path(FOLD, "results_no_tremble_se_bestmodel.csv"),
               show_col_types = FALSE) %>%
  mutate(estimate = as.numeric(estimate), se = as.numeric(se), spec = "No tremble")

tr <- read_csv(file.path(FOLD, "results_tremble_se_bestmodel.csv"),
               show_col_types = FALSE) %>%
  mutate(estimate = as.numeric(estimate), se = as.numeric(se), spec = "With tremble")

#--- the five estimated type shares ------------------------------------------
types <- bind_rows(nt, tr) %>%
  filter(parameter %in% c("pi_1_self", "pi_2_self_mob", "pi_3_self_mer",
                          "pi_4_self_mob_mer", "pi_5_rand")) %>%
  mutate(type = recode(parameter,
                       pi_1_self         = "Self",
                       pi_2_self_mob     = "Mob",
                       pi_3_self_mer     = "Mer",
                       pi_4_self_mob_mer = "MobMer",
                       pi_5_rand         = "Rand"))

#--- aggregate model-based share: M = 1 - pi_Rand, so se(M) = se(pi_Rand) -----
agg <- types %>%
  filter(type == "Rand") %>%
  transmute(spec, type = "M", estimate = 1 - estimate, se = se)

df_fig <- bind_rows(types %>% select(spec, type, estimate, se), agg) %>%
  mutate(
    lo   = pmax(0, estimate - qnorm(0.975) * se),
    hi   = pmin(1, estimate + qnorm(0.975) * se),
    type  = factor(type, levels = c("Self", "Mob", "Mer", "MobMer", "M", "Rand")),
    spec  = factor(spec, levels = c("No tremble", "With tremble")),
    panel = factor(ifelse(type %in% c("M", "Rand"), "Aggregate shares",
                                                    "Model-based types"),
                   levels = c("Model-based types", "Aggregate shares"))
  ) %>%
  arrange(type, spec)

df_fig    # inspect before plotting

#--- plot --------------------------------------------------------------------
# Dimensions and all font sizes match FigureMainResults_R so the two figures
# are typographically consistent in the code. Note the pages differ: the main
# figure is placed at ~6.9in (scale 0.81) while this one sits on a landscape
# page where adjustbox clamps at ~9in, so an 8.5in save is placed at scale 1.
# To match the RENDERED size too, save this one at 8.5 / 0.81 = 10.5 wide and
# 4 / 0.81 = 4.9 high, and let adjustbox shrink it back.
FIG_W <- 8.5    # inches  (same as FigureMainResults_R)
FIG_H <- 4      # inches  (same as FigureMainResults_R)

pt <- function(x) x / .pt   # ggplot2 text `size` is in mm; .pt converts from pt

cols <- c("No tremble"   = "#0072B2",   # blue
          "With tremble" = "#D55E00")   # vermillion

# Two-line tick labels, as plotmath so the Greek renders on the base pdf()
# device (which is WinAnsi-encoded and cannot show literal Greek characters).
# NOTE: supplied in FACTOR-LEVEL order -- Self, Mob, Mer, MobMer, M, Rand.
xlabs <- list(
  Self   = expression(atop("Selfish",     paste("(", delta, " = 0, ", rho, " = 0)"))),
  Mob    = expression(atop("Mobility",    paste("(", delta, " > 0, ", rho, " = 0)"))),
  Mer    = expression(atop("Meritocracy", paste("(", delta, " = 0, ", rho, " > 0)"))),
  MobMer = expression(atop("Mob & Mer",   paste("(", delta, " > 0, ", rho, " > 0)"))),
  M      = expression(atop("All model-based", group("(", italic(k) %in% italic(M), ")"))),
  Rand   = expression(atop("Random",      "(indifferent)"))
)

# Supplied as a FUNCTION: with free x scales each panel asks only for the
# levels it contains, so positional matching would misalign the labels.
xlab_fun <- function(breaks) unlist(xlabs[as.character(breaks)], use.names = FALSE)

dodge <- position_dodge(width = 0.75)

ggplot(df_fig, aes(x = type, y = estimate, colour = spec, group = spec)) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 1,
                 position = dodge, alpha = 0.85) +
  geom_point(size = 2.6, position = dodge) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * estimate)),
            position = dodge, hjust = -0.20, size = pt(9),
            fontface = "bold", family = "serif", show.legend = FALSE) +
  scale_colour_manual(values = cols, name = "Specification") +
  scale_y_continuous(limits = c(-0.05, 1.03),   # headroom below 0; ticks unchanged
                     breaks = seq(0, 1, 0.1),
                     labels = function(v) paste0(round(100 * v), "%"),
                     expand = expansion(mult = c(0.01, 0.02))) +
  facet_grid(~ panel, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = xlab_fun, expand = expansion(add = 0.8)) +
  labs(x = "Type",
       y = expression(paste("Estimated type probability  ", pi[k]))) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    legend.margin      = margin(t = -2),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    strip.background   = element_rect(fill = "grey94", colour = NA),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9),    # as in FigureMainResults_R;
                                                    # atop() renders smaller, raise to 10-11 if needed
    axis.title.y       = element_text(size = 11, margin = margin(r = 6)),
    axis.title.x       = element_text(size = 11, margin = margin(t = 6)),
    plot.margin        = margin(4, 6, 2, 2)
  )

ggsave("FigureStructural_R_v2.pdf", width = FIG_W, height = FIG_H)   # base pdf(): no cairo on Windows
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Beliefs and self-serving meritocracy
# Estimates equations (reg1) and (reg2) and prints the LaTeX table.
#
# Dependent variable: the empirical posterior probability of having a
# meritocracy motive, alone or combined with a mobility motive, taken from the
# model WITH tremble and expressed in percentage points:
#
#   mer = 100 * (type_3 + type_4)
#
# Scaling by 100 leaves R-squared and all t-statistics unchanged; every
# coefficient and standard error is simply 100 times its value on the 0-1
# scale, which makes the table readable without exponents.
#
# Stata's `robust` for OLS is HC1, so vcovHC(type = "HC1") reproduces it.
#---------------------------------------
library(dplyr)
library(readr)
library(sandwich)
library(lmtest)

FOLD <- "Clean version - no overlap - 5 types"

#--- posterior type probabilities (model with tremble) -----------------------
post <- bind_rows(
  read_csv(file.path(FOLD, "posterior_spec_tr_bestmodel.csv"), show_col_types = FALSE) %>%
    rename(participant_code = id_spec)   %>% mutate(role_file = "Spectators"),
  read_csv(file.path(FOLD, "posterior_low_tr_bestmodel.csv"),  show_col_types = FALSE) %>%
    rename(participant_code = id_low)    %>% mutate(role_file = "Low"),
  read_csv(file.path(FOLD, "posterior_med_tr_bestmodel.csv"),  show_col_types = FALSE) %>%
    rename(participant_code = id_medium) %>% mutate(role_file = "Medium"),
  read_csv(file.path(FOLD, "posterior_high_tr_bestmodel.csv"), show_col_types = FALSE) %>%
    rename(participant_code = id_high)   %>% mutate(role_file = "High")
) %>%
  mutate(mer = 100 * (type_3 + type_4))

#--- LABEL CHECK: does each posterior file hold the role it is named for? ----
# This must be a diagonal 4x4 table with 48 on each diagonal cell. Anything
# off-diagonal means the posterior files and df_long disagree about who is who.
lab_chk <- df_long %>%
  distinct(participant_code, parent_status) %>%
  inner_join(post %>% select(participant_code, role_file), by = "participant_code") %>%
  mutate(role_file = factor(role_file, levels = c("Spectators", "High", "Medium", "Low")))

table(lab_chk$role_file, lab_chk$parent_status)

stopifnot(
  nrow(lab_chk) == 192,
  all(as.character(lab_chk$role_file) == as.character(lab_chk$parent_status))
)

#--- merge with beliefs and role ---------------------------------------------
df_reg <- df_long %>%
  distinct(participant_code, parent_status, belief_top, belief_middle, belief_bottom) %>%
  inner_join(post, by = "participant_code") %>%
  mutate(
    Beliefs = belief_top - belief_bottom,
    Members = as.integer(parent_status != "Spectators"),
    Low     = as.integer(parent_status == "Low"),
    Med     = as.integer(parent_status == "Medium"),
    High    = as.integer(parent_status == "High")
  )

stopifnot(
  nrow(df_reg) == 192,
  !any(duplicated(df_reg$participant_code)),
  all(df_reg$Low + df_reg$Med + df_reg$High == df_reg$Members)
)

# label check: means by role, and the saturated-model identity
#   mean(mer_g) = intercept_g + slope_g * mean(Beliefs_g)
df_reg %>%
  group_by(parent_status) %>%
  summarise(n = n(), mean_mer = mean(mer), mean_Beliefs = mean(Beliefs),
            .groups = "drop")

#--- estimation ---------------------------------------------------------------
m1 <- lm(mer ~ Beliefs + Members + Beliefs:Members, data = df_reg)
m2 <- lm(mer ~ Beliefs + Low + Med + High +
               Beliefs:Low + Beliefs:Med + Beliefs:High, data = df_reg)

ct1 <- coeftest(m1, vcov = vcovHC(m1, type = "HC1"))
ct2 <- coeftest(m2, vcov = vcovHC(m2, type = "HC1"))

ct1
ct2

#--- LaTeX table --------------------------------------------------------------
# three significant digits, as in the published table
fmt <- function(x) formatC(x, digits = 3, format = "fg", flag = "#")

tex_cell <- function(ct, term) {
  if (is.null(term) || !(term %in% rownames(ct))) return(c("", ""))
  e <- ct[term, 1]; s <- ct[term, 2]; p <- ct[term, 4]
  st <- if (p < 0.01) "\\sym{***}" else if (p < 0.05) "\\sym{**}" else
        if (p < 0.10) "\\sym{*}"   else ""
  c(paste0("$", fmt(e), st, "$"), paste0("$(", fmt(s), ")$"))
}

rows <- list(
  list(lab = "Beliefs [$b(r_H) - b(r_L)$]",                     t1 = "Beliefs",         t2 = "Beliefs"),
  list(lab = "Members:",                                        t1 = "Members",         t2 = NULL),
  list(lab = "\\hspace{1em} - Low parent's status",              t1 = NULL,              t2 = "Low"),
  list(lab = "\\hspace{1em} - Medium parent's status",           t1 = NULL,              t2 = "Med"),
  list(lab = "\\hspace{1em} - High parent's status",             t1 = NULL,              t2 = "High"),
  list(lab = "Beliefs $\\times$ Members:",                       t1 = "Beliefs:Members", t2 = NULL),
  list(lab = "\\hspace{1em} - Beliefs $\\times$ Low parent's status",            t1 = NULL,              t2 = "Beliefs:Low"),
  list(lab = "\\hspace{1em} - Beliefs $\\times$ Medium parent's status",         t1 = NULL,              t2 = "Beliefs:Med"),
  list(lab = "\\hspace{1em} - Beliefs $\\times$ High parent's status",           t1 = NULL,              t2 = "Beliefs:High"),
  list(lab = "Constant",                                        t1 = "(Intercept)",     t2 = "(Intercept)")
)

body <- character(0)
for (r in rows) {
  a <- tex_cell(ct1, r$t1)
  b <- tex_cell(ct2, r$t2)
  body <- c(body,
            sprintf("%s & %s & %s \\\\", r$lab, a[1], b[1]),
            sprintf(" & %s & %s \\\\",   a[2], b[2]))
}

cat(
"\\begin{table}[ht]\n",
"\\centering\n",
"\\begin{footnotesize}\n",
"\\setlength{\\tabcolsep}{7pt}\n",
"\\begin{tabular}{p{0.42\\linewidth}cc}\n",
"\\toprule\n",
" & \\multicolumn{2}{c}{Posterior that $\\rho_i > 0$ ($\\times 100$)} \\\\\n",
"\\cmidrule(lr){2-3}\n",
" & (1) & (2) \\\\\n",
" & Aggregate & Per status \\\\\n",
"\\midrule\n",
paste(body, collapse = "\n"), "\n",
"\\midrule\n",
sprintf("Observations & %d & %d \\\\\n", nobs(m1), nobs(m2)),
sprintf("R-squared & %.3f & %.3f \\\\\n",
        summary(m1)$r.squared, summary(m2)$r.squared),
"\\bottomrule\n",
"\\end{tabular}\n",
"\\end{footnotesize}\n",
"\n\\vspace{1em}\n",
"\\caption{Beliefs and Meritocracy Motive}\n",
"\\label{tab:beliefs_members_status}\n",
"\\vspace{1em}\n\n",
"\\begin{minipage}{0.95\\textwidth}\n",
"\\footnotesize\n",
"\\textit{Notes:} The dependent variable is the individual posterior probability\n",
"that participant $i$ has a meritocracy motive, $\\rho_i > 0$, i.e.\\ the sum of\n",
"the posteriors for the \\textit{Meritocracy} and \\textit{Mobility \\& Meritocracy}\n",
"types, the two types for which a meritocracy motive is necessary to rationalise\n",
"observed choices. Posteriors are computed from the structural model with\n",
"role-specific trembles reported in Figure~\\ref{fig: structural} and are\n",
"expressed in percentage points. All 192 participants enter both regressions;\n",
"no observation is excluded on the basis of its posterior value.\n",
"$\\text{Beliefs}$ is the difference between optimistic and pessimistic beliefs,\n",
"$b(r_H) - b(r_L)$, in percentage points;\n",
"$\\text{Members} = 1\\{\\text{participant} \\neq \\text{spectator}\\}$;\n",
"$\\text{Low} = 1\\{\\text{parent status = Low}\\}$, and likewise for the other\n",
"statuses. Spectators are the omitted category in both specifications.\n",
"Estimation is by OLS with heteroskedasticity-robust (HC1) standard errors,\n",
"reported in parentheses. \\sym{***} $p<0.01$, \\sym{**} $p<0.05$,\n",
"\\sym{*} $p<0.1$.\n",
"\\end{minipage}\n",
"\\end{table}\n",
sep = "")

#--- Wald test: are the three member slopes equal? ---------------------------
# In (2) the interactions are relative to spectators, so equality of the three
# interactions is equality of the three member slopes.
# H0: Beliefs:Low = Beliefs:Med = Beliefs:High   (2 restrictions)

# restricted model: status-specific intercepts, one common member slope
m3  <- lm(mer ~ Beliefs + Low + Med + High + Beliefs:Members, data = df_reg)
ct3 <- coeftest(m3, vcov = vcovHC(m3, type = "HC1"))
ct3

# robust Wald test of m3 (restricted) against m2 (unrestricted).
# vcov must be passed as a FUNCTION
car::linearHypothesis(m2,
  c("Beliefs:Low = Beliefs:Med", "Beliefs:Med = Beliefs:High"),
  vcov. = vcovHC(m2, type = "HC1"))

#--- Slope for members: beta_1 + beta_3 --------------------------------------
# The interaction beta_3 tests whether the belief slope DIFFERS between members
# and spectators. The member slope itself is beta_1 + beta_3, since spectators
# are the omitted category. The two answer different questions and need not
# agree: beta_3 is a difference that inherits the noise of the 48-spectator
# slope, whereas beta_1 + beta_3 is estimated off the 144 members.
#
# The standard error needs the covariance term. beta_1_hat and beta_3_hat are
# strongly negatively correlated, so sqrt(V11 + V33) alone would badly overstate
# it: the spectator sampling error enters beta_1 positively and beta_3
# negatively and cancels in the sum.
b1 <- coef(m1)
v1 <- vcovHC(m1, type = "HC1")

tot    <- b1["Beliefs"] + b1["Beliefs:Members"]
se_tot <- sqrt(v1["Beliefs", "Beliefs"] + v1["Beliefs:Members", "Beliefs:Members"] +
               2 * v1["Beliefs", "Beliefs:Members"])

c(slope_members = unname(tot), se = unname(se_tot),
  z = unname(tot / se_tot), p = unname(2 * pnorm(-abs(tot / se_tot))),
  lo = unname(tot - qnorm(0.975) * se_tot),
  hi = unname(tot + qnorm(0.975) * se_tot))

# same test via car, as a check on the delta-method arithmetic
car::linearHypothesis(m1, "Beliefs + Beliefs:Members = 0",
                      vcov. = vcovHC(m1, type = "HC1"))

# reparameterisation: members as the omitted category, so the member slope
# appears directly as a coefficient with its own standard error
df_reg$Spectator <- as.integer(df_reg$parent_status == "Spectators")
m1b  <- lm(mer ~ Beliefs + Spectator + Beliefs:Spectator, data = df_reg)
ct1b <- coeftest(m1b, vcov = vcovHC(m1b, type = "HC1"))
ct1b

# the "Beliefs" row of m1b must reproduce the delta-method figures above
stopifnot(isTRUE(all.equal(unname(ct1b["Beliefs", 1]), unname(tot))),
          isTRUE(all.equal(unname(ct1b["Beliefs", 2]), unname(se_tot))))

#--- bottom vs top quintile of Beliefs, members and spectators (model 1) ------
# Compares the AVERAGE belief in the bottom and top quintile of the pooled
# sample, so both roles are evaluated over the same range and the contrast
# reflects the slope difference rather than a difference in spread.
b1 <- coef(m1)
v1 <- vcovHC(m1, type = "HC1")

cut_lo <- quantile(df_reg$Beliefs, 0.20)
cut_hi <- quantile(df_reg$Beliefs, 0.80)
B_lo   <- mean(df_reg$Beliefs[df_reg$Beliefs <= cut_lo])
B_hi   <- mean(df_reg$Beliefs[df_reg$Beliefs >= cut_hi])
gap    <- B_hi - B_lo

c(cut20 = unname(cut_lo), cut80 = unname(cut_hi),
  mean_bottom = B_lo, mean_top = B_hi, gap = gap)

# slope, intercept and contrast for each role
sl_spec <- b1["Beliefs"]
se_spec <- sqrt(v1["Beliefs", "Beliefs"])
in_spec <- b1["(Intercept)"]

sl_memb <- b1["Beliefs"] + b1["Beliefs:Members"]
se_memb <- sqrt(v1["Beliefs", "Beliefs"] + v1["Beliefs:Members", "Beliefs:Members"] +
                2 * v1["Beliefs", "Beliefs:Members"])
in_memb <- b1["(Intercept)"] + b1["Members"]

sl_diff <- b1["Beliefs:Members"]
se_diff <- sqrt(v1["Beliefs:Members", "Beliefs:Members"])

out <- rbind(
  Spectators = c(pred_bottom = in_spec + sl_spec * B_lo,
                 pred_top    = in_spec + sl_spec * B_hi,
                 shift_pp    = gap * sl_spec,
                 lo          = gap * (sl_spec - qnorm(0.975) * se_spec),
                 hi          = gap * (sl_spec + qnorm(0.975) * se_spec)),
  Members    = c(in_memb + sl_memb * B_lo,
                 in_memb + sl_memb * B_hi,
                 gap * sl_memb,
                 gap * (sl_memb - qnorm(0.975) * se_memb),
                 gap * (sl_memb + qnorm(0.975) * se_memb)),
  Difference = c(NA, NA,
                 gap * sl_diff,
                 gap * (sl_diff - qnorm(0.975) * se_diff),
                 gap * (sl_diff + qnorm(0.975) * se_diff))
)
round(out, 2)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------




#---------------------------------------
# Figure v2, UNIFORM assignment rule: distribution of motives
# (no tremble vs with tremble). Two panels: the four model-based types, and
# the two aggregate shares.
#
# This is the robustness counterpart to FigureStructural_R_v2.R. The only
# substantive difference is the source folder: here the within-type
# probabilities p(>= | k, j) are FIXED at 1/|P_k(j)| rather than estimated as
# nuisance parameters. Everything downstream -- which parameters are read, how
# the intervals are built, and every plotting choice -- is identical, so the
# two figures are directly comparable panel by panel.
#
# The `se` column in the results files is a bootstrap STANDARD ERROR, not a CI:
# apply(boot_results, 2, sd, na.rm = TRUE) over 1000 bootstrap replications in
# which participants are resampled with replacement within each role.
# Intervals below are therefore the usual symmetric normal-approximation ones,
#
#   lower = estimate - qnorm(0.975) * se
#   upper = estimate + qnorm(0.975) * se
#
# The aggregate model-based share M = pi_Self + pi_Mob + pi_Mer + pi_MobMer
# equals 1 - pi_Rand exactly, so its standard error is the one estimated for
# pi_Rand. No extra computation is needed.
#---------------------------------------
library(ggplot2)
library(dplyr)
library(readr)

FOLD <- "Clean version - uniform - 5 types"

nt <- read_csv(file.path(FOLD, "results_no_tremble_se_bestmodel.csv"),
               show_col_types = FALSE) %>%
  mutate(estimate = as.numeric(estimate), se = as.numeric(se), spec = "No tremble")

tr <- read_csv(file.path(FOLD, "results_tremble_se_bestmodel.csv"),
               show_col_types = FALSE) %>%
  mutate(estimate = as.numeric(estimate), se = as.numeric(se), spec = "With tremble")

#--- the five estimated type shares ------------------------------------------
types <- bind_rows(nt, tr) %>%
  filter(parameter %in% c("pi_1_self", "pi_2_self_mob", "pi_3_self_mer",
                          "pi_4_self_mob_mer", "pi_5_rand")) %>%
  mutate(type = recode(parameter,
                       pi_1_self         = "Self",
                       pi_2_self_mob     = "Mob",
                       pi_3_self_mer     = "Mer",
                       pi_4_self_mob_mer = "MobMer",
                       pi_5_rand         = "Rand"))

# guard: the uniform folder must contain the same five parameters in both specs
stopifnot(nrow(types) == 10, !any(is.na(types$estimate)), !any(is.na(types$se)))

#--- aggregate model-based share: M = 1 - pi_Rand, so se(M) = se(pi_Rand) -----
agg <- types %>%
  filter(type == "Rand") %>%
  transmute(spec, type = "M", estimate = 1 - estimate, se = se)

df_fig <- bind_rows(types %>% select(spec, type, estimate, se), agg) %>%
  mutate(
    lo   = pmax(0, estimate - qnorm(0.975) * se),
    hi   = pmin(1, estimate + qnorm(0.975) * se),
    type  = factor(type, levels = c("Self", "Mob", "Mer", "MobMer", "M", "Rand")),
    spec  = factor(spec, levels = c("No tremble", "With tremble")),
    panel = factor(ifelse(type %in% c("M", "Rand"), "Aggregate shares",
                                                    "Model-based types"),
                   levels = c("Model-based types", "Aggregate shares"))
  ) %>%
  arrange(type, spec)

df_fig    # inspect before plotting

#--- plot --------------------------------------------------------------------
# Dimensions and all font sizes match FigureMainResults_R and
# FigureStructural_R_v2 so the three figures are typographically consistent.
# NOTE ON PAGE SCALING: if this figure is included with
# adjustbox{width=\linewidth}, it is scaled to exactly \linewidth (~6.5in),
# i.e. 0.765, whereas FigureStructural_R_v2 inside
# adjustwidth{-0.5cm}{-0.5cm} is placed at ~6.9in, i.e. 0.81. To render the two
# at the same size, use the same wrapper for both (see the .tex snippet).
FIG_W <- 8.5    # inches  (same as FigureMainResults_R)
FIG_H <- 4      # inches  (same as FigureMainResults_R)

pt <- function(x) x / .pt   # ggplot2 text `size` is in mm; .pt converts from pt

cols <- c("No tremble"   = "#0072B2",   # blue
          "With tremble" = "#D55E00")   # vermillion

# Two-line tick labels, as plotmath so the Greek renders on the base pdf()
# device (which is WinAnsi-encoded and cannot show literal Greek characters).
xlabs <- list(
  Self   = expression(atop("Selfish",     paste("(", delta, " = 0, ", rho, " = 0)"))),
  Mob    = expression(atop("Mobility",    paste("(", delta, " > 0, ", rho, " = 0)"))),
  Mer    = expression(atop("Meritocracy", paste("(", delta, " = 0, ", rho, " > 0)"))),
  MobMer = expression(atop("Mob & Mer",   paste("(", delta, " > 0, ", rho, " > 0)"))),
  M      = expression(atop("All model-based", group("(", italic(k) %in% italic(M), ")"))),
  Rand   = expression(atop("Random",      "(indifferent)"))
)

# Supplied as a FUNCTION: with free x scales each panel asks only for the
# levels it contains, so positional matching would misalign the labels.
xlab_fun <- function(breaks) unlist(xlabs[as.character(breaks)], use.names = FALSE)

dodge <- position_dodge(width = 0.75)

ggplot(df_fig, aes(x = type, y = estimate, colour = spec, group = spec)) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 1,
                 position = dodge, alpha = 0.85) +
  geom_point(size = 2.6, position = dodge) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * estimate)),
            position = dodge, hjust = -0.20, size = pt(9),
            fontface = "bold", family = "serif", show.legend = FALSE) +
  scale_colour_manual(values = cols, name = "Specification") +
  scale_y_continuous(limits = c(-0.05, 1.03),   # headroom below 0; ticks unchanged
                     breaks = seq(0, 1, 0.1),
                     labels = function(v) paste0(round(100 * v), "%"),
                     expand = expansion(mult = c(0.01, 0.02))) +
  facet_grid(~ panel, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = xlab_fun, expand = expansion(add = 0.8)) +
  labs(x = "Type",
       y = expression(paste("Estimated type probability  ", pi[k]))) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    legend.margin      = margin(t = -2),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    strip.background   = element_rect(fill = "grey94", colour = NA),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9),    # atop() renders smaller,
                                                    # raise to 10-11 if needed
    axis.title.y       = element_text(size = 11, margin = margin(r = 6)),
    axis.title.x       = element_text(size = 11, margin = margin(t = 6)),
    plot.margin        = margin(4, 6, 2, 2)
  )

ggsave("FigureStructuralUniform_R_v2.pdf", width = FIG_W, height = FIG_H)  # base pdf(): no cairo on Windows
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Beliefs and self-serving meritocracy -- UNIFORM ASSIGNMENT RULE
# Robustness counterpart to the beliefs regression table above.
# Estimates equations (reg1) and (reg2) and prints the LaTeX table.
#
# The only substantive difference from the main specification is the source
# folder: here the within-type probabilities p(>= | k, j) are FIXED at
# 1/|P_k(j)| rather than estimated as nuisance parameters. Everything else --
# the dependent variable, the regressors, the estimator, the standard errors --
# is identical, so the two tables are directly comparable row by row.
#
# Dependent variable: the individual posterior probability that participant i
# has a meritocracy motive, i.e. the sum of the posteriors for the Meritocracy
# and Mobility & Meritocracy types, from the model WITH tremble, in pp:
#
#   mer = 100 * (type_3 + type_4)
#
# Scaling by 100 leaves R-squared and all t-statistics unchanged; every
# coefficient and standard error is simply 100 times its value on the 0-1
# scale, which makes the table readable without exponents.
#
# Stata's `robust` for OLS is HC1, so vcovHC(type = "HC1") reproduces it.
#---------------------------------------
library(dplyr)
library(readr)
library(sandwich)
library(lmtest)

FOLD <- "Clean version - uniform - 5 types"

#--- posterior type probabilities (model with tremble) -----------------------
post <- bind_rows(
  read_csv(file.path(FOLD, "posterior_spec_tr_bestmodel.csv"), show_col_types = FALSE) %>%
    rename(participant_code = id_spec)   %>% mutate(role_file = "Spectators"),
  read_csv(file.path(FOLD, "posterior_low_tr_bestmodel.csv"),  show_col_types = FALSE) %>%
    rename(participant_code = id_low)    %>% mutate(role_file = "Low"),
  read_csv(file.path(FOLD, "posterior_med_tr_bestmodel.csv"),  show_col_types = FALSE) %>%
    rename(participant_code = id_medium) %>% mutate(role_file = "Medium"),
  read_csv(file.path(FOLD, "posterior_high_tr_bestmodel.csv"), show_col_types = FALSE) %>%
    rename(participant_code = id_high)   %>% mutate(role_file = "High")
) %>%
  mutate(mer = 100 * (type_3 + type_4))

#--- LABEL CHECK: does each posterior file hold the role it is named for? ----
# This must be a diagonal 4x4 table with 48 on each diagonal cell. Anything
# off-diagonal means the posterior files and df_long disagree about who is who.
lab_chk <- df_long %>%
  distinct(participant_code, parent_status) %>%
  inner_join(post %>% select(participant_code, role_file), by = "participant_code") %>%
  mutate(role_file = factor(role_file, levels = c("Spectators", "High", "Medium", "Low")))

table(lab_chk$role_file, lab_chk$parent_status)

stopifnot(
  nrow(lab_chk) == 192,
  all(as.character(lab_chk$role_file) == as.character(lab_chk$parent_status))
)

#--- merge with beliefs and role ---------------------------------------------
df_reg <- df_long %>%
  distinct(participant_code, parent_status, belief_top, belief_middle, belief_bottom) %>%
  inner_join(post, by = "participant_code") %>%
  mutate(
    Beliefs = belief_top - belief_bottom,
    Members = as.integer(parent_status != "Spectators"),
    Low     = as.integer(parent_status == "Low"),
    Med     = as.integer(parent_status == "Medium"),
    High    = as.integer(parent_status == "High")
  )

stopifnot(
  nrow(df_reg) == 192,
  !any(duplicated(df_reg$participant_code)),
  all(df_reg$Low + df_reg$Med + df_reg$High == df_reg$Members)
)

# means by role, to compare against the main specification
df_reg %>%
  group_by(parent_status) %>%
  summarise(n = n(), mean_mer = mean(mer), mean_Beliefs = mean(Beliefs),
            .groups = "drop")

#--- estimation ---------------------------------------------------------------
m1 <- lm(mer ~ Beliefs + Members + Beliefs:Members, data = df_reg)
m2 <- lm(mer ~ Beliefs + Low + Med + High +
               Beliefs:Low + Beliefs:Med + Beliefs:High, data = df_reg)

ct1 <- coeftest(m1, vcov = vcovHC(m1, type = "HC1"))
ct2 <- coeftest(m2, vcov = vcovHC(m2, type = "HC1"))

ct1
ct2

#--- LaTeX table --------------------------------------------------------------
# three significant digits, as in the main table
fmt <- function(x) formatC(x, digits = 3, format = "fg", flag = "#")

tex_cell <- function(ct, term) {
  if (is.null(term) || !(term %in% rownames(ct))) return(c("", ""))
  e <- ct[term, 1]; s <- ct[term, 2]; p <- ct[term, 4]
  st <- if (p < 0.01) "\\sym{***}" else if (p < 0.05) "\\sym{**}" else
        if (p < 0.10) "\\sym{*}"   else ""
  c(paste0("$", fmt(e), st, "$"), paste0("$(", fmt(s), ")$"))
}

rows <- list(
  list(lab = "Beliefs [$b(r_H) - b(r_L)$]",                            t1 = "Beliefs",         t2 = "Beliefs"),
  list(lab = "Members:",                                               t1 = "Members",         t2 = NULL),
  list(lab = "\\hspace{1em} - Low parent's status",                     t1 = NULL,              t2 = "Low"),
  list(lab = "\\hspace{1em} - Medium parent's status",                  t1 = NULL,              t2 = "Med"),
  list(lab = "\\hspace{1em} - High parent's status",                    t1 = NULL,              t2 = "High"),
  list(lab = "Beliefs $\\times$ Members:",                              t1 = "Beliefs:Members", t2 = NULL),
  list(lab = "\\hspace{1em} - Beliefs $\\times$ Low parent's status",   t1 = NULL,              t2 = "Beliefs:Low"),
  list(lab = "\\hspace{1em} - Beliefs $\\times$ Medium parent's status",t1 = NULL,              t2 = "Beliefs:Med"),
  list(lab = "\\hspace{1em} - Beliefs $\\times$ High parent's status",  t1 = NULL,              t2 = "Beliefs:High"),
  list(lab = "Constant",                                               t1 = "(Intercept)",     t2 = "(Intercept)")
)

body <- character(0)
for (r in rows) {
  col1 <- tex_cell(ct1, r$t1)
  col2 <- tex_cell(ct2, r$t2)
  body <- c(body,
            sprintf("%s & %s & %s \\\\", r$lab, col1[1], col2[1]),
            sprintf(" & %s & %s \\\\",   col1[2], col2[2]))
}

cat(
"\\begin{table}[ht]\n",
"\\centering\n",
"\\begin{footnotesize}\n",
"\\setlength{\\tabcolsep}{7pt}\n",
"\\begin{tabular}{p{0.42\\linewidth}cc}\n",
"\\toprule\n",
" & \\multicolumn{2}{c}{Posterior that $\\rho_i > 0$ ($\\times 100$)} \\\\\n",
"\\cmidrule(lr){2-3}\n",
" & (1) & (2) \\\\\n",
" & Aggregate & Per status \\\\\n",
"\\midrule\n",
paste(body, collapse = "\n"), "\n",
"\\midrule\n",
sprintf("Observations & %d & %d \\\\\n", nobs(m1), nobs(m2)),
sprintf("R-squared & %.3f & %.3f \\\\\n",
        summary(m1)$r.squared, summary(m2)$r.squared),
"\\bottomrule\n",
"\\end{tabular}\n",
"\\end{footnotesize}\n",
"\n\\vspace{1em}\n",
"\\caption{Beliefs and Meritocracy Motive (uniform assignment rule)}\n",
"\\label{tab:beliefs_members_status_uniform}\n",
"\\vspace{1em}\n\n",
"\\begin{minipage}{0.95\\textwidth}\n",
"\\footnotesize\n",
"\\textit{Notes:} The table replicates Table~\\ref{tab:beliefs_members_status}\n",
"using the posteriors of the uniform assignment rule described in this\n",
"appendix, in which the within-type probabilities $p(\\succeq\\mid k,j)$ are\n",
"fixed at $1/|\\mathcal{P}_k(j)|$ rather than estimated. The dependent variable\n",
"is the individual posterior probability that participant $i$ has a meritocracy\n",
"motive, $\\rho_i > 0$, i.e.\\ the sum of the posteriors for the\n",
"\\textit{Meritocracy} and \\textit{Mobility \\& Meritocracy} types. Posteriors\n",
"are computed from the structural model with role-specific trembles reported in\n",
"Figure~\\ref{fig: structural_robust} and are expressed in percentage points.\n",
"All 192 participants enter both regressions; no observation is excluded on the\n",
"basis of its posterior value. $\\text{Beliefs}$ is the difference between\n",
"optimistic and pessimistic beliefs, $b(r_H) - b(r_L)$, in percentage points;\n",
"$\\text{Members} = 1\\{\\text{participant} \\neq \\text{spectator}\\}$;\n",
"$\\text{Low} = 1\\{\\text{parent status = Low}\\}$, and likewise for the other\n",
"statuses. Spectators are the omitted category in both specifications.\n",
"Estimation is by OLS with heteroskedasticity-robust (HC1) standard errors,\n",
"reported in parentheses. \\sym{***} $p<0.01$, \\sym{**} $p<0.05$,\n",
"\\sym{*} $p<0.1$.\n",
"\\end{minipage}\n",
"\\end{table}\n",
sep = "")

#--- Wald test 1: are the three member slopes equal? -------------------------
# In (2) the interactions are relative to spectators, so equality of the three
# interactions is equality of the three member slopes.
# H0: Beliefs:Low = Beliefs:Med = Beliefs:High   (2 restrictions)

# restricted model: status-specific intercepts, one common member slope
m3  <- lm(mer ~ Beliefs + Low + Med + High + Beliefs:Members, data = df_reg)
ct3 <- coeftest(m3, vcov = vcovHC(m3, type = "HC1"))
ct3

# NB: waldtest(m3, m2) fails with "nesting cannot be determined" -- the models
# ARE nested (Members = Low + Med + High) but waldtest matches models by term
# LABEL and "Beliefs:Members" is not a term of m2. Impose it on m2 directly.
car::linearHypothesis(m2,
  c("Beliefs:Low = Beliefs:Med", "Beliefs:Med = Beliefs:High"),
  vcov. = vcovHC(m2, type = "HC1"))

#--- Wald test 2: the belief slope for members, beta_1 + beta_3 --------------
# beta_3 tests whether the belief slope DIFFERS between members and spectators.
# The member slope itself is beta_1 + beta_3, since spectators are the omitted
# category. The standard error needs the covariance term: beta_1_hat and
# beta_3_hat are strongly negatively correlated, so sqrt(V11 + V33) alone would
# badly overstate it.
b1 <- coef(m1)
v1 <- vcovHC(m1, type = "HC1")

tot    <- b1["Beliefs"] + b1["Beliefs:Members"]
se_tot <- sqrt(v1["Beliefs", "Beliefs"] + v1["Beliefs:Members", "Beliefs:Members"] +
               2 * v1["Beliefs", "Beliefs:Members"])

c(slope_members = unname(tot), se = unname(se_tot),
  z = unname(tot / se_tot), p = unname(2 * pnorm(-abs(tot / se_tot))),
  lo = unname(tot - qnorm(0.975) * se_tot),
  hi = unname(tot + qnorm(0.975) * se_tot))

car::linearHypothesis(m1, "Beliefs + Beliefs:Members = 0",
                      vcov. = vcovHC(m1, type = "HC1"))

# reparameterisation: members as the omitted category, so the member slope
# appears directly as a coefficient with its own standard error
df_reg$Spectator <- as.integer(df_reg$parent_status == "Spectators")
m1b  <- lm(mer ~ Beliefs + Spectator + Beliefs:Spectator, data = df_reg)
ct1b <- coeftest(m1b, vcov = vcovHC(m1b, type = "HC1"))
ct1b

# the "Beliefs" row of m1b must reproduce the delta-method figures above
stopifnot(isTRUE(all.equal(unname(ct1b["Beliefs", 1]), unname(tot))),
          isTRUE(all.equal(unname(ct1b["Beliefs", 2]), unname(se_tot))))
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Figure: joint distribution of spectators' choices across meritocracy
# treatments (Appendix: Spectators' choices in No Mobility vs Full Mobility).
#
# Each spectator answers the A1 vs A7 question three times, once per
# meritocracy treatment, so her behaviour is a vector in {0,1}^3. Following the
# appendix text:
#     0 = Full Mobility (A7)      1 = No Mobility (A1)
#     the order of the vector is  (No, Pos, Neg)
#
# For a spectator u is constant across social structures, so
#     V(Full) - V(No) = delta*Dv + rho*Dw_t,
# with Dv > 0 in every treatment and Dw_t = 0 / > 0 / < 0 under No / Pos / Neg.
# A spectator with at least one non-selfish motive (delta > 0 or rho > 0)
# therefore chooses one of three vectors:
#     (0,0,0)  mobility motive only, or both motives
#     (0,0,1)  meritocracy motive only, or both motives
#     (1,0,1)  meritocracy motive only, breaking the tie in No where w is flat
# The other five vectors can arise only from delta = 0 and rho = 0, i.e. from a
# spectator who is indifferent everywhere and chooses at random.
# The dashed line separates the two groups.
#
# Style follows FigureMainResults_R: same dimensions, fonts, palette and theme.
#---------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

update_geom_defaults("text", list(family = "serif"))

#--- one row per spectator, with her 3-choice vector --------------------------
spec_vec <- df_long %>%
  filter(parent_status == "Spectators",
         alternatives  == "A1_A7",
         !is.na(choice)) %>%
  mutate(d = as.integer(choice == "A1")) %>%          # 1 = No Mobility
  select(participant_code, meritocracy, d) %>%
  pivot_wider(names_from = meritocracy, values_from = d) %>%
  drop_na(no, positive, negative) %>%
  mutate(vec = sprintf("(%d,%d,%d)", no, positive, negative))

stopifnot(nrow(spec_vec) == 48, !any(duplicated(spec_vec$participant_code)))

#--- all eight cells, so empty ones still appear as zero-height bars ---------
ADMITTED <- c("(0,0,0)", "(0,0,1)", "(1,0,1)")

df_fig <- expand_grid(no = 0:1, positive = 0:1, negative = 0:1) %>%
  mutate(vec = sprintf("(%d,%d,%d)", no, positive, negative)) %>%
  left_join(count(spec_vec, vec), by = "vec") %>%
  mutate(
    n      = coalesce(n, 0L),
    share  = n / sum(n),
    group  = factor(ifelse(vec %in% ADMITTED, "motive", "no_motive"),
                    levels = c("motive", "no_motive"))
  ) %>%
  # motive block first, each block ordered by decreasing share
  arrange(group, desc(share), vec) %>%
  mutate(vec = factor(vec, levels = vec))

df_fig                                   # inspect before plotting
sum(df_fig$share[df_fig$group == "motive"])   # the 71% in the text

#--- plot --------------------------------------------------------------------
FIG_W <- 8.5    # inches  (same as FigureMainResults_R)
FIG_H <- 4      # inches  (same as FigureMainResults_R)

pt <- function(x) x / .pt   # ggplot2 text `size` is in mm; .pt converts from pt

cols <- c(motive = "#0072B2", no_motive = "grey70")   # blue / grey

# legend labels as plotmath, so the Greek renders on the base pdf() device
leg <- c(motive    = expression(delta > 0 ~ "or" ~ rho > 0),
         no_motive = expression(delta == 0 ~ "and" ~ rho == 0))


ggplot(df_fig, aes(x = vec, y = share, fill = group)) +
  geom_col(width = 0.7) +
  # the line sits between the two groups
  geom_vline(xintercept = length(ADMITTED) + 0.5,
             linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * share)),
            vjust = -0.6, size = pt(9), fontface = "bold",
            show.legend = FALSE) +
  geom_text(aes(y = 0.008, label = sprintf("(n = %d)", n)),
            colour = "grey25", size = pt(8), vjust = 0, show.legend = FALSE) +
  scale_fill_manual(values = cols, labels = leg, name = NULL) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = seq(0, 1, 0.05),
                     labels = function(v) paste0(round(100 * v), "%"),
                     expand = expansion(mult = c(0, 0.10))) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  labs(x = "Vector of choices (No, Pos, Neg)",
       y = "Share of spectators") +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    legend.margin      = margin(t = -2),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9),
    axis.title.y       = element_text(size = 11, margin = margin(r = 6)),
    axis.title.x       = element_text(size = 11, margin = margin(t = 6)),
    plot.margin        = margin(4, 6, 2, 2)
  )

ggsave("FigureJointSpec_R.pdf", width = FIG_W, height = FIG_H)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Selfish motives and asymmetric mobility (Appendix: app: asym)
#
# f^2 (A2) allows mobility between the MIDDLE and LOWER classes; f^3 (A3)
# between the UPPER and MIDDLE classes. Neither v nor w discriminates between
# them: by Propositions (ordermobility) and (ordermeritocracy) they occupy the
# same rank in both orderings, in every meritocracy treatment. The comparison
# is therefore ranked by the SELFISH motive alone, and Proposition
# (orderselfish) gives opposite predictions:
#
#     high- and low-status participants prefer f^2
#     medium-status participants        prefer f^3
#
# Because w is silent here, the meritocracy treatment should make no
# difference. The figure keeps the treatment split visible so that this can be
# read off directly, using the same colour scale as FigureMainResults_R.
#
# Style follows FigureMainResults_R: same dimensions, fonts, palette, theme.
#---------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

update_geom_defaults("text", list(family = "serif"))

#--- shares choosing f^3 over f^2 --------------------------------------------
a23 <- df_long %>%
  filter(alternatives == "A2_A3", !is.na(choice)) %>%
  mutate(a3 = as.integer(choice == "A3"))

stopifnot(nrow(a23) == 192 * 3)

#--- panel 1: members, by parental status ------------------------------------
q1 <- a23 %>%
  filter(parent_status %in% c("Low", "Medium", "High")) %>%
  group_by(parent_status, meritocracy) %>%
  summarise(x = sum(a3), n = n(), .groups = "drop") %>%
  mutate(panel = "Members, by parental status",
         xkey  = paste0("1_", parent_status),
         xlab  = as.character(parent_status),
         merit = as.character(meritocracy))

#--- panel 2: spectators ------------------------------------------------------
q2 <- a23 %>%
  filter(parent_status == "Spectators") %>%
  group_by(meritocracy) %>%
  summarise(x = sum(a3), n = n(), .groups = "drop") %>%
  mutate(panel = "Spectators",
         xkey  = "2_Spectators",
         xlab  = "Spectator",
         merit = as.character(meritocracy),
         parent_status = "Spectators")

#--- combine, with Wald intervals as in FigureMainResults_R ------------------
df_fig <- bind_rows(q1, q2) %>%
  mutate(
    share = x / n,
    se    = sqrt(share * (1 - share) / n),
    lo    = pmax(0, share - qnorm(0.975) * se),
    hi    = pmin(1, share + qnorm(0.975) * se),
    panel = factor(panel, levels = c("Members, by parental status", "Spectators")),
    xkey  = factor(xkey, levels = c("1_Low", "1_Medium", "1_High", "2_Spectators")),
    merit = factor(merit, levels = c("negative", "no", "positive"),
                   labels = c("Negative", "No", "Positive"))
  )

df_fig    # inspect before plotting

#--- plot --------------------------------------------------------------------
FIG_W <- 8.5    # inches  (same as FigureMainResults_R)
FIG_H <- 4      # inches  (same as FigureMainResults_R)

pt <- function(x) x / .pt   # ggplot2 text `size` is in mm; .pt converts from pt

cols <- c("Negative" = "#D55E00",   # vermillion
          "No"       = "#0072B2",   # blue
          "Positive" = "#009E73")   # green

xlabs <- setNames(df_fig$xlab, as.character(df_fig$xkey))

dodge <- position_dodge(width = 0.9)

ggplot(df_fig, aes(x = xkey, y = share, colour = merit)) +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             colour = "grey45", linewidth = 0.5) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 1,
                 position = dodge, alpha = 0.85) +
  geom_point(size = 2.6, position = dodge) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * share)),
            position = dodge, hjust = -0.20, size = pt(8),
            fontface = "bold", show.legend = FALSE) +
  geom_text(aes(y = 0.02, label = sprintf("(n = %d)", n)),
            position = dodge, colour = "grey25", size = pt(8),
            vjust = 0, show.legend = FALSE) +
  facet_grid(~ panel, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = cols, name = "Meritocracy treatment") +
  scale_y_continuous(limits = c(0, 1.03),
                     breaks = seq(0, 1, 0.1),
                     labels = function(v) paste0(round(100 * v), "%"),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_x_discrete(labels = xlabs, expand = expansion(add = 0.8)) +
  labs(x = "Parental status",
       y = expression(paste("Share choosing  ", italic(f)^3, "  over  ", italic(f)^2))) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    legend.margin      = margin(t = -2),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    strip.background   = element_rect(fill = "grey94", colour = NA),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9, lineheight = 0.95),
    axis.title.y       = element_text(size = 11, margin = margin(r = 6)),
    axis.title.x       = element_text(size = 11, margin = margin(t = 6)),
    plot.margin        = margin(4, 6, 2, 2)
  )

ggsave("FigureA2A3_R.pdf", width = FIG_W, height = FIG_H)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# TESTS for Selfish motives and asymmetric mobility (Appendix: app: asym)
#
# The model makes a POINT prediction for each group, not a between-group
# comparison. Proposition (orderselfish) gives f^2 > f^3 for high- and
# low-status participants and f^3 > f^2 for medium-status ones, while v and w
# are silent for this pair. For a spectator u is constant too, so NOTHING ranks
# f^2 against f^3 and the model predicts indifference, i.e. probability 1/2.
#
# Every test therefore compares a share against 1/2, one-sided in the predicted
# direction for members and two-sided for spectators.
#
# Tests are run SEPARATELY BY TREATMENT. Within a treatment each participant
# contributes exactly one observation, so the exact binomial is valid whatever
# the heterogeneity across people. Pooling the 3 x 48 choices would be valid
# only under the null of pure indifference; the members hold stable
# preferences (see the consistency block below), so pooling would overstate
# precision for them. Running by treatment also replicates each test three
# times, which is itself a check: v and w are silent here, so the model
# predicts the meritocracy treatment makes no difference.
#---------------------------------------
PRED <- c(Low = "f2", Medium = "f3", High = "f2", Spectators = NA)

#--- shares choosing f^3 over f^2 --------------------------------------------
a23 <- df_long %>%
  filter(alternatives == "A2_A3", !is.na(choice)) %>%
  mutate(a3 = as.integer(choice == "A3"))

stopifnot(nrow(a23) == 192 * 3)

#--- main table: share choosing the predicted option, by group and treatment --
tests <- a23 %>%
  group_by(parent_status, meritocracy) %>%
  summarise(n = n(), n_f3 = sum(a3), .groups = "drop") %>%
  mutate(
    predicted = PRED[as.character(parent_status)],
    # x counts choices of the PREDICTED option; for spectators we report f3
    x     = ifelse(is.na(predicted), n_f3, ifelse(predicted == "f3", n_f3, n - n_f3)),
    share = x / n,
    alt   = ifelse(is.na(predicted), "two.sided", "greater")
  ) %>%
  # binom.test takes single numbers, so step through the rows one at a time.
  # The CI is always the two-sided one, even where the test is one-sided.
  rowwise() %>%
  mutate(p     = binom.test(x, n, p = 0.5, alternative = alt)$p.value,
         ci_lo = binom.test(x, n, p = 0.5)$conf.int[1],
         ci_hi = binom.test(x, n, p = 0.5)$conf.int[2]) %>%
  ungroup() %>%
  select(parent_status, meritocracy, predicted, x, n, share, p, ci_lo, ci_hi) %>%
  arrange(factor(parent_status, levels = c("Low", "Medium", "High", "Spectators")),
          factor(meritocracy,   levels = c("no", "positive", "negative")))

print(tests, n = Inf)

# the two summary statements used in the text
cat("\nlargest p across the nine member tests :",
    format(max(tests$p[tests$parent_status != "Spectators"]), digits = 3), "\n")
cat("spectator p, each treatment           :",
    paste(format(tests$p[tests$parent_status == "Spectators"], digits = 3),
          collapse = ", "), "\n")

#--- Low vs High: the model predicts no difference (both prefer f2) ----------
lh <- bind_rows(lapply(c("no", "positive", "negative"), function(t) {
  tb <- a23 %>%
    filter(meritocracy == t, parent_status %in% c("Low", "High")) %>%
    mutate(grp = factor(parent_status, levels = c("Low", "High")),
           ch  = factor(ifelse(a3 == 1, "f3", "f2"), levels = c("f3", "f2"))) %>%
    with(table(grp, ch))
  ft <- fisher.test(tb)                                   # two-sided
  tibble(meritocracy = t,
         share_f3_Low  = tb["Low",  "f3"] / sum(tb["Low",  ]),
         share_f3_High = tb["High", "f3"] / sum(tb["High", ]),
         odds_ratio = unname(ft$estimate), p = ft$p.value)
}))
print(lh)

#--- spectators: consistency across the three treatments ---------------------
# Under the model spectators are indifferent, so each of their three choices is
# an independent coin flip and the number choosing identically in all three is
# Binomial(48, p^3 + (1-p)^3) with p the observed rate. Members are shown for
# contrast: there u DOES rank the pair, so consistency should be much higher.
consist <- a23 %>%
  select(participant_code, parent_status, meritocracy, a3) %>%
  pivot_wider(names_from = meritocracy, values_from = a3) %>%
  mutate(same = (no == positive) & (positive == negative)) %>%
  group_by(parent_status) %>%
  summarise(consistent = sum(same), n = n(), .groups = "drop") %>%
  left_join(a23 %>% group_by(parent_status) %>% summarise(p_f3 = mean(a3), .groups = "drop"),
            by = "parent_status") %>%
  mutate(expected_if_random = n * (p_f3^3 + (1 - p_f3)^3),
         p = NA_real_) %>%
  rowwise() %>%
  mutate(p = binom.test(consistent, n, p = p_f3^3 + (1 - p_f3)^3)$p.value) %>%
  ungroup()
print(consist)
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Beliefs and preference types: specification (1) run separately on EACH type
#
# Table (beliefs_members_status) regresses the posterior probability of a
# MERITOCRACY motive -- the sum of the Meritocracy and Mobility & Meritocracy
# posteriors -- on beliefs, membership and their interaction. That sum hides
# where the belief effect actually lands, and in particular cannot say whether
# optimism raises the meritocracy motive specifically or shifts members out of
# the Selfish classification generally. This script estimates the same
# specification (1),
#
#   100 * posterior_k(i) = b0 + b1*Beliefs + b2*Members + b3*Beliefs:Members
#
# once for each of the five types, so the decomposition is visible.
#
# NOTE ON INTERPRETATION. The five posteriors sum to one for every participant,
# so 100 * sum_k posterior_k = 100 is a constant: every coefficient must sum to
# ZERO across the five columns (and the constants to 100). The columns are
# mechanically linked, are not independent tests, and no multiplicity
# correction is meaningful. The adding-up is also a useful check on the merge.
#
# Posteriors come from the model WITH tremble. Stata's `robust` for OLS is HC1,
# so vcovHC(type = "HC1") reproduces it.
#---------------------------------------
library(dplyr)
library(readr)
library(sandwich)
library(lmtest)

FOLD <- "Clean version - no overlap - 5 types"

#--- posterior type probabilities (model with tremble) -----------------------
post <- bind_rows(
  read_csv(file.path(FOLD, "posterior_spec_tr_bestmodel.csv"), show_col_types = FALSE) %>%
    rename(participant_code = id_spec)   %>% mutate(role_file = "Spectators"),
  read_csv(file.path(FOLD, "posterior_low_tr_bestmodel.csv"),  show_col_types = FALSE) %>%
    rename(participant_code = id_low)    %>% mutate(role_file = "Low"),
  read_csv(file.path(FOLD, "posterior_med_tr_bestmodel.csv"),  show_col_types = FALSE) %>%
    rename(participant_code = id_medium) %>% mutate(role_file = "Medium"),
  read_csv(file.path(FOLD, "posterior_high_tr_bestmodel.csv"), show_col_types = FALSE) %>%
    rename(participant_code = id_high)   %>% mutate(role_file = "High")
)

#--- LABEL CHECK: does each posterior file hold the role it is named for? ----
# Must be a diagonal 4x4 table with 48 in each diagonal cell. Anything
# off-diagonal means the posterior files and df_long disagree about who is who.
lab_chk <- df_long %>%
  distinct(participant_code, parent_status) %>%
  inner_join(post %>% select(participant_code, role_file), by = "participant_code") %>%
  mutate(role_file = factor(role_file, levels = c("Spectators", "High", "Medium", "Low")))

table(lab_chk$role_file, lab_chk$parent_status)

stopifnot(
  nrow(lab_chk) == 192,
  all(as.character(lab_chk$role_file) == as.character(lab_chk$parent_status))
)

#--- merge with beliefs and role ---------------------------------------------
df_reg <- df_long %>%
  distinct(participant_code, parent_status, belief_top, belief_middle, belief_bottom) %>%
  inner_join(post, by = "participant_code") %>%
  mutate(
    Beliefs = belief_top - belief_bottom,
    Members = as.integer(parent_status != "Spectators")
  )

stopifnot(
  nrow(df_reg) == 192,
  !any(duplicated(df_reg$participant_code)),
  # the five posteriors must sum to one for every participant
  all(abs(rowSums(df_reg[, paste0("type_", 1:5)]) - 1) < 1e-8)
)

#--- specification (1) by type: table with types in COLUMNS --------------------
types <- c(Selfish = "type_1", Mobility = "type_2", Meritocracy = "type_3",
           `Mob & Mer` = "type_4", Random = "type_5")

fits <- lapply(types, function(v) {
  m  <- lm(100 * df_reg[[v]] ~ Beliefs + Members + Beliefs:Members, data = df_reg)
  ct <- coeftest(m, vcov = vcovHC(m, type = "HC1"))
  V  <- vcovHC(m, type = "HC1")
  # belief slope for members is b1 + b3; its standard error needs the
  # covariance term, since b1_hat and b3_hat are strongly negatively correlated
  sl <- ct["Beliefs", 1] + ct["Beliefs:Members", 1]
  se <- sqrt(V["Beliefs", "Beliefs"] + V["Beliefs:Members", "Beliefs:Members"] +
             2 * V["Beliefs", "Beliefs:Members"])
  list(ct = ct, r2 = summary(m)$r.squared, n = nobs(m),
       sl = c(sl, se, 2 * pnorm(-abs(sl / se))))
})

# adding-up check: each coefficient must sum to zero across the five columns
adding_up <- sapply(c("Beliefs", "Members", "Beliefs:Members"),
                    function(term) sum(sapply(fits, function(f) f$ct[term, 1])))
adding_up                                   # all three must be ~ 0
stopifnot(all(abs(adding_up) < 1e-8),
          abs(sum(sapply(fits, function(f) f$ct["(Intercept)", 1])) - 100) < 1e-8)

# compact summary, for reading before the LaTeX is generated
t(sapply(fits, function(f)
  c(b1_spec = f$ct["Beliefs", 1], b3_int = f$ct["Beliefs:Members", 1],
    se_b3 = f$ct["Beliefs:Members", 2], p_b3 = f$ct["Beliefs:Members", 4],
    slope_mb = f$sl[1], se_mb = f$sl[2], p_mb = f$sl[3]))) %>% round(4)

#--- LaTeX table --------------------------------------------------------------
fmt <- function(x) formatC(x, digits = 3, format = "fg", flag = "#")
star <- function(p) if (p < 0.01) "\\sym{***}" else if (p < 0.05) "\\sym{**}" else
                    if (p < 0.10) "\\sym{*}"   else ""

# one row of estimates + one row of standard errors, across the five columns
mk <- function(lab, get) {
  e <- sapply(fits, function(f) { z <- get(f); paste0("$", fmt(z[1]), star(z[3]), "$") })
  s <- sapply(fits, function(f) { z <- get(f); paste0("$(", fmt(z[2]), ")$") })
  c(paste(lab, paste(e, collapse = " & "), sep = " & ") |> paste("\\\\"),
    paste("",  paste(s, collapse = " & "), sep = " & ") |> paste("\\\\"))
}

body <- c(
  mk("Beliefs [$b(r_H) - b(r_L)$]", function(f) f$ct["Beliefs", c(1, 2, 4)]),
  mk("Members",                     function(f) f$ct["Members", c(1, 2, 4)]),
  mk("Beliefs $\\times$ Members",   function(f) f$ct["Beliefs:Members", c(1, 2, 4)]),
  mk("Constant",                    function(f) f$ct["(Intercept)", c(1, 2, 4)]),
  "\\midrule",
  mk("Belief slope for members",    function(f) f$sl)
)

cat(
"\\begin{table}[ht]\n\\centering\n\\begin{footnotesize}\n",
"\\setlength{\\tabcolsep}{5pt}\n",
"\\begin{tabular}{lccccc}\n\\toprule\n",
" & \\multicolumn{5}{c}{Posterior probability of each type ($\\times 100$)} \\\\\n",
"\\cmidrule(lr){2-6}\n",
" & (1) & (2) & (3) & (4) & (5) \\\\\n",
" & Selfish & Mobility & Meritocracy & Mob \\& Mer & Random \\\\\n",
"\\midrule\n",
paste(body, collapse = "\n"), "\n\\midrule\n",
sprintf("Observations & %s \\\\\n", paste(rep(fits[[1]]$n, 5), collapse = " & ")),
sprintf("R-squared & %s \\\\\n",
        paste(sprintf("%.3f", sapply(fits, `[[`, "r2")), collapse = " & ")),
"\\bottomrule\n\\end{tabular}\n\\end{footnotesize}\n",
"\n\\vspace{1em}\n",
"\\caption{Beliefs and preference types}\n",
"\\label{tab:beliefs_by_type}\n",
"\\vspace{1em}\n\n",
"\\begin{minipage}{0.95\\textwidth}\n\\footnotesize\n",
"\\textit{Notes:} Each column reports specification~(1) of\n",
"Table~\\ref{tab:beliefs_members_status} estimated on the individual posterior\n",
"probability of a single type, in percentage points. Posteriors are computed\n",
"from the structural model with role-specific trembles reported in\n",
"Figure~\\ref{fig: structural}. Because the five posteriors sum to one for every\n",
"participant, the coefficients in each row sum to zero across columns (and the\n",
"constants to 100); the five columns are therefore mechanically linked and are\n",
"not independent tests. $\\text{Beliefs}=b(r_H)-b(r_L)$ is expressed in\n",
"percentage points and $\\text{Members}=1\\{\\text{participant}\\neq\\text{spectator}\\}$;\n",
"spectators are the omitted category. The belief slope for members is\n",
"$\\beta_1+\\beta_3$, with standard error obtained by the delta method.\n",
"Estimation is by OLS with heteroskedasticity-robust (HC1) standard errors,\n",
"reported in parentheses. \\sym{***} $p<0.01$, \\sym{**} $p<0.05$,\n",
"\\sym{*} $p<0.1$.\n",
"\\end{minipage}\n\\end{table}\n", sep = "")
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Table: estimated tremble parameters by role (tab: tremble)
#
# First row: role-specific tremble epsilon_j of the structural model with
# trembles (Section 5), with bootstrap standard errors. Second row: implied
# probability of choosing the strictly preferred option, 1 - epsilon_j / 2.
# The transformation is linear, so its standard error is SE(epsilon_j) / 2.
#
# Reads results_tremble_se_bestmodel.csv, produced by
# 02_stima_tremble_se_bestmodel.R in the main estimation folder.
#---------------------------------------
FOLD <- "Clean version - no overlap - 5 types"

res <- read.csv(file.path(FOLD, "results_tremble_se_bestmodel.csv"), stringsAsFactors = FALSE)
est <- setNames(res$estimate, res$parameter)
se  <- setNames(res$se,       res$parameter)

roles   <- c("spec", "high", "medium", "low")
tr_est  <- est[paste0("tremble_", roles)]
tr_se   <- se [paste0("tremble_", roles)]
pref_est <- 1 - tr_est / 2
pref_se  <- tr_se / 2

fmt3 <- function(x) sprintf("%.3f", x)
cell <- function(e, s) paste0("& ", fmt3(e), " \\hspace{0.2em} (", fmt3(s), ")")

cat(
"\\begin{table}[t]\n\\centering\n\\tabsetup\n",
"\\begin{tabular}{c|cccc}\n",
"& \\text{Spectators} & \\text{High} & \\text{Medium} & \\text{Low} \\\\ \\hline\n",
"\\text{Tremble}\n",
paste(cell(tr_est, tr_se), collapse = "\n"), " \\\\\n",
"\\text{Implied prob. of preferred choice}\n",
paste(cell(pref_est, pref_se), collapse = "\n"), " \\\\\n",
"\\end{tabular}\n",
"\n\\vspace{1em}\n",
"\\caption{Estimated tremble parameters by role}\n",
"\\label{tab: tremble}\n",
"\\tabnotes{\\textit{Notes:} The first row reports the estimated role-specific tremble parameter $\\varepsilon_j$. The second row reports the corresponding probability of choosing the option strictly preferred by the underlying preference relation, $1-\\varepsilon_j/2$. Columns correspond to spectators and participants with high, medium, and low parental status. Standard errors are reported in parentheses; standard errors in the second row are obtained from the corresponding transformation of the tremble parameter.}\n",
"\\end{table}\n", sep = "")
rm(list = setdiff(ls(), KEEP))
#---------------------------------------





#---------------------------------------
# Table: estimated type shares and restricted-model comparisons (tab:shares)
#
# Column (1) is the unrestricted model with trembles. Columns (2)-(3) drop a
# motive (all types requiring mobility, resp. meritocracy); columns (4)-(6)
# drop a single type (Mobility, Meritocracy, Mobility & Meritocracy).
#
# Type shares and bootstrap SE come from the results_*.csv files written by the
# 02_ scripts. Excluded types have share and SE fixed at exactly 0 and are left
# blank. The results files store the minimised objective, i.e. the NEGATIVE
# log-likelihood.
#
# LRT p-values come from the bootstrap_LRT_*.rds files written by the 03_
# scripts (parametric bootstrap, B = 1,000, estimator (1 + #{LR* >= LR})/(1 + B)).
# With B = 1,000 the smallest attainable value is 1/1001, printed as "<0.001".
# The block checks that each stored observed LR equals twice the log-likelihood
# gap computed from the results files.
#---------------------------------------
dir_main    <- "Clean version - no overlap - 5 types"
dir_motives <- "Clean version - no overlap - 5 types drop motives"
dir_types   <- "Clean version - no overlap - 5 types drop types"

specs <- list(
  list(res = file.path(dir_main,    "results_tremble_se_bestmodel.csv"),       lrt = NULL),
  list(res = file.path(dir_motives, "results_tremble_se_drop_mobmer_mob.csv"), lrt = file.path(dir_motives, "bootstrap_LRT_drop_mobmer_mob.rds")),
  list(res = file.path(dir_motives, "results_tremble_se_drop_mobmer_mer.csv"), lrt = file.path(dir_motives, "bootstrap_LRT_drop_mobmer_mer.rds")),
  list(res = file.path(dir_types,   "results_tremble_se_drop_mob.csv"),        lrt = file.path(dir_types,   "bootstrap_LRT_drop_mob.rds")),
  list(res = file.path(dir_types,   "results_tremble_se_drop_mer.csv"),        lrt = file.path(dir_types,   "bootstrap_LRT_drop_mer.rds")),
  list(res = file.path(dir_types,   "results_tremble_se_drop_mobmer.csv"),     lrt = file.path(dir_types,   "bootstrap_LRT_drop_mobmer.rds"))
)
read_results <- function(path) {
  r <- read.csv(path, stringsAsFactors = FALSE)
  list(est = setNames(r$estimate, r$parameter), se = setNames(r$se, r$parameter))
}
res <- lapply(specs, function(s) read_results(s$res))
lrt <- lapply(specs, function(s) if (is.null(s$lrt)) NULL else readRDS(s$lrt))

#--- consistency checks -------------------------------------------------------
ll_full <- -res[[1]]$est[["Log_lik"]]
for (j in 2:length(specs)) {
  lr_csv <- 2 * (ll_full + res[[j]]$est[["Log_lik"]])
  stopifnot(abs(lrt[[j]]$LR_obs - lr_csv) < 1e-6)
}
n_valid <- vapply(lrt[-1], function(x) x$n_valid, numeric(1))
if (any(n_valid != 1000))
  warning("Some LRTs do not have 1,000 valid replications: ", paste(n_valid, collapse = ", "))

#--- LaTeX table --------------------------------------------------------------
fmt3 <- function(x) sprintf("%.3f", x)
pad  <- function(x) formatC(x, width = -10)
is_dropped <- function(r, p) r$est[[p]] == 0 && r$se[[p]] == 0

pi_par   <- c("pi_1_self", "pi_2_self_mob", "pi_3_self_mer", "pi_4_self_mob_mer", "pi_5_rand")
pi_label <- c("$\\pi_{\\text{Selfish}}$", "$\\pi_{\\text{Mobility}}$", "$\\pi_{\\text{Meritocracy}}$",
              "$\\pi_{\\text{Mob \\& Mer}}$", "$\\pi_{\\text{Random}}$")

pi_rows <- vapply(seq_along(pi_par), function(k) {
  p  <- pi_par[k]
  e  <- vapply(res, function(r) if (is_dropped(r, p)) "" else fmt3(r$est[[p]]), "")
  s  <- vapply(res, function(r) if (is_dropped(r, p)) "" else paste0("(", fmt3(r$se[[p]]), ")"), "")
  sep <- if (k < length(pi_par)) " \\\\[4pt]" else " \\\\"
  paste0(pi_label[k], "\n",
         "    & ", paste(pad(e), collapse = " & "), " \\\\\n",
         "    & ", paste(pad(s), collapse = " & "), sep)
}, "")

ll    <- vapply(res, function(r) sprintf("$%.3f$", -r$est[["Log_lik"]]), "")
pvals <- vapply(lrt, function(x) {
  if (is.null(x)) "--" else if (x$p_value < 0.001) "$<0.001$" else sprintf("$%.3f$", x$p_value)
}, "")

cat(
"\\begin{table}[t]\n\\centering\n\\tabsetup\n",
"\\begin{tabular}{lcccccc}\n\\toprule\n",
" & & \\multicolumn{2}{c}{Motives} & \\multicolumn{3}{c}{Types} \\\\\n",
"\\cmidrule(lr){3-4} \\cmidrule(lr){5-7}\n",
" & (1) & (2) & (3) & (4) & (5) & (6) \\\\\n",
"\\midrule\n",
paste(pi_rows, collapse = "\n"), "\n",
"\\midrule\n",
"Log-lik.\n    & ", paste(ll, collapse = " & "), " \\\\\n",
"LRT ($p$-value)\n    & ", paste(pvals, collapse = " & "), " \\\\\n",
"\\bottomrule\n\\end{tabular}\n",
"\n\\vspace{1em}\n",
"\\caption{Estimated type shares and restricted-model comparisons}\n",
"\\label{tab:shares}\n",
"\\tabnotes{\\textit{Notes:} Specification (1) is the unrestricted model reported in Figure~\\ref{fig: structural}. All specifications allow for role-specific trembles. Under \\textit{Motives}, Specification (2) excludes the mobility motive and Specification (3) excludes the meritocracy motive. Under \\textit{Types}, Specifications (4)--(6) exclude, respectively, the \\textit{Mobility}, \\textit{Meritocracy}, and \\textit{Mobility \\& Meritocracy} types. Standard errors, reported in parentheses, are obtained from 1,000 bootstrap replications, resampling participants within role. The likelihood-ratio tests compare each restricted specification with Specification (1). Because the restrictions lie on the boundary of the parameter space, the LR statistic does not follow a $\\chi^2$ distribution; $p$-values are obtained by parametric bootstrap ($B=1{,}000$ datasets simulated under the restricted MLE), re-estimating both models on each dataset using multi-start optimization and using the $(1+\\#\\{LR^*\\geq LR\\})/(1+B)$ estimator.}\n",
"\\end{table}\n", sep = "")
rm(list = setdiff(ls(), KEEP))
#---------------------------------------
