#---------------------------------------
# Parametric bootstrap LRT
# H0 (restricted): drop Self_Mer only   (pi[3] = 0)
# H1 (full)      : unrestricted 5-motive tremble model
#
# Statistic: LR = 2 * (loglik_full - loglik_restricted)  >= 0
#
# Why parametric bootstrap: the restriction sets one mixture weight to the
# boundary (pi = 0) and the corresponding conditional parameters are
# unidentified under H0, so the asymptotic null is a chi-bar-squared mixture
# (~ 1/2 chi^2_0 + 1/2 chi^2_1) AND n is small (48 per subsample). The
# parametric bootstrap sidesteps both issues: simulate data under the
# restricted MLE, refit both models, read the p-value off the empirical LR.
#---------------------------------------
source("stima_pre_process.R")


#=======================================================================
# 0. Restricted-model machinery (must match the 01/02 drop_mer files)
#=======================================================================
drop_pi <- c(3)      # Self_Mer
keep_pi <- c(1, 2, 4, 5)

# Pin only blocks that feed Self_Mer ALONE. The low 'a' block (15:19) is shared
# with Self_Mob_Mer, which survives here, so it stays free (identified through
# Self_Mob_Mer). Spec has no Self_Mer parameter.
pin_idx <- c(20:23,   # medium Self_Mer (a)
             38:42)   # high   Self_Mer (a)

negloglik_restricted <- function(par, C_spec, C_high, C_medium, C_low) {
  pi <- numeric(5)
  pi[keep_pi] <- softmax(par[keep_pi])
  pi[drop_pi] <- 0

  tremble_spec   <- par[6]
  tremble_high   <- par[7]
  tremble_medium <- par[8]
  tremble_low    <- par[9]

  pmt_spec   <- softmax(par[10:14])
  pmt_low    <- softmax(par[15:19])
  pmt_medium <- c(softmax(par[20:23]), softmax(par[24:25]), softmax(par[26:37]))
  pmt_high   <- c(softmax(par[38:42]), softmax(par[43:47]), softmax(par[48:72]))

  m_comp_spec   <- prob_profile_given_type_spec(pmt_spec,   tremble_lik(tremble_spec,   C_spec))
  m_comp_low    <- prob_profile_given_type_low(pmt_low,     tremble_lik(tremble_low,    C_low))
  m_comp_medium <- prob_profile_given_type_medium(pmt_medium, tremble_lik(tremble_medium, C_medium))
  m_comp_high   <- prob_profile_given_type_high(pmt_high,   tremble_lik(tremble_high,   C_high))

  log_pi <- log(pi)
  log_terms_spec   <- sweep(log(m_comp_spec),   2, log_pi, "+")
  log_terms_low    <- sweep(log(m_comp_low),    2, log_pi, "+")
  log_terms_medium <- sweep(log(m_comp_medium), 2, log_pi, "+")
  log_terms_high   <- sweep(log(m_comp_high),   2, log_pi, "+")

  li_spec   <- apply(log_terms_spec,   1, logsumexp)
  li_low    <- apply(log_terms_low,    1, logsumexp)
  li_medium <- apply(log_terms_medium, 1, logsumexp)
  li_high   <- apply(log_terms_high,   1, logsumexp)

  ll <- sum(li_spec) + sum(li_high) + sum(li_low) + sum(li_medium)
  if (!is.finite(ll)) return(1e100)
  -ll
}

n_par <- 72

# Bounds for the FULL model (only trembles constrained)
lower_full <- rep(-Inf, n_par); upper_full <- rep(Inf, n_par)
lower_full[6:9] <- 0; upper_full[6:9] <- 0.5

# Bounds for the RESTRICTED model (trembles + pinned conditionals)
lower_restr <- lower_full; upper_restr <- upper_full
lower_restr[pin_idx] <- 0; upper_restr[pin_idx] <- 0


#=======================================================================
# 1. Load the two fitted models -> get observed LR and the H0 DGP
#=======================================================================
# Full model fit (unrestricted 5-motive). Produced by 01_stima_tremble_1000starts.
full_env <- new.env()
load("stima_tremble_1000starts.Rdata", envir = full_env)
par_full_hat <- if (!is.null(full_env$best_par)) full_env$best_par else full_env$fit$par
ll_full_obs  <- -full_env$fit$objective

# Restricted model fit. Produced by 01_stima_tremble_drop_mer_1000starts.
restr_env <- new.env()
load("stima_tremble_drop_mer_1000starts.Rdata", envir = restr_env)
par_restr_hat <- if (!is.null(restr_env$best_par)) restr_env$best_par else restr_env$fit$par
ll_restr_obs  <- -restr_env$fit$objective

LR_obs <- 2 * (ll_full_obs - ll_restr_obs)
cat(sprintf("Observed loglik  full = %.4f\n", ll_full_obs))
cat(sprintf("Observed loglik restr = %.4f\n", ll_restr_obs))
cat(sprintf("Observed LR statistic = %.4f\n", LR_obs))
if (LR_obs < -1e-6) warning("LR_obs < 0: the restricted fit beat the full fit -> rerun the full multistart, it missed the optimum.")


#=======================================================================
# 2. Build the H0 data-generating process from the restricted MLE
#=======================================================================
# Recover the H0 parameters in interpretable form.
h0_pi <- numeric(5)
h0_pi[keep_pi] <- softmax(par_restr_hat[keep_pi])
h0_pi[drop_pi] <- 0

h0_tremble <- c(spec = par_restr_hat[6], high = par_restr_hat[7],
                medium = par_restr_hat[8], low = par_restr_hat[9])

# Profile-distribution (per type, per context) as a matrix:
#   rows = profile index, cols = type (1..5). This is exactly the M matrix
#   evaluated at the fitted pmt, i.e. P(profile | type) for that context.
pmt_spec   <- softmax(par_restr_hat[10:14])
pmt_low    <- softmax(par_restr_hat[15:19])
pmt_medium <- c(softmax(par_restr_hat[20:23]), softmax(par_restr_hat[24:25]), softmax(par_restr_hat[26:37]))
pmt_high   <- c(softmax(par_restr_hat[38:42]), softmax(par_restr_hat[43:47]), softmax(par_restr_hat[48:72]))

# Build P(profile | type) matrices by applying each context's M to an
# identity in "count space": prob_profile_given_type_*(pmt, I) returns, for a
# one-hot profile vector, the column of M. We instead reconstruct M directly
# by passing identity rows through the same function used in the likelihood.
profile_dist_matrix <- function(prob_fun, pmt, n_prof) {
  # Feed the function a set of "pseudo-likelihood" rows = identity, so the
  # output row p corresponds to M[p, ] (the per-type weight on profile p).
  I <- diag(n_prof)
  prob_fun(pmt, I)   # rows = profiles, cols = types ; columns may not sum to 1
}

# NOTE: columns of M are P(profile | type) and DO sum to 1 by construction
# (each type's column is either a one-hot, a shared softmax block, or a full
# softmax over its profiles). We normalise defensively to guard rounding.
Pspec   <- profile_dist_matrix(prob_profile_given_type_spec,   pmt_spec,   nrow(profiles_spec))
Plow    <- profile_dist_matrix(prob_profile_given_type_low,    pmt_low,    nrow(profiles_low))
Pmedium <- profile_dist_matrix(prob_profile_given_type_medium, pmt_medium, nrow(profiles_medium))
Phigh   <- profile_dist_matrix(prob_profile_given_type_high,   pmt_high,   nrow(profiles_high))

norm_cols <- function(M) sweep(M, 2, pmax(colSums(M), .Machine$double.eps), "/")
Pspec <- norm_cols(Pspec); Plow <- norm_cols(Plow)
Pmedium <- norm_cols(Pmedium); Phigh <- norm_cols(Phigh)

# Number of participants per subsample (kept fixed across bootstrap reps).
n_spec   <- dim(dati_stima_spec)[1]
n_high   <- dim(dati_stima_high)[1]
n_medium <- dim(dati_stima_medium)[1]
n_low    <- dim(dati_stima_low)[1]

# Simulate ONE subsample's C/S/I array under H0.
#   profiles_mat : the profiles_* matrix (rows = profile, cols = 12 choices, 0/1/NA)
#   Pmat         : P(profile | type) matrix (rows = profile, cols = type)
#   pi           : 5-vector mixture weights
#   tremble      : scalar tremble for this context
#   n            : number of participants to simulate
simulate_subsample <- function(profiles_mat, Pmat, pi, tremble, n) {
  n_prof    <- nrow(profiles_mat)
  n_choices <- ncol(profiles_mat)

  # 1. draw a type for each participant
  types <- sample.int(5, n, replace = TRUE, prob = pi)

  # 2. draw a profile given type, then 3. generate observed 0/1 choices
  X <- matrix(NA_integer_, nrow = n, ncol = n_choices)
  for (i in seq_len(n)) {
    k    <- types[i]
    prof <- sample.int(n_prof, 1, prob = Pmat[, k])
    intended <- profiles_mat[prof, ]   # vector of 0/1/NA over the 12 choices
    for (j in seq_len(n_choices)) {
      a <- intended[j]
      if (is.na(a)) {
        X[i, j] <- rbinom(1, 1, 0.5)                     # indifferent -> 50/50
      } else if (a == 0) {
        X[i, j] <- rbinom(1, 1, tremble / 2)             # intended 0, slip to 1 w.p. t/2
      } else {
        X[i, j] <- rbinom(1, 1, 1 - tremble / 2)         # intended 1, stay 1 w.p. 1 - t/2
      }
    }
  }

  # 4. recompute C/S/I counts exactly as the real pipeline (evaluate_choice)
  cnts <- t(apply(X, 1, evaluate_choice, mat = profiles_mat))
  array(cnts, c(n, n_prof, 3))
}

simulate_all <- function() {
  list(
    spec   = simulate_subsample(profiles_spec,   Pspec,   h0_pi, h0_tremble["spec"],   n_spec),
    high   = simulate_subsample(profiles_high,   Phigh,   h0_pi, h0_tremble["high"],   n_high),
    medium = simulate_subsample(profiles_medium, Pmedium, h0_pi, h0_tremble["medium"], n_medium),
    low    = simulate_subsample(profiles_low,    Plow,    h0_pi, h0_tremble["low"],    n_low)
  )
}


#=======================================================================
# 3. Refit BOTH models on a simulated dataset -> one bootstrap LR
#=======================================================================
# Warm starts: use the observed MLEs (the DGP is centred near par_restr_hat,
# so both are excellent starting points). We add a handful of random restarts
# per model to avoid local optima distorting the LR upward.
n_restarts <- 8   # random restarts per model, per replicate (plus warm starts)

make_random_start <- function(base, lower, upper) {
  s <- base + rnorm(n_par, 0, 0.5)
  s[6:9] <- runif(4, 0, 0.15)
  # respect pinned bounds where lower==upper
  fixed <- which(is.finite(lower) & is.finite(upper) & (lower == upper))
  s[fixed] <- lower[fixed]
  s
}

fit_one <- function(objfun, C, start, lower, upper) {
  tryCatch(
    nlminb(start = start, objective = objfun,
           C_spec = C$spec, C_high = C$high, C_medium = C$medium, C_low = C$low,
           lower = lower, upper = upper,
           control = list(iter.max = 10000, eval.max = 20000, rel.tol = 1e-12, trace = 0)),
    error = function(e) list(par = rep(NA_real_, n_par), objective = Inf, convergence = NA_integer_))
}

# Fit from an explicit list of starting vectors; return the best fit object
# (not just the objective) so we can reuse the winning par as a downstream start.
best_fit <- function(objfun, C, starts, lower, upper) {
  best <- list(par = rep(NA_real_, n_par), objective = Inf)
  for (s in starts) {
    f <- fit_one(objfun, C, s, lower, upper)
    if (is.finite(f$objective) && f$objective < best$objective) best <- f
  }
  best
}

# Translate a RESTRICTED solution into a valid FULL-model starting vector.
# The full model reads pi via softmax(par[1:5]); the restricted solution has
# pi = 0 on drop_pi. We reproduce that by pushing the dropped raw-pi slots far
# below the kept ones (softmax -> ~0) while copying all other coordinates.
# The dropped motives' conditional blocks were pinned at 0 in the restricted
# fit; 0 is a fine interior start for the full model (softmax(0,..)=uniform).
# This guarantees the full model is offered a point that attains the restricted
# log-likelihood, so the full optimum can only be >= the restricted one.
restricted_to_full_start <- function(par_restr) {
  s <- par_restr
  s[drop_pi] <- min(s[keep_pi]) - 15   # pi[drop_pi] ~ 0 after softmax
  s
}

boot_worker <- function(b) {
  set.seed(10000 + b)            # reproducible per-replicate seed
  C <- simulate_all()

  # --- Restricted model first (data simulated under its MLE -> warm start nails it).
  restr_starts <- c(
    list(par_restr_hat),
    lapply(seq_len(n_restarts), function(j) make_random_start(par_restr_hat, lower_restr, upper_restr))
  )
  fit_restr <- best_fit(negloglik_restricted, C, restr_starts, lower_restr, upper_restr)
  negll_restr <- fit_restr$objective

  # --- Full model: seed it with (a) the full MLE, (b) the restricted solution
  #     re-expressed as a full-model point [GUARANTEES ll_full >= ll_restr],
  #     and (c) random restarts around the full MLE.
  full_seed_from_restr <- restricted_to_full_start(fit_restr$par)
  full_starts <- c(
    list(par_full_hat, full_seed_from_restr),
    lapply(seq_len(n_restarts), function(j) make_random_start(par_full_hat, lower_full, upper_full))
  )
  fit_full <- best_fit(negloglik, C, full_starts, lower_full, upper_full)
  negll_full <- fit_full$objective

  ll_full  <- -negll_full
  ll_restr <- -negll_restr
  lr <- 2 * (ll_full - ll_restr)
  # Numerical guard: full must fit at least as well as restricted.
  if (is.finite(lr) && lr < 0) lr <- 0

  # --- Estimated pi (5) and trembles (4) for each model ---
  # Full model: pi = softmax over all 5 raw slots.
  if (all(is.finite(fit_full$par))) {
    pi_full <- softmax(fit_full$par[1:5])
    tr_full <- fit_full$par[6:9]
  } else {
    pi_full <- rep(NA_real_, 5); tr_full <- rep(NA_real_, 4)
  }
  # Restricted model: pi reconstructed with dropped motives = exact 0.
  if (all(is.finite(fit_restr$par[keep_pi]))) {
    pi_restr <- numeric(5)
    pi_restr[keep_pi] <- softmax(fit_restr$par[keep_pi])
    pi_restr[drop_pi] <- 0
    tr_restr <- fit_restr$par[6:9]
  } else {
    pi_restr <- rep(NA_real_, 5); tr_restr <- rep(NA_real_, 4)
  }

  pi_lab <- c("self", "self_mob", "self_mer", "self_mob_mer", "rand")
  tr_lab <- c("spec", "high", "medium", "low")

  out <- c(lr = lr, ll_full = ll_full, ll_restr = ll_restr,
           pi_full, tr_full, pi_restr, tr_restr)
  names(out) <- c("lr", "ll_full", "ll_restr",
                  paste0("pi_full_",    pi_lab),
                  paste0("tr_full_",    tr_lab),
                  paste0("pi_restr_",   pi_lab),
                  paste0("tr_restr_",   tr_lab))
  out
}


#=======================================================================
# 4. Run the bootstrap in parallel
#=======================================================================
inizio  <- Sys.time()
n_boot  <- 1000
n_cores <- max(1, detectCores())

boot_list <- mclapply(seq_len(n_boot), boot_worker, mc.cores = n_cores, mc.preschedule = FALSE)
boot_mat  <- do.call(rbind, boot_list)   # cols: lr, ll_full, ll_restr, pi_*/tr_* for both models
LR_boot       <- boot_mat[, "lr"]
ll_full_boot  <- boot_mat[, "ll_full"]
ll_restr_boot <- boot_mat[, "ll_restr"]

cat("Finished", n_boot, "parametric bootstrap replicates using", n_cores, "cores.\n")
print(Sys.time() - inizio)


#=======================================================================
# 5. p-value and summary
#=======================================================================
ok <- is.finite(LR_boot)
LR_boot_ok <- LR_boot[ok]
n_ok <- length(LR_boot_ok)

# Bootstrap p-value with the standard +1 correction.
p_value <- (1 + sum(LR_boot_ok >= LR_obs)) / (1 + n_ok)

cat(sprintf("\nValid bootstrap replicates: %d / %d\n", n_ok, n_boot))
cat(sprintf("Observed LR: %.4f\n", LR_obs))
cat(sprintf("Bootstrap LR quantiles (90/95/99%%): %.3f / %.3f / %.3f\n",
            quantile(LR_boot_ok, 0.90), quantile(LR_boot_ok, 0.95), quantile(LR_boot_ok, 0.99)))
cat(sprintf("Parametric bootstrap p-value: %.4f\n", p_value))

# Per-replicate draws: LR, the two log-likelihoods, and the estimated
# pi (5) and trembles (4) for BOTH the full and restricted fits.
boot_draws <- data.frame(replicate = seq_len(n_boot), boot_mat,
                         row.names = NULL, check.names = FALSE)

results_lrt <- list(
  model         = "drop_mer",
  LR_obs        = LR_obs,
  ll_full_obs   = ll_full_obs,
  ll_restr_obs  = ll_restr_obs,
  p_value       = p_value,
  n_boot        = n_boot,
  n_valid       = n_ok,
  LR_boot       = LR_boot,
  ll_full_boot  = ll_full_boot,
  ll_restr_boot = ll_restr_boot,
  boot_mat      = boot_mat        # full per-replicate matrix incl. pi & trembles
)

write.csv(boot_draws, "bootstrap_LRT_drop_mer_draws.csv", row.names = FALSE)
saveRDS(results_lrt, "bootstrap_LRT_drop_mer.rds")
save(list = ls(), file = "bootstrap_LRT_drop_mer.Rdata")
#---------------------------------------
