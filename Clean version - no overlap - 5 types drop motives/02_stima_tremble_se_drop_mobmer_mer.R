#---------------------------------------
# Load libraries and prepare data
#---------------------------------------
source("stima_pre_process.R")
#---------------------------------------


#---------------------------------------
# RESTRICTED MODEL B: drop Self_Mob_Mer AND Self_Mer
# (see 01_stima_tremble_drop_mobmer_mer_1000starts_parallel.R for the trick)
#---------------------------------------
drop_pi <- c(3, 4)   # Self_Mer, Self_Mob_Mer
keep_pi <- c(1, 2, 5)

pin_idx <- c(10:14,   # spec   Self_Mob_Mer (a)
             15:19,   # low    Self_Mer + Self_Mob_Mer (shared a)
             20:23,   # medium Self_Mer (a)
             38:42,   # high   Self_Mer (a)
             26:37,   # medium Self_Mob_Mer (c)
             48:72)   # high   Self_Mob_Mer (c)

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


#---------------------------------------
# Load best model from the restricted 1000-start estimation
#---------------------------------------
inizio <- Sys.time()

load("stima_tremble_drop_mobmer_mer_1000starts.Rdata")

if (exists("best_par") && all(is.finite(best_par))) {
  best_model_par <- best_par
} else if (exists("fit") && !is.null(fit$par) && all(is.finite(fit$par))) {
  best_model_par <- fit$par
} else {
  stop("Could not find a valid best parameter vector in stima_tremble_drop_mobmer_mer_1000starts.Rdata")
}

best_objective <- if (exists("fit") && !is.null(fit$objective)) fit$objective else negloglik_restricted(
  best_model_par,
  dati_stima_spec,
  dati_stima_high,
  dati_stima_medium,
  dati_stima_low
)

best_convergence <- if (exists("fit") && !is.null(fit$convergence)) fit$convergence else NA_integer_
best_message <- if (exists("fit") && !is.null(fit$message)) fit$message else NA_character_
best_start_id <- if (exists("fit") && !is.null(fit$start_id)) fit$start_id else NA_integer_

fit <- list(
  par = best_model_par,
  objective = best_objective,
  convergence = best_convergence,
  message = best_message,
  start_id = best_start_id
)

cat("Loaded best model from stima_tremble_drop_mobmer_mer_1000starts.Rdata\n")
cat("Best start:", fit$start_id, "\n")
cat("Best objective:", fit$objective, "\n")
cat("Convergence code:", fit$convergence, "\n")
cat("Message:", fit$message, "\n")
pi_best <- numeric(5); pi_best[keep_pi] <- softmax(fit$par[keep_pi]); pi_best[drop_pi] <- 0
print(c(pi_best, fit$par[6], fit$par[7], fit$par[8], fit$par[9], fit$objective))
#---------------------------------------


#---------------------------------------
# ML bootstrap SE -- parallel version (macOS)
#---------------------------------------
inizio <- Sys.time()

# --- Bootstrap Setup ---
n_boot <- 1000
set.seed(123)

n_par <- 72

# Number of participants in each sample
n_spec   <- dim(dati_stima_spec)[1]
n_high   <- dim(dati_stima_high)[1]
n_medium <- dim(dati_stima_medium)[1]
n_low    <- dim(dati_stima_low)[1]

# Use best 1000-start estimate as warm start
start_boot <- fit$par

lower <- rep(-Inf, n_par)
upper <- rep( Inf, n_par)

# constrain trembling parameters
lower[6:9] <- 0
upper[6:9] <- 0.5

# pin the dropped motives' conditional parameters to 0
lower[pin_idx] <- 0
upper[pin_idx] <- 0

# Number of cores to use
n_cores <- max(1, detectCores() - 2)

# Pre-generate bootstrap indices for reproducibility
boot_indices <- lapply(seq_len(n_boot), function(b) {
  list(
    idx_spec   = sample.int(n_spec,   n_spec,   replace = TRUE),
    idx_high   = sample.int(n_high,   n_high,   replace = TRUE),
    idx_medium = sample.int(n_medium, n_medium, replace = TRUE),
    idx_low    = sample.int(n_low,    n_low,    replace = TRUE)
  )
})

boot_worker <- function(b) {
  idx <- boot_indices[[b]]

  boot_spec   <- dati_stima_spec[idx$idx_spec,   , , drop = FALSE]
  boot_high   <- dati_stima_high[idx$idx_high,   , , drop = FALSE]
  boot_medium <- dati_stima_medium[idx$idx_medium, , , drop = FALSE]
  boot_low    <- dati_stima_low[idx$idx_low,    , , drop = FALSE]

  boot_fit <- tryCatch(
    nlminb(
      start     = start_boot,
      objective = negloglik_restricted,
      C_spec    = boot_spec,
      C_high    = boot_high,
      C_medium  = boot_medium,
      C_low     = boot_low,
      lower     = lower,
      upper     = upper,
      control   = list(
        iter.max = 10000,
        eval.max = 20000,
        rel.tol  = 1e-12,
        trace    = 0
      )
    ),
    error = function(e) NULL
  )

  if (!is.null(boot_fit) && boot_fit$convergence %in% c(0, 1)) {
    pi_b <- numeric(5)
    pi_b[keep_pi] <- softmax(boot_fit$par[keep_pi])
    pi_b[drop_pi] <- 0
    c(
      pi_b,                # pi (dropped motives = exact 0)
      boot_fit$par[6],     # tremble_spec
      boot_fit$par[7],     # tremble_high
      boot_fit$par[8],     # tremble_medium
      boot_fit$par[9]      # tremble_low
    )
  } else {
    rep(NA_real_, 9)
  }
}

# Run bootstrap in parallel
boot_list <- mclapply(
  X        = seq_len(n_boot),
  FUN      = boot_worker,
  mc.cores = n_cores
)

boot_results <- do.call(rbind, boot_list)

cat("Finished", n_boot, "bootstrap iterations using", n_cores, "cores.\n")

fine <- Sys.time()
print(fine - inizio)
#---------------------------------------


#---------------------------------------
# Prepare tables with estimates and SE
#---------------------------------------
pi_hat <- numeric(5); pi_hat[keep_pi] <- softmax(fit$par[keep_pi]); pi_hat[drop_pi] <- 0
estimates <- c(
  pi_hat,
  fit$par[6],
  fit$par[7],
  fit$par[8],
  fit$par[9]
)

param_names <- c(
  paste0("pi_", 1:5, c("_self","_self_mob","_self_mer","_self_mob_mer","_rand")),
  "tremble_spec",
  "tremble_high",
  "tremble_medium",
  "tremble_low"
)

boot_se <- apply(boot_results, 2, sd, na.rm = TRUE)

results_tbl <- tibble(
  parameter = param_names,
  estimate = estimates,
  se = boot_se
)

results_tbl <- rbind(
  results_tbl,
  c("best_start_id", fit$start_id, NA),
  c("Log_lik", fit$objective, NA),
  c("convergence", fit$convergence, NA)
)

write.csv(results_tbl, "results_tremble_se_drop_mobmer_mer.csv", row.names = FALSE)
write.csv(as.data.frame(boot_results), "bootstrap_tremble_se_drop_mobmer_mer.csv", row.names = FALSE)

library(gt)
gt::gtsave(gt::gt(results_tbl), "results_tremble_se_drop_mobmer_mer.pdf")
#---------------------------------------


#---------------------------------------
# Compute empirical posterior and save them
#---------------------------------------
posterior_types <- function(par, C_spec, C_high, C_medium, C_low) {
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

  m_comp_spec   <- prob_profile_given_type_spec(pmt_spec, tremble_lik(tremble_spec, C_spec))
  m_comp_low    <- prob_profile_given_type_low(pmt_low, tremble_lik(tremble_low, C_low))
  m_comp_medium <- prob_profile_given_type_medium(pmt_medium, tremble_lik(tremble_medium, C_medium))
  m_comp_high   <- prob_profile_given_type_high(pmt_high, tremble_lik(tremble_high, C_high))

  log_pi <- log(pi)
  log_terms_spec   <- sweep(log(m_comp_spec), 2, log_pi, "+")
  log_terms_low    <- sweep(log(m_comp_low), 2, log_pi, "+")
  log_terms_medium <- sweep(log(m_comp_medium), 2, log_pi, "+")
  log_terms_high   <- sweep(log(m_comp_high), 2, log_pi, "+")

  log_denom_spec   <- apply(log_terms_spec, 1, logsumexp)
  log_denom_low    <- apply(log_terms_low, 1, logsumexp)
  log_denom_medium <- apply(log_terms_medium, 1, logsumexp)
  log_denom_high   <- apply(log_terms_high, 1, logsumexp)

  post_spec   <- exp(log_terms_spec - log_denom_spec)
  post_low    <- exp(log_terms_low - log_denom_low)
  post_medium <- exp(log_terms_medium - log_denom_medium)
  post_high   <- exp(log_terms_high - log_denom_high)

  colnames(post_spec)   <- paste0("type_", 1:5)
  colnames(post_low)    <- paste0("type_", 1:5)
  colnames(post_medium) <- paste0("type_", 1:5)
  colnames(post_high)   <- paste0("type_", 1:5)

  list(
    spec   = post_spec,
    low    = post_low,
    medium = post_medium,
    high   = post_high
  )
}

post <- posterior_types(
  par      = fit$par,
  C_spec   = dati_stima_spec,
  C_high   = dati_stima_high,
  C_medium = dati_stima_medium,
  C_low    = dati_stima_low
)

write.csv(data.frame(id_high, post$high), "posterior_high_tr_drop_mobmer_mer.csv", row.names = FALSE)
write.csv(data.frame(id_medium, post$medium), "posterior_med_tr_drop_mobmer_mer.csv", row.names = FALSE)
write.csv(data.frame(id_low, post$low), "posterior_low_tr_drop_mobmer_mer.csv", row.names = FALSE)
write.csv(data.frame(id_spec, post$spec), "posterior_spec_tr_drop_mobmer_mer.csv", row.names = FALSE)

save(list = ls(), file = "stima_tremble_se_drop_mobmer_mer.Rdata")
#---------------------------------------
