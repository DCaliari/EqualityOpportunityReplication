#---------------------------------------
# Load libraries and prepare data
#---------------------------------------
source("stima_pre_process.R")
#---------------------------------------


#---------------------------------------
# RESTRICTED MODEL A: drop Self_Mob_Mer AND Self_Mob
#
# Surviving motives: Self (pi[1]), Self_Mer (pi[3]), Rand (pi[5]).
# Dropped motives:   Self_Mob (pi[2]) and Self_Mob_Mer (pi[4]) -> pi fixed at 0.
#
# Trick (same as the project's no_self_mob script):
#   * pi is reparametrised so that the dropped motives are EXACTLY 0
#     (softmax can never return an exact 0; reparametrisation can).
#       pi[c(1,3,5)] = softmax(par[c(1,3,5)]) ;  pi[c(2,4)] = 0
#     We keep the full 72-vector so NOTHING downstream needs renumbering;
#     par[2] and par[4] simply become inert (their value is ignored).
#   * The conditional profile-distribution parameters belonging to the two
#     dropped motives are unidentified once their pi = 0, so we pin them to 0
#     via the bounds (and force the random starts to respect that).
#
#   IMPORTANT: which raw block feeds which motive is context-specific.
#   From the M matrices in stima_pre_process.R:
#     - spec  a (par[10:14]) fills the Self_Mob_Mer column  (Self_Mer is a
#       deterministic column, no free par)
#     - low   a (par[15:19]) fills BOTH Self_Mer AND Self_Mob_Mer (shared) ***
#     - medium a (par[20:23]) -> Self_Mer ; b (par[24:25]) -> Self_Mob ;
#               c (par[26:37]) -> Self_Mob_Mer
#     - high   a (par[38:42]) -> Self_Mer ; b (par[43:47]) -> Self_Mob ;
#               c (par[48:72]) -> Self_Mob_Mer
#
#   Here we drop Self_Mob and Self_Mob_Mer. So we pin:
#       Self_Mob_Mer : spec par[10:14], medium par[26:37], high par[48:72]
#       Self_Mob     : medium par[24:25], high par[43:47]
#   *** We do NOT pin low par[15:19]: in low that block is SHARED with the
#       surviving Self_Mer column, so it is still identified through Self_Mer
#       and must stay free. (Numerically the dropped Self_Mob_Mer column in low
#       just gets pi = 0 and drops out of the likelihood; the same a still
#       parametrises the kept Self_Mer column.)
#
# Effective free parameters:
#   pinned conditionals = 5 (spec a) + 12 (medium c) + 25 (high c)
#                        + 2 (medium b) + 5 (high b) = 49 ; inert raw pi = 2
#   free = 72 - 2 - 49 = 21
#---------------------------------------

# Indices of the dropped motives in the pi vector
drop_pi <- c(2, 4)   # Self_Mob, Self_Mob_Mer
keep_pi <- c(1, 3, 5)

# Conditional-parameter indices to pin at 0.
# (low par[15:19] deliberately NOT pinned: shared with surviving Self_Mer.)
pin_idx <- c(10:14,   # spec   Self_Mob_Mer (a)
             26:37,   # medium Self_Mob_Mer (c)
             48:72,   # high   Self_Mob_Mer (c)
             24:25,   # medium Self_Mob (b)
             43:47)   # high   Self_Mob (b)

# Restricted negative log-likelihood: same 72-vector interface as negloglik(),
# but pi is forced to 0 on the dropped motives via reparametrisation.
negloglik_restricted <- function(par, C_spec, C_high, C_medium, C_low) {
  pi <- numeric(5)
  pi[keep_pi] <- softmax(par[keep_pi])
  pi[drop_pi] <- 0                       # exact zeros

  tremble_spec   <- par[6]
  tremble_high   <- par[7]
  tremble_medium <- par[8]
  tremble_low    <- par[9]

  pmt_spec   <- softmax(par[10:14])                                           # 5 par
  pmt_low    <- softmax(par[15:19])                                           # 5 par
  pmt_medium <- c(softmax(par[20:23]), softmax(par[24:25]), softmax(par[26:37])) # 18 par
  pmt_high   <- c(softmax(par[38:42]), softmax(par[43:47]), softmax(par[48:72])) # 35 par

  m_comp_spec   <- prob_profile_given_type_spec(pmt_spec,   tremble_lik(tremble_spec,   C_spec))
  m_comp_low    <- prob_profile_given_type_low(pmt_low,     tremble_lik(tremble_low,    C_low))
  m_comp_medium <- prob_profile_given_type_medium(pmt_medium, tremble_lik(tremble_medium, C_medium))
  m_comp_high   <- prob_profile_given_type_high(pmt_high,   tremble_lik(tremble_high,   C_high))

  # log(pi)=-Inf on dropped cols is fine: logsumexp drops them in exp-space.
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
# ML estimation with 1000 starting points -- parallel version (macOS)
#---------------------------------------
inizio <- Sys.time()

set.seed(100)

n_starts <- 1000
n_par <- 72

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

make_start <- function() {
  s <- rnorm(n_par)
  s[6] <- runif(1, 0, 0.15)
  s[7] <- runif(1, 0, 0.15)
  s[8] <- runif(1, 0, 0.15)
  s[9] <- runif(1, 0, 0.15)
  s[pin_idx] <- 0          # respect the pinned-at-0 bounds
  s
}

# Pre-generate all starting values for reproducibility
start_list <- lapply(seq_len(n_starts), function(i) make_start())

multi_start_worker <- function(i) {
  start_i <- start_list[[i]]

  fit_i <- tryCatch(
    nlminb(
      start = start_i,
      objective = negloglik_restricted,
      C_spec = dati_stima_spec,
      C_high = dati_stima_high,
      C_medium = dati_stima_medium,
      C_low = dati_stima_low,
      lower = lower,
      upper = upper,
      control = list(
        iter.max = 10000,
        eval.max = 20000,
        rel.tol  = 1e-12,
        trace    = 0
      )
    ),
    error = function(e) {
      list(
        par = rep(NA_real_, n_par),
        objective = Inf,
        convergence = NA_integer_,
        message = paste("error:", conditionMessage(e))
      )
    }
  )

  fit_message <- if (is.null(fit_i$message)) NA_character_ else as.character(fit_i$message)
  # report pi with the dropped motives shown as exact 0
  if (all(is.finite(fit_i$par[keep_pi]))) {
    fit_pi <- numeric(5)
    fit_pi[keep_pi] <- softmax(fit_i$par[keep_pi])
    fit_pi[drop_pi] <- 0
  } else {
    fit_pi <- rep(NA_real_, 5)
  }

  list(
    summary = tibble(
      start_id = i,
      objective = fit_i$objective,
      convergence = fit_i$convergence,
      message = fit_message,
      pi_self = fit_pi[1],
      pi_self_mob = fit_pi[2],
	  pi_self_mer = fit_pi[3],
	  pi_self_mob_mer = fit_pi[4],
      pi_rand = fit_pi[5],
      tremble_spec = fit_i$par[6],
      tremble_high = fit_i$par[7],
      tremble_medium = fit_i$par[8],
      tremble_low = fit_i$par[9]
    ),
    par = fit_i$par
  )
}

# Run multi-start optimization in parallel
multi_start_list <- mclapply(
  X = seq_len(n_starts),
  FUN = multi_start_worker,
  mc.cores = n_cores
)

# Collect results
starts_tbl <- bind_rows(lapply(multi_start_list, `[[`, "summary")) %>%
  arrange(objective)
all_par <- lapply(multi_start_list, `[[`, "par")

cat("Finished", n_starts, "starting points using", n_cores, "cores.\n")

valid_idx <- which(is.finite(starts_tbl$objective))
if (length(valid_idx) == 0) {
  stop("All optimization runs failed.")
}

best_start_id <- starts_tbl$start_id[1]
best_par <- all_par[[best_start_id]]

fit <- list(
  par = best_par,
  objective = starts_tbl$objective[1],
  convergence = starts_tbl$convergence[1],
  message = starts_tbl$message[1],
  start_id = best_start_id
)

cat("Best start:", fit$start_id, "\n")
cat("Best objective:", fit$objective, "\n")
cat("Convergence code:", fit$convergence, "\n")
cat("Message:", fit$message, "\n")
pi_best <- numeric(5); pi_best[keep_pi] <- softmax(fit$par[keep_pi]); pi_best[drop_pi] <- 0
print(c(pi_best, fit$par[6], fit$par[7], fit$par[8], fit$par[9], fit$objective))

fine <- Sys.time()
print(fine - inizio)

save(list = ls(), file = "stima_tremble_drop_mobmer_mob_1000starts.Rdata")
#---------------------------------------
