#---------------------------------------
# Load libraries and prepare data
#---------------------------------------
source("stima_pre_process.R")
#---------------------------------------


#---------------------------------------
# ML estimation with 1000 starting points -- parallel version (macOS)
#---------------------------------------
inizio <- Sys.time()

set.seed(100)

n_starts <- 1000
n_par <- 89

lower <- rep(-Inf, n_par)
upper <- rep( Inf, n_par)

# constrain trembling parameters
lower[6:9] <- 0
upper[6:9] <- 0.5

# Uniform conditional probabilities
lower[10:n_par] <- 0
upper[10:n_par] <- 0

# Number of cores to use
n_cores <- max(1, detectCores() - 2)

make_start <- function() {
  s <- rnorm(n_par)
  s[6] <- runif(1, 0, 0.15)
  s[7] <- runif(1, 0, 0.15)
  s[8] <- runif(1, 0, 0.15)
  s[9] <- runif(1, 0, 0.15)
  s[10:n_par] <- 0
  s
}

# Pre-generate all starting values for reproducibility
start_list <- lapply(seq_len(n_starts), function(i) make_start())

multi_start_worker <- function(i) {
  start_i <- start_list[[i]]

  fit_i <- tryCatch(
    nlminb(
      start = start_i,
      objective = negloglik,
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
  fit_pi <- if (all(is.finite(fit_i$par[1:5]))) softmax(fit_i$par[1:5]) else rep(NA_real_, 5)

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
print(c(softmax(fit$par[1:5]), fit$par[6], fit$par[7], fit$par[8], fit$par[9], fit$objective))

fine <- Sys.time()
print(fine - inizio)

save(list = ls(), file = "stima_tremble_1000starts.Rdata")
#---------------------------------------
