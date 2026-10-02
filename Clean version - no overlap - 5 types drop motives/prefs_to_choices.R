# ----------------------------------------------------------------------------
# Evaluate binary choices against a string of preference relations.
#
# Notation:
#   ">"  -> strictly preferred to
#   "="  -> indifferent (same equivalence class)
#   "-"  -> joins parts of one preference relation. Parts that share items
#           are merged into a single partial order through the shared items.
#   ";"  -> separates different preference relations (e.g. respondents)
#
# A binary choice (a, b) returns:
#   * 0 if a (first) is preferred to b under the merged partial order,
#   * 1 if b (second) is preferred to a,
#   * NA if a and b are in the same indifference class, or
#   * NA if they remain incomparable after merging (no path links them).
# ----------------------------------------------------------------------------
 
 
# Parse one preference relation into a partial-order structure:
#   $reps      : named character vector item -> equivalence-class representative
#   $reachable : named list class -> classes strictly less preferred than it
parse_preference <- function(pref_str) {
  pref_str <- gsub("\\s+", "", pref_str)
  parts <- strsplit(pref_str, "-", fixed = TRUE)[[1]]
 
  # Parse each part into a list of tiers; each tier = character vector of items
  parsed_parts <- lapply(parts, function(part) {
    tiers <- strsplit(part, ">", fixed = TRUE)[[1]]
    lapply(tiers, function(tier) strsplit(tier, "=", fixed = TRUE)[[1]])
  })
 
  all_items <- unique(unlist(parsed_parts))
  if (length(all_items) == 0) return(list(reps = character(0), reachable = list()))
 
  # ---- Union-Find: group items declared equivalent (same tier anywhere) ----
  parent <- setNames(all_items, all_items)
  find <- function(x) {
    while (parent[[x]] != x) {
      parent[[x]] <<- parent[[parent[[x]]]]   # path compression
      x <- parent[[x]]
    }
    x
  }
  do_union <- function(x, y) {
    rx <- find(x); ry <- find(y)
    if (rx != ry) parent[[rx]] <<- ry
  }
  for (part in parsed_parts) {
    for (tier in part) {
      if (length(tier) > 1) {
        for (k in seq.int(2, length(tier))) do_union(tier[1], tier[k])
      }
    }
  }
 
  reps <- vapply(all_items, find, character(1))
  names(reps) <- all_items
  classes <- unique(reps)
 
  # ---- Direct preference edges between equivalence classes -----------------
  adj <- setNames(lapply(classes, function(x) character(0)), classes)
  for (part in parsed_parts) {
    nt <- length(part)
    if (nt < 2) next
    tier_reps <- vapply(part, function(tier) reps[[tier[1]]], character(1))
    for (i in seq_len(nt - 1)) {
      for (j in seq.int(i + 1, nt)) {
        from <- tier_reps[i]; to <- tier_reps[j]
        if (from != to && !(to %in% adj[[from]])) {
          adj[[from]] <- c(adj[[from]], to)
        }
      }
    }
  }
 
  # ---- Transitive closure: BFS from every class ----------------------------
  reachable <- setNames(lapply(classes, function(start) {
    seen <- character(0)
    queue <- adj[[start]]
    while (length(queue) > 0) {
      n <- queue[1]; queue <- queue[-1]
      if (!(n %in% seen)) {
        seen <- c(seen, n)
        queue <- c(queue, adj[[n]])
      }
    }
    seen
  }), classes)
 
  list(reps = reps, reachable = reachable)
}
 
 
# Decide which of two items is preferred under a parsed relation.
# Returns 0L if `a` (first) is preferred, 1L if `b` (second) is preferred,
# and NA_integer_ if indifferent or incomparable.
compare_items <- function(parsed_pref, a, b) {
  a <- as.character(a); b <- as.character(b)
  reps      <- parsed_pref$reps
  reachable <- parsed_pref$reachable
 
  if (!(a %in% names(reps)) || !(b %in% names(reps))) return(NA_integer_)
  ca <- reps[[a]]; cb <- reps[[b]]
  if (ca == cb)                     return(NA_integer_)   # indifferent
  if (cb %in% reachable[[ca]])      return(0L)            # a > b
  if (ca %in% reachable[[cb]])      return(1L)            # b > a
  NA_integer_                                              # incomparable
}
 
 
# Main entry point.
# pref_string    : a character string with one or more relations split by ";"
# binary_choices : a list of length-2 vectors, e.g. list(c(1, 7), c(2, 3))
# Returns a named integer vector, one entry per (relation, binary choice):
# 0 = first item of the pair chosen, 1 = second item chosen, NA = indifferent.
evaluate_choices <- function(pref_string, binary_choices) {
  relations   <- trimws(strsplit(pref_string, ";", fixed = TRUE)[[1]])
  parsed_list <- lapply(relations, parse_preference)
 
  result <- matrix(
    NA_integer_,
    nrow = length(parsed_list),
    ncol = length(binary_choices),
    dimnames = list(
      paste0("Pref_", seq_along(parsed_list)),
      vapply(binary_choices,
             function(bc) paste0(bc[1], " vs ", bc[2]),
             character(1))
    )
  )
 
  for (i in seq_along(parsed_list)) {
    for (j in seq_along(binary_choices)) {
      bc <- binary_choices[[j]]
      result[i, j] <- compare_items(parsed_list[[i]], bc[1], bc[2])
    }
  }
 
  # Flatten the matrix into a single vector
  as.numeric(c(t(result)))
}
 
 
# ---- Prefs by status ---------------------------------------------------------
prefs_high <- c(
  "1=4=7 - 2=3;  1=4=7 - 2=3;  1=4=7 - 2=3", #Rand
  "1>4>7 - 2>3;  1>4>7 - 2>3;  1>4>7 - 2>3", #Self
  "1>4>7 - 2>3;  1>4>7 - 2>3;  1>7>4 - 2>3", #Self_Mer
  "1>4>7 - 2>3;  1>4>7 - 2>3;  7>1>4 - 2>3", #Self_Mer
  "1>4>7 - 2>3;  1>4>7 - 2>3;  4>1>7 - 2>3", #Self_Mer
  "1>4>7 - 2>3;  1>4>7 - 2>3;  4>7>1 - 2>3", #Self_Mer
  "1>4>7 - 2>3;  1>4>7 - 2>3;  7>4>1 - 2>3", #Self_Mer
  "1>7>4 - 2>3;  1>7>4 - 2>3;  1>7>4 - 2>3", #Self_Mob
  "1>7>4 - 2>3;  1>7>4 - 2>3;  7>1>4 - 2>3", #Self_Mob_Mer
  "1>7>4 - 2>3;  1>7>4 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "1>4>7 - 2>3;  1>7>4 - 2>3;  1>7>4 - 2>3", #Self_Mob_Mer
  "1>4>7 - 2>3;  1>7>4 - 2>3;  7>1>4 - 2>3", #Self_Mob_Mer
  "1>4>7 - 2>3;  1>7>4 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "7>1>4 - 2>3;  7>1>4 - 2>3;  7>1>4 - 2>3", #Self_Mob
  "7>1>4 - 2>3;  7>1>4 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "1>4>7 - 2>3;  7>1>4 - 2>3;  7>1>4 - 2>3", #Self_Mob_Mer
  "1>4>7 - 2>3;  7>1>4 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "1>7>4 - 2>3;  7>1>4 - 2>3;  7>1>4 - 2>3", #Self_Mob_Mer
  "1>7>4 - 2>3;  7>1>4 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "7>4>1 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mob
  "1>4>7 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "1>7>4 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "7>1>4 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "4>1>7 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "4>7>1 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "4>1>7 - 2>3;  4>1>7 - 2>3;  4>1>7 - 2>3", #Self_Mob
  "1>4>7 - 2>3;  4>1>7 - 2>3;  4>1>7 - 2>3", #Self_Mob_Mer
  "4>1>7 - 2>3;  4>1>7 - 2>3;  4>7>1 - 2>3", #Self_Mob_Mer
  "1>4>7 - 2>3;  4>1>7 - 2>3;  4>7>1 - 2>3", #Self_Mob_Mer
  "4>1>7 - 2>3;  4>1>7 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "1>4>7 - 2>3;  4>1>7 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "4>7>1 - 2>3;  4>7>1 - 2>3;  4>7>1 - 2>3", #Self_Mob
  "1>4>7 - 2>3;  4>7>1 - 2>3;  4>7>1 - 2>3", #Self_Mob_Mer
  "4>1>7 - 2>3;  4>7>1 - 2>3;  4>7>1 - 2>3", #Self_Mob_Mer
  "4>7>1 - 2>3;  4>7>1 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "1>4>7 - 2>3;  4>7>1 - 2>3;  7>4>1 - 2>3", #Self_Mob_Mer
  "4>1>7 - 2>3;  4>7>1 - 2>3;  7>4>1 - 2>3"  #Self_Mob_Mer
)

prefs_low <- c(
  "1=4=7 - 2=3;  1=4=7 - 2=3;  1=4=7 - 2=3", #Rand
  "7>4>1 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self + Self_Mob
  "1>4>7 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mer + Self_Mob_Mer 
  "1>7>4 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mer + Self_Mob_Mer
  "7>1>4 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mer + Self_Mob_Mer
  "4>1>7 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3", #Self_Mer + Self_Mob_Mer
  "4>7>1 - 2>3;  7>4>1 - 2>3;  7>4>1 - 2>3"  #Self_Mer + Self_Mob_Mer
)

prefs_medium <- c(
  "1=4=7 - 2=3;  1=4=7 - 2=3;  1=4=7 - 2=3", #Rand
  "3>7>2 - 7=1;  3>7>2 - 7=1;  3>7>2 - 7=1", #Self
  "3>7>2 - 1>7;  3>7>2 - 7=1;  3>7>2 - 7>1", #Self_Mer
  "3>7>2 - 1>7;  3>7>2 - 7=1;  7>3>2 - 7>1", #Self_Mer
  "3>2>7 - 1>7;  3>7>2 - 7=1;  3>7>2 - 7>1", #Self_Mer
  "3>2>7 - 1>7;  3>7>2 - 7=1;  7>3>2 - 7>1", #Self_Mer
  "3>7>2 - 7>1;  3>7>2 - 7>1;  3>7>2 - 7>1", #Self_Mob
  "3>7>2 - 1>7;  3>7>2 - 7>1;  3>7>2 - 7>1", #Self_Mob_Mer
  "3>2>7 - 7>1;  3>7>2 - 7>1;  3>7>2 - 7>1", #Self_Mob_Mer
  "3>2>7 - 1>7;  3>7>2 - 7>1;  3>7>2 - 7>1", #Self_Mob_Mer
  "3>7>2 - 7>1;  3>7>2 - 7>1;  7>3>2 - 7>1", #Self_Mob_Mer
  "3>7>2 - 1>7;  3>7>2 - 7>1;  7>3>2 - 7>1", #Self_Mob_Mer
  "3>2>7 - 7>1;  3>7>2 - 7>1;  7>3>2 - 7>1", #Self_Mob_Mer
  "3>2>7 - 1>7;  3>7>2 - 7>1;  7>3>2 - 7>1", #Self_Mob_Mer
  "7>3>2 - 7>1;  7>3>2 - 7>1;  7>3>2 - 7>1", #Self_Mob
  "3>2>7 - 7>1;  7>3>2 - 7>1;  7>3>2 - 7>1", #Self_Mob_Mer
  "3>2>7 - 1>7;  7>3>2 - 7>1;  7>3>2 - 7>1", #Self_Mob_Mer
  "3>7>2 - 1>7;  7>3>2 - 7>1;  7>3>2 - 7>1", #Self_Mob_Mer
  "3>7>2 - 7>1;  7>3>2 - 7>1;  7>3>2 - 7>1", #Self_Mob_Mer
  "7>3>2 - 1>7;  7>3>2 - 7>1;  7>3>2 - 7>1"  #Self_Mob_Mer
)

prefs_spectator <- c(
  "7=4=1 - 2=3;  7=4=1 - 2=3;  7=4=1 - 2=3", #Rand + Self
  "1>4>7 - 2=3;  7=4=1 - 2=3;  7>4>1 - 2=3", #Self_Mer
  "1>4>7 - 2=3;  7>4>1 - 2=3;  7>4>1 - 2=3", #Self_Mob_Mer
  "1>7>4 - 2=3;  7>4>1 - 2=3;  7>4>1 - 2=3", #Self_Mob_Mer
  "7>1>4 - 2=3;  7>4>1 - 2=3;  7>4>1 - 2=3", #Self_Mob_Mer
  "4>1>7 - 2=3;  7>4>1 - 2=3;  7>4>1 - 2=3", #Self_Mob_Mer
  "4>7>1 - 2=3;  7>4>1 - 2=3;  7>4>1 - 2=3", #Self_Mob_Mer
  "7>4>1 - 2=3;  7>4>1 - 2=3;  7>4>1 - 2=3"  #Self_Mob
)

