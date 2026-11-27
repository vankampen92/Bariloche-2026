#!/usr/bin/env Rscript
# Logistic fishery: exact Gillespie simulations and deterministic approximation.
# Illustrates the model used in the course's third simple application.
# No external R packages are required.
#
# IMPORTANT DIFFERENCE BETWEEN THE TWO IMMIGRATION MODELS:
#  notes:    g_n = beta*n*(1-n/K) + m*(K-n)      (space-limited immigration)
#  constant: g_n = beta*n*(1-n/K) + I        (constant immigration, as in Cuenda et al.)
# With immigration_mode = "constant", the default I=m*K is chosen so
# the immigration rates agree at n=0; this is NOT a full calibration of
# all the model parameters in the Cuenda et al. paper.

## ---- EDIT PARAMETERS HERE --------------------------------------------------
beta <- 1.0             # per-capita recruitment coefficient, year^-1
delta <- 0.2            # per-capita natural death rate, year^-1
K <- 400L              # maximum number of fish (finite state space 0,...,K)
m <- 0.002             # immigrant arrival rate per vacant space, year^-1
I <- m * K             # constant-immigration alternative, fish/year
a <- 0.08              # attack coefficient, (fish * year)^-1 per boat
tau <- 0.2             # handling time (years)
boats_after <- 18.3     # fishing pressure from year 10 onward
t_switch <- 10.0        # fishing begins (years)
t_end <- 70.0           # duration (years)
n_sim <- 10L            # number of independent stochastic realizations
master_seed <- 2026L   # reproducibility
grid_step <- 0.05       # spacing for CSV output and deterministic integration
immigration_mode <- "notes"  # choose "notes" or "constant"
output_prefix <- "fishery_R"
## ---------------------------------------------------------------------------

stopifnot(immigration_mode %in% c("notes", "constant"), K > 0,
          t_switch > 0, t_switch < t_end, n_sim > 0, grid_step > 0)

# Treat the ecological state as a finite-capacity count (0 <= n <= K).
# At the upper boundary, no more recruitment events can occur.
rates <- function(n, boats) {
  immigration <- if (immigration_mode == "notes") m * (K - n) else I
  recruitment <- beta * n * (1 - n / K)
  g <- max(0, recruitment + immigration)
  if (n >= K) g <- 0
  r <- delta * n + a * boats * n / (1 + a * tau * n)
  if (n <= 0) r <- 0
  c(birth = g, death = r)
}

# Deterministic mean-field approximation. Keep the same rates, but allow
# continuous state values: the clipping at K is only a state-space boundary.
drift <- function(x, boats) {
  immigration <- if (immigration_mode == "notes") m * (K - x) else I
  beta * x * (1 - x / K) + immigration - delta * x -
    a * boats * x / (1 + a * tau * x)
}

# Initialize near the stable no-harvest equilibrium.
# For both default immigration forms the root lies in (0,K).
n_equilibrium <- uniroot(function(x) drift(x, 0), c(0, K))$root
n0 <- as.integer(round(n_equilibrium))

# One exact path, with an externally imposed change at t_switch.
# Switching harvesting on is NOT a stochastic population event.
gillespie_path <- function(seed) {
  set.seed(seed)
  t <- 0.0
  n <- n0
  times <- t
  states <- n
  while (t < t_end) {
    boats <- if (t < t_switch) 0 else boats_after
    rr <- rates(n, boats)
    total <- sum(rr)
    if (total <= 0) {
      if (t < t_switch) {
        t <- t_switch; times <- c(times,t); states <- c(states,n)
        next
      }
      break
    }
    dt <- rexp(1, rate = total)
    if (t < t_switch && t + dt >= t_switch) {
      t <- t_switch
      times <- c(times, t); states <- c(states, n)
      next
    }
    if (t + dt > t_end) break
    t <- t + dt
    n <- if (runif(1) < rr["birth"] / total) n + 1L else n - 1L
    times <- c(times, t)
    states <- c(states, n)
  }
  times <- c(times, t_end)
  states <- c(states, n)
  list(time = times, fish = states)
}

paths <- lapply(seq_len(n_sim), function(i) gillespie_path(master_seed + i - 1L))

# 4th-order Runge-Kutta method, with the change at t_switch handled exactly.
rk4 <- function(x, h, boats) {
  k1 <- drift(x, boats)
  k2 <- drift(x + h*k1/2, boats)
  k3 <- drift(x + h*k2/2, boats)
  k4 <- drift(x + h*k3, boats)
  x + h*(k1 + 2*k2 + 2*k3 + k4)/6
}
t_grid <- sort(unique(c(seq(0, t_end, by = grid_step), t_switch, t_end)))
t_grid <- t_grid[t_grid <= t_end]
det <- numeric(length(t_grid))
det[1] <- n0
for (j in 2:length(t_grid)) {
  t_left <- t_grid[j-1]
  h <- t_grid[j] - t_left
  det[j] <- rk4(det[j-1], h, if (t_left < t_switch) 0 else boats_after)
}

# Sample each right-continuous path at common times for convenient CSV use.
get_grid_state <- function(path, grid) {
  idx <- findInterval(grid, path$time, rightmost.closed = TRUE)
  path$fish[pmax(idx, 1L)]
}
samples <- sapply(paths, get_grid_state, grid = t_grid)
colnames(samples) <- paste0("simulation_", seq_len(n_sim))
out <- data.frame(year = t_grid, samples, deterministic = det, check.names = FALSE)
write.csv(out, paste0(output_prefix, "_trajectories.csv"), row.names = FALSE)

plot_simulations <- function() {
  cols <- grDevices::hcl.colors(n_sim, "Dark 3")
  plot(NA, xlim = c(0, t_end), ylim = c(0, K * 1.02),
       xlab = "Time (years)", ylab = "Number of fish",
       main = paste0("Logistic fishery: ", n_sim,
                     " stochastic trajectories vs deterministic dynamics"))
  abline(v = t_switch, col = "gray40", lty = 2, lwd = 1.5)
  for (i in seq_along(paths)) {
    lines(paths[[i]]$time, paths[[i]]$fish, type = "s",
          col = cols[i], lwd = 0.9)
  }
  lines(t_grid, det, col = "black", lwd = 3)
  text(t_switch / 2, K * 0.94, "No harvesting", cex = 0.9)
  text(t_switch + 1, K * 0.94,
       paste("Harvesting on: boats =", boats_after), adj = 0, cex = 0.9)
  legend("topright", legend = c(paste("Simulation", seq_len(n_sim)),
                                   "Deterministic"),
         col = c(cols, "black"), lwd = c(rep(1, n_sim), 3),
         ncol = 3, cex = 0.7, bty = "n")
}

pdf(paste0(output_prefix, ".pdf"), width = 12, height = 6.8)
plot_simulations()
dev.off()
png(paste0(output_prefix, ".png"), width = 1800, height = 1020, res = 150)
plot_simulations()
dev.off()

metadata <- c(
  paste("immigration_mode:", immigration_mode),
  paste("beta:", beta), paste("delta:", delta), paste("K:", K),
  paste("m:", m), paste("I:", I), paste("a:", a), paste("tau:", tau),
  paste("boats_after:", boats_after), paste("t_switch:", t_switch),
  paste("t_end:", t_end), paste("n_sim:", n_sim),
  paste("master_seed:", master_seed),
  paste("no-harvest deterministic equilibrium:", round(n_equilibrium, 4)),
  paste("initial abundance:", n0),
  paste("final deterministic abundance:", round(tail(det, 1), 3)),
  paste("final stochastic abundances:",
        paste(vapply(paths, function(z) tail(z$fish, 1), numeric(1)), collapse = ", "))
)
writeLines(metadata, paste0(output_prefix, "_parameters.txt"))
cat(paste(metadata, collapse = "\n"), "\n")
cat("Files: ", paste0(output_prefix, c(".png", ".pdf", "_trajectories.csv", "_parameters.txt")), "\n")
