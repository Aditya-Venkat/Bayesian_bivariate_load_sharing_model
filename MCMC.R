#-------------------------------------------------------------------------------
# 1. Data Simulation (Corrected)
#-------------------------------------------------------------------------------
draw_biv_data <- function(n, lambdas, alpha) {
  out <- matrix(nrow = n, ncol = 2)
  u1 <- runif(n, 0, 1)
  u2 <- runif(n, 0, 1)
  u3 <- runif(n, 0, 1)
  
  v1 <- -log(u1^(1/lambdas[1])) / alpha
  v2 <- -log(u2^(1/lambdas[2])) / alpha
  v3 <- -log(u3^(1/lambdas[3])) / alpha
  
  x <- pmin(v1, v3)
  y <- pmin(v2, v3)
  
  final_draw <- function(v3_i, x_i, y_i, lambdas, alpha) {
    if (x_i == y_i) {
      return(c(x_i, y_i))
    } else if (x_i > y_i) {
      y2 <- y_i
      u <- runif(1, 0, 1)
      y1 <- -log(u^(1/lambdas[4]) * exp(-alpha * y2)) / alpha
      return(c(y1, y2))
    } else {
      y1 <- x_i
      u <- runif(1, 0, 1)
      y2 <- -log(u^(1/lambdas[5]) * exp(-alpha * y1)) / alpha
      return(c(y1, y2))
    }
  }
  
  out <- t(sapply(1:n, function(i) final_draw(v3[i], x[i], y[i], lambdas, alpha)))
  return(out)
}

#-------------------------------------------------------------------------------
# 2. Log-Posterior Calculation (Corrected Indexing)
#-------------------------------------------------------------------------------
log_posterior <- function(lambdas, data, hyperparams) {
  set1 <- data$set1  # I1: y1 = y2
  set2 <- data$set2  # I2: y1 < y2
  set3 <- data$set3  # I3: y1 > y2
  n1 <- nrow(set1)
  n2 <- nrow(set2)
  n3 <- nrow(set3)
  
  # Compute S(lambda) (Eq. 4.2 with corrected indexing)
  S <- (lambdas[1] + lambdas[2] + lambdas[3]) * sum(set1[, 1]) +
    (lambdas[1] + lambdas[2] + lambdas[3] - lambdas[5]) * sum(set2[, 1]) +
    lambdas[5] * sum(set2[, 2]) +
    (lambdas[1] + lambdas[2] + lambdas[3] - lambdas[4]) * sum(set3[, 2]) +
    lambdas[4] * sum(set3[, 1])
  
  # Gamma prior terms (corrected n_j assignments)
  log_prior <- sum(
    (n2 + hyperparams$cs[1] - 1) * log(lambdas[1]) - hyperparams$ds[1] * lambdas[1],
    (n3 + hyperparams$cs[2] - 1) * log(lambdas[2]) - hyperparams$ds[2] * lambdas[2],
    (n1 + hyperparams$cs[3] - 1) * log(lambdas[3]) - hyperparams$ds[3] * lambdas[3],
    (n3 + hyperparams$cs[4] - 1) * log(lambdas[4]) - hyperparams$ds[4] * lambdas[4],
    (n2 + hyperparams$cs[5] - 1) * log(lambdas[5]) - hyperparams$ds[5] * lambdas[5]
  )
  
  log_likelihood <- -(hyperparams$a + n1 + 2*n2 + 2*n3) * log(hyperparams$b + S)
  
  return(log_prior + log_likelihood)
}

#-------------------------------------------------------------------------------
# 3. Gradient Computation (Fully Corrected)
#-------------------------------------------------------------------------------
compute_gradient <- function(lambda, data, hyperparams) {
  set1 <- data$set1  # I1: y1 = y2
  set2 <- data$set2  # I2: y1 < y2
  set3 <- data$set3  # I3: y1 > y2
  n1 <- nrow(set1)
  n2 <- nrow(set2)
  n3 <- nrow(set3)
  
  # Compute S and its gradient
  S <- (lambda[1] + lambda[2] + lambda[3]) * sum(set1[, 1]) +
    (lambda[1] + lambda[2] + lambda[3] - lambda[5]) * sum(set2[, 1]) +
    lambda[5] * sum(set2[, 2]) +
    (lambda[1] + lambda[2] + lambda[3] - lambda[4]) * sum(set3[, 2]) +
    lambda[4] * sum(set3[, 1])
  
  # Partial derivatives of S
  dS_dlambda <- c(
    sum(set1[, 1]) + sum(set2[, 1]) + sum(set3[, 2]),  # dS/dlambda1
    sum(set1[, 1]) + sum(set2[, 1]) + sum(set3[, 2]),  # dS/dlambda2
    sum(set1[, 1]) + sum(set2[, 1]) + sum(set3[, 2]),  # dS/dlambda3
    -sum(set3[, 2]) + sum(set3[, 1]),                  # dS/dlambda4
    -sum(set2[, 1]) + sum(set2[, 2])                   # dS/dlambda5
  )
  
  # Correct Gamma prior terms
  prior_terms <- c(
    (n2 + hyperparams$cs[1] - 1)/lambda[1] - hyperparams$ds[1],
    (n3 + hyperparams$cs[2] - 1)/lambda[2] - hyperparams$ds[2],
    (n1 + hyperparams$cs[3] - 1)/lambda[3] - hyperparams$ds[3],
    (n3 + hyperparams$cs[4] - 1)/lambda[4] - hyperparams$ds[4],
    (n2 + hyperparams$cs[5] - 1)/lambda[5] - hyperparams$ds[5]
  )
  
  # Full gradient
  grad <- prior_terms - (hyperparams$a + n1 + 2*n2 + 2*n3) * dS_dlambda / (hyperparams$b + S)
  return(grad)
}

#-------------------------------------------------------------------------------
# 4. MALA MCMC Sampler
#-------------------------------------------------------------------------------
mala_sample_lambda <- function(data, hyperparams, n_iter = 20000, burn_in = 5000, step_size = 0.01) {
  # Initialize at prior mean
  lambda_current <- hyperparams$cs / hyperparams$ds
  samples <- matrix(NA, nrow = n_iter, ncol = 5)
  accept <- 0
  
  for (iter in 1:n_iter) {
    # Compute gradient at current point
    grad_current <- compute_gradient(lambda_current, data, hyperparams)
    
    # Langevin proposal (on log scale for positivity)
    log_lambda_proposed <- log(lambda_current) + 0.5 * step_size^2 * grad_current * lambda_current + 
      rnorm(5, 0, step_size)
    lambda_proposed <- exp(log_lambda_proposed)
    
    # Compute gradient at proposed point
    grad_proposed <- compute_gradient(lambda_proposed, data, hyperparams)
    
    # Log acceptance ratio (with Jacobian)
    log_alpha <- log_posterior(lambda_proposed, data, hyperparams) - 
      log_posterior(lambda_current, data, hyperparams) +
      sum(dnorm(log(lambda_current), 
                mean = log(lambda_proposed) + 0.5 * step_size^2 * grad_proposed * lambda_proposed,
                sd = step_size, log = TRUE)) -
      sum(dnorm(log(lambda_proposed), 
                mean = log(lambda_current) + 0.5 * step_size^2 * grad_current * lambda_current,
                sd = step_size, log = TRUE))
    
    if (log(runif(1)) < log_alpha) {
      lambda_current <- lambda_proposed
      accept <- accept + 1
    }
    samples[iter, ] <- lambda_current
    
    # Adaptive tuning
    if (iter %% 500 == 0) {
      if (accept/500 < 0.4) step_size <- step_size * 0.9
      if (accept/500 > 0.6) step_size <- step_size * 1.1
      accept <- 0
    }
  }
  
  samples <- samples[(burn_in + 1):n_iter, ]
  return(samples)
}

#-------------------------------------------------------------------------------
# 5. Run Full Pipeline
#-------------------------------------------------------------------------------
# Simulate data
library(coda)

set.seed(123)
n_replicates <- 100
n_iter <- 100000
burn_in <- 20000
step_size <- 0.05
lambdas_true <- c(2.0, 1.5, 0.5, 1.2, 1.8)
alpha_true <- 1.0
avg_ci_lower <- matrix(0, nrow = n_replicates, ncol = 5)
avg_ci_upper <- matrix(0, nrow = n_replicates, ncol = 5)
coverage_counts <- rep(0, 5)

for (rep in 1:n_replicates) {
  
  avg_ci_lower[rep, ] <- ci_lower
  avg_ci_upper[rep, ] <- ci_upper
  dat <- draw_biv_data(n = 500, lambdas_true, alpha_true)
  set1 <- dat[dat[,1] == dat[,2], , drop = FALSE]
  set2 <- dat[dat[,1] < dat[,2], , drop = FALSE]
  set3 <- dat[dat[,1] > dat[,2], , drop = FALSE]
  if (nrow(set1) == 0) set1 <- matrix(1e-6, nrow = 1, ncol = 2)
  y_mean <- mean(dat)
  data <- list(set1 = set1 / y_mean, set2 = set2 / y_mean, set3 = set3 / y_mean)
  
  hyperparams <- list(
    a = 1, b = 1,
    cs = c(2, 2, 1, 2, 2),
    ds = c(1, 1, 1, 1, 1)
  )
  
  # MCMC
  samples <- mala_sample_lambda(data, hyperparams, n_iter = n_iter, burn_in = burn_in, step_size = step_size)
  
  # Compute 95% credible intervals
  ci_lower <- apply(samples, 2, quantile, probs = 0.025)
  ci_upper <- apply(samples, 2, quantile, probs = 0.975)
  
  # Count how many true lambdas fall in their intervals
  coverage_counts <- coverage_counts + ((lambdas_true >= ci_lower) & (lambdas_true <= ci_upper))
  
  if (rep %% 10 == 0) cat("Completed:", rep, "\n")
}
mean_ci_lower <- colMeans(avg_ci_lower)
mean_ci_upper <- colMeans(avg_ci_upper)

# Print them
cat("Average 95% Credible Intervals:\n")
for (i in 1:5) {
  cat(sprintf("lambda%d: [%.3f, %.3f]\n", i, mean_ci_lower[i], mean_ci_upper[i]))
}

# Report empirical coverage
coverage_rate <- coverage_counts / n_replicates
names(coverage_rate) <- paste0("lambda", 1:5)
print(coverage_rate)
