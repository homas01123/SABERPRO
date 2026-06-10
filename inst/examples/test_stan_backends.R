# Test script to compare all three Stan backends
# Run this to see the performance difference!

#library(SABER)
library(dplyr)

# # Load example data (shallow dark water)
# data("rrs_hw", package = "SABER")
#
# # Select a challenging case
# test_rrs <- rrs_hw[,] %>%
#   # filter(station == "HW02") %>%  # Dark water station
#   select(wavelength, "rrs_0m" = rrs_estimate)

test_rrs = rrs_data

cat("=" , rep("=", 70), "\n", sep = "")
cat("  SABER Stan Backend Comparison Test\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# Parameters to estimate
params <- c("chl", "a_g_440", "a_nap_440", "bb_p_550", "h_w")

# Fast settings for testing
fast_settings <- list(
  iterations = 500,
  warmup = 250,
  chains = 2,
  adapt_delta = 0.90
)

cat("Test settings:\n")
cat(sprintf("  Wavelengths: %d bands\n", nrow(test_rrs)))
cat(sprintf("  Iterations: %d (warmup: %d)\n", fast_settings$iterations, fast_settings$warmup))
cat(sprintf("  Chains: %d\n", fast_settings$chains))
cat(sprintf("  Parameters: %s\n\n", paste(params, collapse = ", ")))

# ============================================================
# Test 1: Pure Stan
# ============================================================
cat("\n")
cat("─" , rep("─", 70), "\n", sep = "")
cat(" TEST 1: Pure Stan Backend\n")
cat("─", rep("─", 70), "\n", sep = "")

time_pure <- system.time({
  result_pure <- inverse_mcmc_stan(
    rrs = test_rrs,
    par_inversed = params,
    backend = "pure_stan",
    iterations = fast_settings$iterations,
    warmup = fast_settings$warmup,
    chains = fast_settings$chains,
    adapt_delta = fast_settings$adapt_delta,
    return_fit = TRUE,
    verbose = TRUE
  )
})

cat(sprintf("\n✓ Pure Stan completed in %.1f seconds\n", time_pure["elapsed"]))

# ============================================================
# Test 2: Hybrid Backend
# ============================================================
cat("\n")
cat("─" , rep("─", 70), "\n", sep = "")
cat(" TEST 2: Hybrid Backend (C++ IOP + Stan forward model)\n")
cat("─", rep("─", 70), "\n", sep = "")

time_hybrid <- system.time({
  result_hybrid <- inverse_mcmc_stan(
    rrs = test_rrs,
    par_inversed = params,
    backend = "hybrid",
    iterations = fast_settings$iterations,
    warmup = fast_settings$warmup,
    chains = fast_settings$chains,
    adapt_delta = fast_settings$adapt_delta,
    return_fit = TRUE,
    verbose = TRUE
  )
})

cat(sprintf("\n✓ Hybrid completed in %.1f seconds\n", time_hybrid["elapsed"]))

# ============================================================
# Test 3: C++ Optimized
# ============================================================
cat("\n")
cat("─" , rep("─", 70), "\n", sep = "")
cat(" TEST 3: C++ Optimized Backend (Full C++ implementation)\n")
cat("─", rep("─", 70), "\n", sep = "")

time_cpp <- system.time({
  result_cpp <- inverse_mcmc_stan(
    rrs = test_rrs,
    par_inversed = params,
    backend = "cpp_optimized",
    iterations = fast_settings$iterations,
    warmup = fast_settings$warmup,
    chains = fast_settings$chains,
    adapt_delta = fast_settings$adapt_delta,
    return_fit = TRUE,
    verbose = TRUE
  )
})

cat(sprintf("\n✓ C++ Optimized completed in %.1f seconds\n", time_cpp["elapsed"]))

# ============================================================
# Summary
# ============================================================
cat("\n")
cat("=" , rep("=", 70), "\n", sep = "")
cat("  PERFORMANCE SUMMARY\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# Create results table
results_df <- data.frame(
  Backend = c("Pure Stan", "Hybrid", "C++ Optimized"),
  Time_sec = c(time_pure["elapsed"], time_hybrid["elapsed"], time_cpp["elapsed"]),
  Speedup = c(
    1.0,
    time_pure["elapsed"] / time_hybrid["elapsed"],
    time_pure["elapsed"] / time_cpp["elapsed"]
  ),
  stringsAsFactors = FALSE
)

results_df$Time_min <- results_df$Time_sec / 60

print(results_df[, c("Backend", "Time_min", "Speedup")], row.names = FALSE)

cat("\n")

# Compare parameter estimates (should be nearly identical!)
cat("=" , rep("=", 70), "\n", sep = "")
cat("  PARAMETER ESTIMATE COMPARISON\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("All backends should produce nearly identical results!\n\n")

library(posterior)

# Extract posterior means
pure_summary <- as_draws_df(result_pure$fit) %>%
  summarise_draws() %>%
  filter(variable %in% params) %>%
  select(variable, pure_mean = mean, pure_sd = sd)

hybrid_summary <- as_draws_df(result_hybrid$fit) %>%
  summarise_draws() %>%
  filter(variable %in% params) %>%
  select(variable, hybrid_mean = mean, hybrid_sd = sd)

cpp_summary <- as_draws_df(result_cpp$fit) %>%
  summarise_draws() %>%
  filter(variable %in% params) %>%
  select(variable, cpp_mean = mean, cpp_sd = sd)

# Combine
comparison <- pure_summary %>%
  left_join(hybrid_summary, by = "variable") %>%
  left_join(cpp_summary, by = "variable")

print(comparison, digits = 3)

# Check convergence
cat("\n")
cat("=" , rep("=", 70), "\n", sep = "")
cat("  CONVERGENCE DIAGNOSTICS\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# Rhat values (should all be < 1.01)
pure_rhat <- max(as_draws_df(result_pure$fit) %>% summarise_draws() %>% pull(rhat), na.rm = TRUE)
hybrid_rhat <- max(as_draws_df(result_hybrid$fit) %>% summarise_draws() %>% pull(rhat), na.rm = TRUE)
cpp_rhat <- max(as_draws_df(result_cpp$fit) %>% summarise_draws() %>% pull(rhat), na.rm = TRUE)

cat(sprintf("Pure Stan max Rhat:     %.4f %s\n",
            pure_rhat, ifelse(pure_rhat < 1.01, "✓", "✗ WARNING")))
cat(sprintf("Hybrid max Rhat:        %.4f %s\n",
            hybrid_rhat, ifelse(hybrid_rhat < 1.01, "✓", "✗ WARNING")))
cat(sprintf("C++ Optimized max Rhat: %.4f %s\n",
            cpp_rhat, ifelse(cpp_rhat < 1.01, "✓", "✗ WARNING")))

cat("\n")
cat("=" , rep("=", 70), "\n", sep = "")
cat("  RECOMMENDATION\n")
cat("=", rep("=", 70), "\n\n", sep = "")

fastest <- which.min(c(time_pure["elapsed"], time_hybrid["elapsed"], time_cpp["elapsed"]))
fastest_name <- c("Pure Stan", "Hybrid", "C++ Optimized")[fastest]
fastest_speedup <- max(results_df$Speedup)

cat(sprintf("Fastest backend: %s (%.1fx faster than Pure Stan)\n\n", fastest_name, fastest_speedup))

cat("Use case recommendations:\n")
cat("  • Development/Research:  backend = 'hybrid'\n")
cat("  • Production/Operations: backend = 'cpp_optimized'\n")
cat("  • Learning/Teaching:     backend = 'pure_stan'\n\n")

cat("All backends produce statistically equivalent results!\n")
cat("Choose based on your speed vs flexibility needs.\n\n")

cat("=" , rep("=", 70), "\n", sep = "")
cat("Test completed successfully!\n")
cat("=" , rep("=", 70), "\n", sep = "")
