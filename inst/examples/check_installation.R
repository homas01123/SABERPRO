# Pre-Testing Checklist for SABER Stan Backends
# Run this before testing the new Stan implementations

cat("
╔════════════════════════════════════════════════════════════════╗
║          SABER Stan Backend - Pre-Testing Checklist            ║
╚════════════════════════════════════════════════════════════════╝
")

# Step 1: Load package with new changes
cat("\n[1/5] Loading SABER package with devtools...\n")
if (!requireNamespace("devtools", quietly = TRUE)) {
  cat("  ✗ devtools not found. Installing...\n")
  install.packages("devtools")
}

library(devtools)
devtools::load_all(".")
cat("  ✓ SABER package loaded\n")

# Step 2: Check if Stan files exist
cat("\n[2/5] Checking Stan model files...\n")
stan_files <- c(
  "am03_shallow.stan" = "Pure Stan backend",
  "am03_hybrid.stan" = "Hybrid backend (C++ IOP + Stan forward model)",
  "am03_cpp_optimized.stan" = "C++ Optimized backend (full C++)"
)

all_exist <- TRUE
for (file in names(stan_files)) {
  path <- system.file("stan", file, package = "SABER")
  if (file.exists(path) && path != "") {
    cat(sprintf("  ✓ %s (%s)\n", file, stan_files[file]))
  } else {
    cat(sprintf("  ✗ %s NOT FOUND\n", file))
    all_exist <- FALSE
  }
}

# Step 3: Check C++ header for external functions
cat("\n[3/5] Checking C++ header for external functions...\n")
cpp_header <- system.file("stan/include/rtm_stan_funcs.hpp", package = "SABER")
if (file.exists(cpp_header) && cpp_header != "") {
  cat(sprintf("  ✓ rtm_stan_funcs.hpp found at:\n    %s\n", cpp_header))
  
  # Check file size (should be ~180 lines, ~6KB)
  file_info <- file.info(cpp_header)
  cat(sprintf("    Size: %.1f KB\n", file_info$size / 1024))
  
  if (file_info$size < 1000) {
    cat("  ⚠ WARNING: File seems too small. May be incomplete.\n")
  }
} else {
  cat("  ✗ rtm_stan_funcs.hpp NOT FOUND\n")
  cat("    This is required for 'hybrid' and 'cpp_optimized' backends!\n")
  all_exist <- FALSE
}

# Step 4: Check Stan/cmdstanr installation
cat("\n[4/5] Checking Stan installation...\n")
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("  ✗ cmdstanr not installed\n")
  cat("    Install with: install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))\n")
} else {
  library(cmdstanr)
  
  # Check CmdStan
  cmdstan_ver <- tryCatch(
    cmdstan_version(error_on_NA = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(cmdstan_ver) || is.na(cmdstan_ver)) {
    cat("  ✗ CmdStan not found\n")
    cat("    Install with: cmdstanr::install_cmdstan()\n")
  } else {
    cat(sprintf("  ✓ cmdstanr installed\n"))
    cat(sprintf("  ✓ CmdStan version: %s\n", cmdstan_ver))
    cat(sprintf("    Path: %s\n", cmdstan_path()))
  }
}

# Step 5: Quick function availability check
cat("\n[5/5] Checking function availability...\n")
functions_to_check <- c(
  "inverse_mcmc_stan",
  "compare_mcmc_methods",
  "validate_stan_installation"
)

for (func in functions_to_check) {
  if (exists(func, mode = "function")) {
    cat(sprintf("  ✓ %s() available\n", func))
  } else {
    cat(sprintf("  ✗ %s() NOT FOUND\n", func))
    all_exist <- FALSE
  }
}

# Summary
cat("\n")
cat("═", rep("═", 70), "\n", sep = "")
if (all_exist) {
  cat("  ✓ ALL CHECKS PASSED - Ready to test!\n")
  cat("═", rep("═", 70), "\n", sep = "")
  cat("\n")
  cat("Next steps:\n")
  cat("  1. Test pure_stan:      result <- inverse_mcmc_stan(rrs, backend='pure_stan')\n")
  cat("  2. Test hybrid:         result <- inverse_mcmc_stan(rrs, backend='hybrid')\n")
  cat("  3. Test cpp_optimized:  result <- inverse_mcmc_stan(rrs, backend='cpp_optimized')\n")
  cat("  4. Run full benchmark:  source('inst/examples/test_stan_backends.R')\n")
} else {
  cat("  ✗ SOME CHECKS FAILED - See messages above\n")
  cat("═", rep("═", 70), "\n", sep = "")
  cat("\n")
  cat("Action required:\n")
  cat("  - Make sure you're in the SABER package directory\n")
  cat("  - Re-run devtools::load_all()\n")
  cat("  - Check file paths above\n")
}

cat("\n")
cat("═", rep("═", 70), "\n", sep = "")
cat("  About saber-lib-main directory:\n")
cat("═", rep("═", 70), "\n", sep = "")
cat("\n")
cat("Q: Do I still need the saber-lib-main folder?\n")
cat("A: NO - it's optional now!\n\n")
cat("  • We copied rtm_stan_funcs.hpp from saber-lib-main/src/\n")
cat("  • It's now in SABER at: inst/stan/include/rtm_stan_funcs.hpp\n")
cat("  • The header is self-contained (no linking needed)\n")
cat("  • saber-lib-main was just the SOURCE\n\n")
cat("You can:\n")
cat("  ✓ Keep it for reference/documentation\n")
cat("  ✓ Delete it if you want to save space\n")
cat("  ✓ Update it separately if colleague makes changes\n\n")
cat("The SABER package is now independent!\n")
cat("\n")
