# Test SABER_fast vs SABERstan comparison
devtools::load_all(".")  # Load SABER_fast with devtools
#devtools::load_all("C:/R/SABER-main/SABER-main")  # Load SABERstan with devtools
library(dplyr)
library(tidyr)

# Define wavelength first
wavelength <- seq(400, 700, by = 10)

# Load same benthic classes in both packages
select_benthic_classes(c("Sand_2019", "Eelgrass_2019", "Mud_2019"))  # Will apply to both since both loaded

# Create synthetic data
set.seed(123)

oac <- tibble(
  name = c("chl", "a_g_440", "bb_p_550", "a_g_s_g", "a_g_s_d", "bb_p_gamma"),
  value = c(15, 1.2, 0.005, 0.014, 0.003, 0.5)
)
h_w <- 3

true_params <- stats::setNames(oac$value, oac$name)

# Generate synthetic IOPs
iop_synthetic <- iop_from_oac(wavelength, true_params)

# Create benthic reflectance mixture
benthic_fractions <- c(
  "r_rs_b_Sand_2019" = 0.6,
  "r_rs_b_Eelgrass_2019" = 0.3,
  "r_rs_b_Mud_2019" = 0.1
)
r_b_synthetic <- compute_r_rs_b_lmm(benthic_fractions)

# Generate synthetic Rrs using forward model
rrs_true <- forward_am03(
  wavelength = wavelength,
  iop = iop_synthetic,
  water_type = 2,
  theta_sun = 30,
  theta_view = 0,
  h_w = h_w,
  r_b = r_b_synthetic
)

# Add Gaussian noise
noise_level <- 0.0005
rrs_noisy <- rrs_true + rnorm(length(rrs_true), mean = 0, sd = noise_level)

# Prepare data frame for SABER
rrs_data <- data.frame(
  wavelength = wavelength,
  rrs_0m = rrs_noisy
)

cat("\n=== Synthetic Data Generated ===\n")
cat("True parameters:\n")
print(true_params)
cat("\nRrs range:", range(rrs_noisy), "\n")

# Test 1: SABER_fast pure_stan backend
cat("\n=== Testing SABER_fast (pure_stan backend) ===\n")

par_inversed <- c("chl", "a_g_440", "bb_p_550", "h_w",
                  "r_b_fractions[1]" ,  # Sand
                  "r_b_fractions[2]",   # Eelgrass
                  "r_b_fractions[3]"# Mud
)
par_fixed <- list(
  water_type = 2,
  theta_sun = 30,
  theta_view = 0,
  a_g_s_g = 0.014,
  a_g_s_d = 0.003,
  bb_p_gamma = 0.5
)

lower = c(0.5, 0.1, #0.005, #0.0002,
          0.002, #0.2,
           0.5, 0.01, 0.01, 0.01,
          0.0001
)

best = c(3, 0.75, #0.017, #0.0015,
        0.005, #0.46,
           4, 0.5, 0.1, 0.25,
         0.1
)

upper = c(100, 2,# 0.022, #0.0025,
          0.015,# 1,
            15, 1, 1, 1,
          10
)

system.time({
  result_saber <- inverse_mcmc_stan(
    rrs = rrs_data,
    forward_model = "am03",
    par_inversed = par_inversed,
    par_fixed = par_fixed,
    backend = "pure_stan",
    iterations = 1000,
    warmup = 500,
    chains = 4,
    adapt_delta = 0.90,
    return_fit = TRUE,
    verbose = FALSE
  )
}) -> time_saber

cat("\nSABER_fast completed in", time_saber["elapsed"], "seconds\n")

# Test 2: SABER_fast Bayesiantools backend
system.time({
  result_saber_bt <- inverse_mcmc(
    rrs = rrs_data,
    forward_model = "am03",
    par_inversed = c(par_inversed, "sd"),
    par_fixed = par_fixed,
    prior = NULL,
  lower = lower,
  best = NULL,
  upper = upper,
  iterations = 30000,
  burnin = 10000,
  sampler = "DEzs"
  )
}) -> time_saber_bt

# Test 2: SABERstan model_r_b_mix_amplitude
cat("\n=== Testing SABERstan (model_r_b_mix_amplitude) ===\n")

# Prepare SABERstan data using internal functions (now accessible via devtools::load_all)
wavelength_grid <- unique(rrs_data$wavelength)

devtools::load_all("C:/R/SABER-main/SABER-main")  # Load SABERstan with devtools
select_benthic_classes(c("Sand_2019", "Eelgrass_2019", "Mud_2019"))  # Will apply to both since both loaded


# Build base Stan data (includes all LUTs)
stan_data_base <- make_saber_stan_data_base(
  wavelength = wavelength_grid,
  water_type = 2,
  theta_sun_deg = 30,
  theta_view_deg = 0,
  shallow = 1,
  bottom_class_names = c("Sand_2019", "Eelgrass_2019", "Mud_2019"),
  K = 3,
  delta_scale = 0.5,
  basis = "cosine",
  pkgname = "SABERstan"
)

# Add observation data (use_measured_sigma = FALSE for synthetic data)
obs_inputs <- prepare_obs_inputs(
  df = rrs_data,
  stan_data_base = stan_data_base,
  use_measured_sigma = FALSE,
  sigma_fallback = noise_level
)

# Combine data (modifyList overwrites duplicates)
stan_data <- modifyList(stan_data_base, obs_inputs)

# Compile model
model_file <- system.file("stan", "model_r_b_mix_amplitude.stan", package = "SABERstan")
header_file <- system.file("stan", "rtm_stan_funcs.hpp", package = "SABERstan")

model_saberstan <- cmdstanr::cmdstan_model(
  stan_file = model_file,
  user_header = header_file,
  quiet = TRUE
)

# Run sampling
system.time({
  fit_saberstan <- model_saberstan$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    adapt_delta = 0.90,
    refresh = 0,
    show_messages = FALSE
  )
}) -> time_saberstan


# See all variable names
fit_saberstan$metadata()$stan_variables

# Or get summary of ALL variables (not just specific ones)
result_saberstan_all <- fit_saberstan$summary()
print(result_saberstan_all)

result_saberstan <- fit_saberstan$summary(variables = c("chl", "a_g_440", "bb_p_550", "h_w",
                                                        "r_b_mix[1]", "r_b_mix[2]", "r_b_mix[3]"))

cat("\nSABERstan completed in", time_saberstan["elapsed"], "seconds\n")

# Compare results
cat("\n=== COMPARISON ===\n")
cat("SABER_fast runtime:", round(time_saber["elapsed"], 1), "s\n")
cat("SABERstan runtime:", round(time_saberstan["elapsed"], 1), "s\n")
cat("Speedup:", round(time_saber["elapsed"] / time_saberstan["elapsed"], 2), "x\n")

cat("\nTrue values:\n")
print(true_params)

cat("\nSABER_fast estimates:\n")
print(result_saber$estimates[c("chl", "a_g_440", "bb_p_550", "h_w")])

cat("\nSABERstan estimates:\n")
print(result_saberstan[, c("variable", "mean", "sd")])

cat("\n=== Test Complete ===\n")
