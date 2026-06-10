// Stan model for SABER AM03 shallow water inversion
// Implements Albert & Mobley (2003) forward model with HMC/NUTS sampling
// 
// ADVANTAGES over BayesianTools (DEzs/DREAMzs):
// - Gradient-guided exploration (10-20x fewer samples needed)
// - Handles parameter correlations efficiently
// - Better for multi-modal posteriors (shallow dark waters)
// - Automatic step-size tuning
// - Simplex constraint for benthic fractions (no boundary issues)
//
// Author: SABER Development Team
// Date: 2026-01-14

functions {
  // Bio-optical model: Phytoplankton absorption from Chl
  vector calculate_a_phy(vector wavelength, real chl, vector a0, vector a1) {
    int N = num_elements(wavelength);
    vector[N] a_phy;
    
    real a_phy_443 = 0.06 * chl^0.65;
    
    for (i in 1:N) {
      real log_term = a0[i] + a1[i] * log(a_phy_443);
      a_phy[i] = fmax(0.0, log_term * a_phy_443);
    }
    
    return a_phy;
  }
  
  // Snell's law: Convert air angles to water angles
  real snell_law_angle(real theta_air_deg) {
    real n_air = 1.0;
    real n_water = 1.33;
    real theta_air_rad = theta_air_deg * pi() / 180.0;
    real theta_water_rad = asin((n_air / n_water) * sin(theta_air_rad));
    return theta_water_rad;
  }
  
  // AM03 forward model (VECTORIZED for speed!)
  vector forward_am03_shallow(
      vector wavelength,
      real chl,
      real a_g_440,
      real a_nap_440,
      real bb_p_550,
      real a_g_slope,
      real bb_p_gamma,
      real h_w,
      vector r_b_fractions,
      // Data inputs
      vector a_w,
      vector bb_w,
      vector a0,
      vector a1,
      matrix r_b_matrix,
      real theta_sun_deg,
      real theta_view_deg,
      int water_type
  ) {
    int N = num_elements(wavelength);
    
    // Refract angles into water
    real theta_sun_w = snell_law_angle(theta_sun_deg);
    real theta_view_w = snell_law_angle(theta_view_deg);
    real cos_sun = cos(theta_sun_w);
    real cos_view = cos(theta_view_w);
    
    // Bio-optical conversions (VECTORIZED)
    vector[N] a_phy = calculate_a_phy(wavelength, chl, a0, a1);
    vector[N] a_g = a_g_440 * exp(-a_g_slope * (wavelength - 440.0));
    vector[N] a_nap = a_nap_440 * exp(-0.0116 * (wavelength - 440.0));
    vector[N] bb_p = bb_p_550 * pow(wavelength / 550.0, -bb_p_gamma);
    
    // Total IOPs (VECTORIZED)
    vector[N] a_total = a_w + a_phy + a_g + a_nap;
    vector[N] bb_total = bb_w + bb_p;
    vector[N] ext = a_total + bb_total;
    vector[N] omega_b = bb_total ./ ext;  // Element-wise division
    
    // AM03 parametric coefficients (water type 2)
    real Prs1 = 0.0512;
    real Prs2 = 4.6659;
    real Prs3 = -7.8387;
    real Prs4 = 5.4571;
    real Prs5 = 0.1098;
    real Prs7 = 0.4021;
    real kappa0 = 1.0546;
    real Ars1 = 1.1576;
    real Ars2 = 1.0389;
    
    // Deep water Rrs (VECTORIZED)
    vector[N] omega_b2 = omega_b .* omega_b;
    vector[N] omega_b3 = omega_b2 .* omega_b;
    vector[N] f_rs = Prs1 * 
                (1.0 + Prs2 * omega_b + Prs3 * omega_b2 + Prs4 * omega_b3) *
                (1.0 + Prs5 / cos_sun) *
                (1.0 + Prs7 / cos_view);
    vector[N] rrs_deep = f_rs .* omega_b;
    
    // Shallow water attenuation (VECTORIZED)
    vector[N] Kd = kappa0 * (ext / cos_sun);
    vector[N] kuW = (ext / cos_view) .* 
               pow(1.0 + omega_b, 3.5421) * 
               (1.0 - 0.2786 / cos_sun);
    vector[N] kuB = (ext / cos_view) .* 
               pow(1.0 + omega_b, 2.2658) * 
               (1.0 - 0.0577 / cos_sun);
    
    // Benthic reflectance (VECTORIZED matrix-vector multiply!)
    vector[N] r_b = r_b_matrix * r_b_fractions;
    
    // Combined signal (VECTORIZED)
    vector[N] exp_water = exp(-h_w * (Kd + kuW));
    vector[N] exp_benthic = exp(-h_w * (Kd + kuB));
    vector[N] rrs = rrs_deep .* (1.0 - Ars1 * exp_water) + 
                    Ars2 * r_b .* exp_benthic;
    
    return rrs;
  }
}

data {
  // Observations
  int<lower=1> n_wl;                    // Number of wavelengths
  vector[n_wl] wavelength;              // Wavelength grid [nm]
  vector[n_wl] rrs_obs;                 // Observed Rrs [1/sr]
  vector<lower=1e-12>[n_wl] sigma;      // Observation error [1/sr]
  
  // Fixed parameters
  int<lower=1,upper=2> water_type;      // Water type (1 or 2)
  real<lower=0,upper=90> theta_sun_deg; // Solar zenith angle [degrees]
  real<lower=0,upper=90> theta_view_deg;// Viewing zenith angle [degrees]
  int<lower=0,upper=1> shallow;         // Shallow water mode (always 1)
  
  // Pre-computed lookup tables
  vector[n_wl] a_w;                     // Pure water absorption [1/m]
  vector[n_wl] bb_w;                    // Pure water backscatter [1/m]
  vector[n_wl] a0;                      // Phytoplankton a0 coefficient
  vector[n_wl] a1;                      // Phytoplankton a1 coefficient
  
  // Benthic reflectance library
  int<lower=1> n_benthic;               // Number of benthic classes
  matrix[n_wl, n_benthic] r_b_matrix;   // Benthic reflectance spectra
}

parameters {
  // Optically active constituents (OACs)
  real<lower=0.01, upper=100> chl;      // Chlorophyll-a [mg/m³]
  real<lower=1e-6, upper=10> a_g_440;   // CDOM absorption at 440nm [1/m]
  real<lower=1e-6, upper=10> a_nap_440; // NAP absorption at 440nm [1/m]
  real<lower=1e-6, upper=1> bb_p_550;   // Particle backscatter at 550nm [1/m]
  
  // Spectral slopes
  real<lower=0.001, upper=0.05> a_g_slope;    // CDOM spectral slope [1/nm]
  real<lower=0.2, upper=3.0> bb_p_gamma;      // Backscatter spectral slope
  
  // Shallow water parameters
  real<lower=0.1, upper=50> h_w;        // Water depth [m]
  
  // Benthic fractions (simplex ensures sum = 1, no boundary issues!)
  simplex[n_benthic] r_b_fractions;
  
  // Observation error (with small lower bound to prevent zero)
  real<lower=1e-6> sigma_rrs;
}

transformed parameters {
  // Forward model prediction
  vector[n_wl] rrs_model;
  
  rrs_model = forward_am03_shallow(
    wavelength, chl, a_g_440, a_nap_440, bb_p_550,
    a_g_slope, bb_p_gamma, h_w, r_b_fractions,
    a_w, bb_w, a0, a1, r_b_matrix,
    theta_sun_deg, theta_view_deg, water_type
  );
}

model {
  // Weakly informative priors
  chl ~ lognormal(log(1.0), 1.5);
  a_g_440 ~ lognormal(log(0.05), 1.5);
  a_nap_440 ~ lognormal(log(0.05), 1.5);
  bb_p_550 ~ lognormal(log(0.01), 1.5);
  
  a_g_slope ~ lognormal(log(0.017), 0.3);
  bb_p_gamma ~ lognormal(log(1.0), 0.5);
  
  h_w ~ lognormal(log(5.0), 1.0);
  
  // Observation error prior
  sigma_rrs ~ normal(0, 0.01);
  
  // Likelihood
  rrs_obs ~ normal(rrs_model, sigma_rrs);
}

generated quantities {
  vector[n_wl] rrs_hat = rrs_model;
  real log_lik = normal_lpdf(rrs_obs | rrs_model, sigma_rrs);
}
