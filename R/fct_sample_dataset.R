#' pepare_samples
#'
#' Helper function: truncate and sample

prepare_samples <- function(vec, bounds = NULL, n = 100) {
  vec <- vec[!is.na(vec)]
  if (!is.null(bounds)) {
    vec <- vec[vec >= bounds[1] & vec <= bounds[2]]
  }
  sample(vec, size = n)
}

#' sample_nomad
#'
#' sample the NOMAD dataset
#'
#' @author Soham Mukherjee
#'
#' @import dplyr
#'
#' @returns input data prior distribution to be used by ?

sample_nomad <- function(
    truncate_chl = NULL,
    truncate_ag = NULL,
    truncate_ad = NULL,
    truncate_bbp = NULL,
    sample_count = 100) {
  # Load NOMAD data
  bgc_data <- nomad

  bgc_data <- bgc_data %>%
    mutate(across(where(is.numeric), ~ na_if(., -999)))

  # Prepare samples
  chl_sample <- prepare_samples(bgc_data$chl, truncate_chl, sample_count)
  ag_sample <- prepare_samples(bgc_data$ag443, truncate_ag, sample_count)
  ad_sample <- prepare_samples(bgc_data$ad443, truncate_ad, sample_count)
  bbp_sample <- prepare_samples(bgc_data$bb555, truncate_bbp, sample_count)

  return(list(
    chl = chl_sample,
    ag443 = ag_sample,
    ad443 = ad_sample,
    bbp555 = bbp_sample
  ))
}
