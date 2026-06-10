#' Transmit rrs(0+) to rrs(0-)
#'
#' @export
rrs_0p_to_0m <- function(rrs_0p) {
  rrs_0m <- rrs_0p / (0.52 + 1.7 * rrs_0p)
  return(rrs_0m)
}

#' Transmit rrs(0-) to Rrs(0+)
#'
#' @export
rrs_0m_to_0p <- function(rrs_0m) {
  rrs_0p <- -0.52 * rrs_0m / (1.7 * rrs_0m - 1)
  return(rrs_0p)
}

#' snell law, compute refracted angles and reflection
#' @param theta_view sensor zenith angle \[deg\]
#' @param theta_sun solar zenith angle \[deg\]
#'
#' @return list with `view_w`, `sun_w`, `rho_L`
#' @export
snell_law <- function(theta_view, theta_sun) {
  stopifnot(is.numeric(theta_view), length(theta_view) == 1)
  stopifnot(is.numeric(theta_sun), length(theta_sun) == 1)
  .Call("c_snell_law", as.numeric(theta_view), as.numeric(theta_sun))
}
# snell_law <- memoise(function(theta_view, theta_sun) { # Function to convert above water to under water geometry
#
#   # Index of refrations (real)
#   n_air <- 1 # air index of refration (real part)
#   n_w <- 1.33 # water index of refration (real part)
#
#   # Angles from the water
#
#   # from deg to rad
#   theta_view <- theta_view * (180 / pi) # rad
#   theta_sun <- theta_sun * (180 / pi) # rad
#
#   # angles inside the water in rad
#   view_w <- asin((n_air / n_w) * sin(theta_view)) # rad
#   sun_w <- asin((n_air / n_w) * sin(theta_sun)) # rad
#
#   # Fresnel Law
#
#   rho_L <- (1 / 2) * abs(((sin(theta_view - view_w)^2) / (sin(theta_view + view_w)^2)) + ((tan(theta_view - view_w)^2) / (tan(theta_view + view_w)^2)))
#   return(data.frame("view_w" = view_w, "sun_w" = sun_w, "rho_L" = rho_L))
# })
