#' pure_water_iop
#' compute pure water iop at requested wavelength
#'
#' @param wavelength vector of requested wavelength
#' @param water_type vector of requested wavelength
#'
#' @import memoise
#'
#' @export
pure_water_iop <- function(wavelength) {
  .Call("c_pure_water_iop", as.numeric(wavelength))
}

# pure_water_iop <- memoise(function(
#     wavelength,
#     water_type = 2) {
#   # absorption (1/m)
#   a_w <- approx(a_w$wavelength, a_w$a_w, wavelength)$y # abs. of pure water [1/m]
#
#   ## backscattering (1/m)
#   if (water_type == 1) {
#     b1 <- 0.00144 #  [1/m]
#   } else if (water_type == 2) {
#     b1 <- 0.00111 #  [1/m]
#   } else {
#     rlang::abort("Water type are limited to case 1 and 2")
#   }
#   lambda1 <- 500 # [nm]
#   bb_w <- b1 * (wavelength / lambda1)^(-4.32) # [1/m]
#
#   return(
#     list(
#       "a_w" = a_w,
#       "bb_w" = bb_w
#     )
#   )
# })
