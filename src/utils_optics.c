#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>

// Last-used cached geometry
static double cached_theta_view = -9999;
static double cached_theta_sun  = -9999;
static double cached_view_w     = 0;
static double cached_sun_w      = 0;
// static double cached_rho_L      = 0;

// [[register]]
SEXP c_snell_law(SEXP theta_view_sexp, SEXP theta_sun_sexp) {
  if (!Rf_isReal(theta_view_sexp) || !Rf_isReal(theta_sun_sexp))
    Rf_error("Both theta_view and theta_sun must be numeric scalars");

  double theta_view = REAL(theta_view_sexp)[0];
  double theta_sun  = REAL(theta_sun_sexp)[0];

  // Return cache if angles match exactly
  if (theta_view == cached_theta_view && theta_sun == cached_theta_sun) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(cached_view_w));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(cached_sun_w));
    // SET_VECTOR_ELT(out, 2, Rf_ScalarReal(cached_rho_L));

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("view_w"));
    SET_STRING_ELT(names, 1, Rf_mkChar("sun_w"));
    // SET_STRING_ELT(names, 2, Rf_mkChar("rho_L"));
    Rf_setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(2);
    return out;
  }

  // Convert from degrees to radians
  double theta_view_rad = theta_view * M_PI / 180.0;
  double theta_sun_rad  = theta_sun  * M_PI / 180.0;

  double n_air = 1.0;
  double n_w   = 1.33;

  // Snell's law to get underwater angles (in radians)
  double view_w = asin((n_air / n_w) * sin(theta_view_rad));
  double sun_w  = asin((n_air / n_w) * sin(theta_sun_rad));

  // Fresnel reflectance (scalar rho_L)
  // double num1 = pow(sin(theta_view_rad - view_w), 2);
  // double den1 = pow(sin(theta_view_rad + view_w), 2);
  // double term1 = (den1 > 0) ? (num1 / den1) : 0.0;
  //
  // double num2 = pow(tan(theta_view_rad - view_w), 2);
  // double den2 = pow(tan(theta_view_rad + view_w), 2);
  // double term2 = (den2 > 0) ? (num2 / den2) : 0.0;
  //
  // double rho_L = 0.5 * (term1 + term2);

  // Cache the result
  cached_theta_view = theta_view;
  cached_theta_sun  = theta_sun;
  cached_view_w     = view_w;
  cached_sun_w      = sun_w;
  // cached_rho_L      = rho_L;

  // Return result
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(view_w));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(sun_w));
  // SET_VECTOR_ELT(out, 2, Rf_ScalarReal(rho_L));

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, Rf_mkChar("view_w"));
  SET_STRING_ELT(names, 1, Rf_mkChar("sun_w"));
  // SET_STRING_ELT(names, 2, Rf_mkChar("rho_L"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);
  return out;
}
