#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>

// extern declarations
extern SEXP get_cached_a_w(void);
extern SEXP get_cached_bb_w(void);
extern int wavelength_match(SEXP wavelengths);
extern SEXP c_build_cache(SEXP wavelengths);

// -- Main pure_water_iop function --
// [[register]]
SEXP c_pure_water_iop(SEXP wavelength_sexp) {
  int nprotect = 0; // Track PROTECTs

  if (!Rf_isReal(wavelength_sexp))
    Rf_error("wavelength must be a numeric (real) vector");

  // Check if cache matches, otherwise rebuild
  if (!wavelength_match(wavelength_sexp)) {
    c_build_cache(wavelength_sexp);
  }

  int n = LENGTH(wavelength_sexp);

  double* aw_ptr = REAL(get_cached_a_w());
  double* bb_w_ptr = REAL(get_cached_bb_w());

  SEXP a_w_out = PROTECT(Rf_allocVector(REALSXP, n)); nprotect++;
  SEXP bb_w_out = PROTECT(Rf_allocVector(REALSXP, n)); nprotect++;

  double* a_out = REAL(a_w_out);
  double* bb_out = REAL(bb_w_out);

  for (int i = 0; i < n; i++) {
    a_out[i] = aw_ptr[i];
    bb_out[i] = bb_w_ptr[i];
  }

  // Create list output
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2)); nprotect++;
  SET_VECTOR_ELT(out, 0, a_w_out);
  SET_VECTOR_ELT(out, 1, bb_w_out);

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2)); nprotect++;
  SET_STRING_ELT(names, 0, Rf_mkChar("a_w"));
  SET_STRING_ELT(names, 1, Rf_mkChar("bb_w"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return out;
}
