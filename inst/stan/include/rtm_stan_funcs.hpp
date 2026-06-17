#pragma once
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <ostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

// SABER STAN C++ FUNCTIONS
// External C++ functions for Stan models with automatic differentiation

// ---------- geometry helper (double only; angles are data) ----------
inline void snell_law_double(double theta_view_deg, double theta_sun_deg,
                             double* view_w_rad, double* sun_w_rad) {
    constexpr double deg2rad = 3.14159265358979323846 / 180.0;
    constexpr double n_air = 1.0;
    constexpr double n_water = 1.34;

    const double tv = theta_view_deg * deg2rad;
    const double ts = theta_sun_deg  * deg2rad;

    const double sin_tv_w = (n_air / n_water) * std::sin(tv);
    const double sin_ts_w = (n_air / n_water) * std::sin(ts);

    const double stv = std::max(-1.0, std::min(1.0, sin_tv_w));
    const double sts = std::max(-1.0, std::min(1.0, sin_ts_w));

    *view_w_rad = std::asin(stv);
    *sun_w_rad  = std::asin(sts);
}


// ============================================================
// EXPORTED STAN FUNCTIONS (global namespace)
// Use Eigen::Ref so Eigen::Map from Stan binds correctly
// ============================================================

// Returns [n_wl, 2] where col1=a, col2=bb
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 2>
iop_from_oac_all(const Eigen::Ref<const Eigen::VectorXd>& wavelength,
                 const Eigen::Ref<const Eigen::VectorXd>& a_w,
                 const Eigen::Ref<const Eigen::VectorXd>& a0,
                 const Eigen::Ref<const Eigen::VectorXd>& a1,
                 const Eigen::Ref<const Eigen::VectorXd>& bb_w,
                 const T& chl,
                 const T& a_g_440,
                 const T& a_nap_440,
                 const T& a_g_s,
                 const T& a_nap_s,
                 const T& bb_p_550,
                 const T& bb_p_gamma,
                 std::ostream* pstream__) {
    using stan::math::exp;
    using stan::math::log;
    using stan::math::pow;

    const int n = wavelength.size();
    if (a_w.size() != n || a0.size() != n || a1.size() != n || bb_w.size() != n)
        throw std::domain_error("iop_from_oac_all: LUT size mismatch vs wavelength");

    Eigen::Matrix<T, Eigen::Dynamic, 2> out(n, 2);

    const T aph_440 = T(0.06) * pow(chl, 0.65);

    for (int i = 0; i < n; ++i) {
        const double wl_d = wavelength(i);
        const T wl = T(wl_d);

        // phytoplankton absorption
        T a_phy = (T(a0(i)) + T(a1(i)) * log(aph_440)) * aph_440;
        if (stan::math::value_of(a_phy) < 0.0) a_phy = T(0.0);

        const T a_g   = a_g_440   * exp(-a_g_s   * (wl - 440.0));
        const T a_nap = a_nap_440 * exp(-a_nap_s * (wl - 440.0));
        const T bb_p  = bb_p_550  * pow(wl / 550.0, -bb_p_gamma);

        out(i, 0) = T(a_w(i))  + a_phy + a_g + a_nap;
        out(i, 1) = T(bb_w(i)) + bb_p;
    }

    return out;
}


// Forward model (returns vector[n_wl])
template <typename DerivedA, typename DerivedBB, typename DerivedRB>
inline Eigen::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, 1>
forward_am03_ad(const Eigen::Ref<const Eigen::VectorXd>& wavelength,  // data
                const Eigen::MatrixBase<DerivedA>& a,
                const Eigen::MatrixBase<DerivedBB>& bb,
                int water_type,
                double theta_sun_deg,
                double theta_view_deg,
                int shallow,
                const typename DerivedA::Scalar& h_w,
                const Eigen::MatrixBase<DerivedRB>& r_b,
                std::ostream* pstream__) {
    using T = typename DerivedA::Scalar;
    using stan::math::exp;
    using stan::math::pow;

    const int n = wavelength.size();
    if (a.size() != n || bb.size() != n)
        throw std::domain_error("forward_am03_ad: size mismatch (a/bb vs wavelength)");
    if (shallow && r_b.size() != n)
        throw std::domain_error("forward_am03_ad: r_b size mismatch for shallow water");

    double view_w_rad = 0.0, sun_w_rad = 0.0;
    snell_law_double(theta_view_deg, theta_sun_deg, &view_w_rad, &sun_w_rad);
    const double cos_sun  = std::cos(sun_w_rad);
    const double cos_view = std::cos(view_w_rad);
    if (cos_sun <= 0.0 || cos_view <= 0.0)
        throw std::domain_error("forward_am03_ad: invalid geometry (cos <= 0)");

    Eigen::Matrix<T, Eigen::Dynamic, 1> rrs(n);

    for (int i = 0; i < n; ++i) {
        const T ext = a(i) + bb(i);
        if (stan::math::value_of(ext) <= 0.0) { rrs(i) = T(0.0); continue; }

        const T omega_b = bb(i) / ext;

        T f_rs;
        if (water_type == 1) {
            f_rs = T(0.095);
        } else if (water_type == 2) {
            f_rs = T(0.0512)
                   * (T(1.0) + T(4.6659) * omega_b
                      - T(7.8387) * omega_b * omega_b
                      + T(5.4571) * omega_b * omega_b * omega_b)
                   * (T(1.0) + T(0.1098) / cos_sun)
                   * (T(1.0) + T(0.4021) / cos_view);
        } else {
            throw std::domain_error("forward_am03_ad: water_type must be 1 or 2");
        }

        const T rrs_deep = f_rs * omega_b;

        if (shallow) {
            const double k0 = (water_type == 1) ? 1.0395 : 1.0546;

            const T Kd  = T(k0) * (ext / cos_sun);
            const T kuW = (ext / cos_view)
                          * pow(T(1.0) + omega_b, 3.5421)
                          * (T(1.0) - T(0.2786) / cos_sun);

            const T kuB = (ext / cos_view)
                          * pow(T(1.0) + omega_b, 2.2658)
                          * (T(1.0) - T(0.0577) / cos_sun);

            const T Ars1 = T(1.1576);
            const T Ars2 = T(1.0389);

            rrs(i) = rrs_deep * (T(1.0) - Ars1 * exp(-h_w * (Kd + kuW)))
                     + Ars2 * r_b(i) * exp(-h_w * (Kd + kuB));
        } else {
            rrs(i) = rrs_deep;
        }

        // Guard against NAN
        T val = rrs(i);
        if (!std::isfinite(stan::math::value_of(val))) val = T(0.0);
        rrs(i) = val;
    }

    return rrs;
}
// ============================================================
// STAN WRAPPER FUNCTIONS
// Explicit overloads that Stan's stanc compiler can find
// ============================================================

// Wrapper for iop_from_oac_all with var type (for autodiff)
inline Eigen::Matrix<stan::math::var, Eigen::Dynamic, 2>
iop_from_oac_all(const Eigen::Matrix<double, Eigen::Dynamic, 1>& wavelength,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1>& a_w,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1>& a0,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1>& a1,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1>& bb_w,
                 const stan::math::var& chl,
                 const stan::math::var& a_g_440,
                 const stan::math::var& a_nap_440,
                 const stan::math::var& a_g_s,
                 const stan::math::var& a_nap_s,
                 const stan::math::var& bb_p_550,
                 const stan::math::var& bb_p_gamma,
                 std::ostream* pstream__) {
    return iop_from_oac_all<stan::math::var>(wavelength, a_w, a0, a1, bb_w,
                                              chl, a_g_440, a_nap_440,
                                              a_g_s, a_nap_s, bb_p_550, bb_p_gamma,
                                              pstream__);
}

// Wrapper for forward_am03_ad with var type (for autodiff)
inline Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>
forward_am03_ad(const Eigen::Matrix<double, Eigen::Dynamic, 1>& wavelength,
                const Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>& a,
                const Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>& bb,
                int water_type,
                double theta_sun_deg,
                double theta_view_deg,
                int shallow,
                const stan::math::var& h_w,
                const Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>& r_b,
                std::ostream* pstream__) {
    return forward_am03_ad<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>,
                           Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>,
                           Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>(
        wavelength, a, bb, water_type, theta_sun_deg, theta_view_deg,
        shallow, h_w, r_b, pstream__);
}