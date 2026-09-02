#ifndef COSMOLOGY_H
#define COSMOLOGY_H


/**
 * \brief Compute the curvature density parameter Omega_k.
 * 
 * \param omega_m Matter density parameter.
 * \param omega_lambda Dark energy density parameter.
 * 
 * \return Curvature density parameter Omega_k.
 */
double compute_omega_k(
    const double omega_m,
    const double omega_lambda
);

/**
 * \brief Compute the Hubble parameter H(a).
 * 
 * \param a Scale factor.
 * \param H0 Hubble constant.
 * \param omega_m Matter density parameter.
 * \param omega_lambda Dark energy density parameter.
 * 
 * \return Hubble parameter H(a).
 */
double compute_H(
    const double a,
    const double h0,
    const double omega_m,
    const double omega_lambda
);

#endif
