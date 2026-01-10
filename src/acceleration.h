/**
 * \file acceleration.h
 * \brief Header files for acceleration and acceleration-related
 *        functions for gravity-simulator
 *
 * \author Ching Yin Ng
 */

#ifndef ACCELERATION_H
#define ACCELERATION_H

#include "common.h"
#include "system.h"

#define ACCELERATION_METHOD_PAIRWISE 1
#define ACCELERATION_METHOD_MASSLESS 2
#define ACCELERATION_METHOD_BARNES_HUT 3

typedef struct AccelerationParam
{
    int method;
    double opening_angle;
    double softening_length;
    int max_num_particles_per_leaf;
} AccelerationParam;

/**
 * \brief Get a new acceleration parameter struct
 *
 * \return AccelerationParam
 */
AccelerationParam get_new_acceleration_param(void);

/**
 * \brief Finalize the acceleration parameter
 *
 * \param acceleration_param Pointer to the acceleration parameters
 */
void finalize_acceleration_param(AccelerationParam *restrict acceleration_param);

/**
 * \brief Wrapper function for computing acceleration
 *
 * \param[out] a Array of acceleration vectors to be modified
 * \param[in] system Pointer to the gravitational system
 * \param[in] acceleration_param Pointer to the acceleration parameters
 */
void acceleration(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
);

/**
 * \brief Compute acceleration with Barnes-Hut algorithm
 *
 * \param[out] a Array of acceleration vectors to be modified
 * \param[in] system Pointer to the gravitational system
 * \param[in] acceleration_param Pointer to the acceleration parameters
 */
void acceleration_barnes_hut(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
);

#ifdef USE_FFTW3
void acceleration_PM(
    double *restrict a,
    const CosmologicalSystem *restrict system,
    const double G,
    const int pm_grid_size
);
#endif

/**
 * \brief Benchmark acceleration
 *
 * \param system Pointer to the gravitational system
 * \param acceleration_params Array of acceleration parameters
 * \param num_acceleration_params Number of acceleration parameters
 * \param num_times_acceleration_param Array of number of times to run for each
 * acceleration parameter
 */
void benchmark_acceleration(
    const System *restrict system,
    const AccelerationParam *acceleration_params,
    const int num_acceleration_params,
    const int *restrict num_times_acceleration_param
);

#endif
