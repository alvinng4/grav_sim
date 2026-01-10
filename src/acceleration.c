/**
 * \file acceleration.c
 * \brief Functions for computing gravitational acceleration
 *
 * \author Ching-Yin Ng
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "c_traceback.h"

#include "acceleration.h"
#include "common.h"
#include "math_functions.h"
#include "system.h"
#include "utils.h"

/**
 * \brief Check the acceleration method
 *
 * \param acceleration_method Acceleration method
 */
static void check_acceleration_method(const int acceleration_method);

/**
 * \brief Compute acceleration with direct pairwise method
 *
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 */
static void acceleration_pairwise(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
);

/**
 * \brief Compute acceleration with direct pairwise method,
 *        ignoring the contribution of massless particles
 *
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 */
static void acceleration_massless(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
);

AccelerationParam get_new_acceleration_param(void)
{
    AccelerationParam acceleration_param = {
        .method = ACCELERATION_METHOD_PAIRWISE,
        .opening_angle = 1.0,
        .softening_length = 0.0,
        .max_num_particles_per_leaf = -1
    };
    return acceleration_param;
}

void finalize_acceleration_param(AccelerationParam *restrict acceleration_param)
{
    /* Check the acceleration method */
    TRY_GOTO(check_acceleration_method(acceleration_param->method), error);

    /* Check the softening length */
    if (acceleration_param->softening_length < 0.0)
    {
        THROW_FMT(
            CTB_VALUE_ERROR,
            "Softening length is negative. Got: %.3g",
            acceleration_param->softening_length
        );
        goto error;
    }

    /* Check the opening angle */
    if (acceleration_param->method == ACCELERATION_METHOD_BARNES_HUT &&
        acceleration_param->opening_angle < 0.0)
    {
        THROW_FMT(
            CTB_VALUE_ERROR,
            "Opening angle is negative. Got: %.3g",
            acceleration_param->opening_angle
        );
        goto error;
    }

    /* Check the maximum number of particles per leaf */
    if (acceleration_param->method == ACCELERATION_METHOD_BARNES_HUT)
    {
        if (acceleration_param->max_num_particles_per_leaf == -1)
        {
            acceleration_param->max_num_particles_per_leaf = 1;
        }
        else if (acceleration_param->max_num_particles_per_leaf < 1)
        {
            THROW(
                CTB_VALUE_ERROR,
                "Maximum number of particles per leaf must be positive. Got: %d",
                acceleration_param->max_num_particles_per_leaf
            );
            goto error;
        }
    }

    return;

error:
    return;
}

void acceleration(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
)
{
    switch (acceleration_param->method)
    {
        case ACCELERATION_METHOD_PAIRWISE:
            acceleration_pairwise(a, system, acceleration_param);
            break;
        case ACCELERATION_METHOD_MASSLESS:
            acceleration_massless(a, system, acceleration_param);
            break;
        case ACCELERATION_METHOD_BARNES_HUT:
            acceleration_barnes_hut(a, system, acceleration_param);
            break;
        default:
        {
            THROW_FMT(
                CTB_VALUE_ERROR,
                "Unknown acceleration method. Got: %d",
                acceleration_param->method
            );
        }
    }
}

static void check_acceleration_method(const int acceleration_method)
{
    switch (acceleration_method)
    {
        case ACCELERATION_METHOD_PAIRWISE:
        case ACCELERATION_METHOD_MASSLESS:
        case ACCELERATION_METHOD_BARNES_HUT:
            break;
        default:
        {
            THROW_FMT(
                CTB_VALUE_ERROR,
                "Unknown acceleration method. Got: %d",
                acceleration_method
            );
        }
    }

    return;
}

static void acceleration_pairwise(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
)
{
    const int num_particles = system->num_particles;
    const double *x = system->x;
    const double *m = system->m;
    const double G = system->G;
    const double softening_length = acceleration_param->softening_length;

    /* Empty the input array */
    for (int i = 0; i < num_particles; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Compute the pairwise acceleration */
    for (int i = 0; i < num_particles; i++)
    {
        const double m_i = m[i];
        for (int j = i + 1; j < num_particles; j++)
        {
            // Calculate \vec{R} and its norm
            const double R[3] = {
                x[i * 3 + 0] - x[j * 3 + 0],
                x[i * 3 + 1] - x[j * 3 + 1],
                x[i * 3 + 2] - x[j * 3 + 2]
            };
            const double R_norm = sqrt(
                R[0] * R[0] + R[1] * R[1] + R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            const double temp_value = G / (R_norm * R_norm * R_norm);
            const double m_j = m[j];
            double temp_vec[3] = {
                temp_value * R[0], temp_value * R[1], temp_value * R[2]
            };
            a[i * 3 + 0] -= temp_vec[0] * m_j;
            a[i * 3 + 1] -= temp_vec[1] * m_j;
            a[i * 3 + 2] -= temp_vec[2] * m_j;
            a[j * 3 + 0] += temp_vec[0] * m_i;
            a[j * 3 + 1] += temp_vec[1] * m_i;
            a[j * 3 + 2] += temp_vec[2] * m_i;
        }
    }

    return;
}

static void acceleration_massless(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
)
{
    const int num_particles = system->num_particles;
    const double *x = system->x;
    const double *m = system->m;
    const double G = system->G;
    const double softening_length = acceleration_param->softening_length;

    /* Empty the input array */
    for (int i = 0; i < num_particles; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Find the numbers of massive and massless particles */
    int massive_objects_count = 0;
    int massless_objects_count = 0;
    for (int i = 0; i < num_particles; i++)
    {
        if (m[i] != 0.0)
        {
            massive_objects_count++;
        }
        else
        {
            massless_objects_count++;
        }
    }

    /* Find the indices of massive and massless particles */
    int *restrict massive_indices = malloc(massive_objects_count * sizeof(int));
    int *restrict massless_indices = malloc(massless_objects_count * sizeof(int));
    massive_objects_count = 0;
    massless_objects_count = 0;

    if (massive_indices == NULL || massless_indices == NULL)
    {
        free(massive_indices);
        free(massless_indices);
        THROW(
            CTB_MEMORY_ERROR,
            "Failed to allocate memory for massive and massless indices"
        );
    }

    for (int i = 0; i < num_particles; i++)
    {
        if (m[i] != 0.0)
        {
            massive_indices[massive_objects_count] = i;
            massive_objects_count++;
        }
        else
        {
            massless_indices[massless_objects_count] = i;
            massless_objects_count++;
        }
    }

    /* Pairwise acceleration calculation for massive particles */
    for (int i = 0; i < massive_objects_count; i++)
    {
        const int idx_i = massive_indices[i];
        const double m_i = m[idx_i];
        for (int j = i + 1; j < massive_objects_count; j++)
        {
            const int idx_j = massive_indices[j];
            const double m_j = m[idx_j];
            double temp_vec[3];
            double R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            const double R_norm = sqrt(
                R[0] * R[0] + R[1] * R[1] + R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            double temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[idx_i * 3 + 0] -= temp_vec[0] * m_j;
            a[idx_i * 3 + 1] -= temp_vec[1] * m_j;
            a[idx_i * 3 + 2] -= temp_vec[2] * m_j;
            a[idx_j * 3 + 0] += temp_vec[0] * m_i;
            a[idx_j * 3 + 1] += temp_vec[1] * m_i;
            a[idx_j * 3 + 2] += temp_vec[2] * m_i;
        }
    }

    /* Acceleration calculation for massless particles due to massive particles */
    for (int i = 0; i < massive_objects_count; i++)
    {
        for (int j = 0; j < massless_objects_count; j++)
        {
            int idx_i = massive_indices[i];
            int idx_j = massless_indices[j];
            double R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            double R_norm = sqrt(
                R[0] * R[0] + R[1] * R[1] + R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            double temp_value = G / (R_norm * R_norm * R_norm);
            a[idx_j * 3 + 0] += temp_value * R[0] * m[i];
            a[idx_j * 3 + 1] += temp_value * R[1] * m[i];
            a[idx_j * 3 + 2] += temp_value * R[2] * m[i];
        }
    }

    free(massive_indices);
    free(massless_indices);

    return;
}

void benchmark_acceleration(
    const System *restrict system,
    const AccelerationParam *acceleration_params,
    const int num_acceleration_params,
    const int *restrict num_times_acceleration_param
)
{
    double *restrict reference_a = malloc(system->num_particles * 3 * sizeof(double));
    double *restrict a = malloc(system->num_particles * 3 * sizeof(double));
    if (!reference_a || !a)
    {
        THROW(CTB_MEMORY_ERROR, "Failed to allocate memory for acceleration arrays");
        goto err_malloc;
    }

    fputs("Benchmarking acceleration...\n", stdout);

    for (int i = 0; i < num_acceleration_params; i++)
    {
        const AccelerationParam *acceleration_param = &(acceleration_params[i]);
        const int num_times = num_times_acceleration_param[i];

        if (num_times <= 0)
        {
            printf("Test %d:    Skipped since num_times: %d <= 0\n\n", i, num_times);
            continue;
        }

        double *restrict run_time = calloc(num_times, sizeof(double));
        double mae = 0.0;

        if (!run_time)
        {
            THROW(CTB_MEMORY_ERROR, "Failed to allocate memory for runtime array");
            goto err_malloc;
        }

        for (int j = 0; j < num_times; j++)
        {
            if (i == 0 && j == 0)
            {
                double start_time = grav_get_current_time();
                TRY_GOTO(acceleration(reference_a, system, acceleration_param), error);
                double end_time = grav_get_current_time();
                run_time[j] += (end_time - start_time);
            }
            else
            {
                double start_time = grav_get_current_time();
                TRY_GOTO(acceleration(a, system, acceleration_param), error);
                double end_time = grav_get_current_time();
                run_time[j] += (end_time - start_time);
            }

            // Calculate the MAE
            if (i != 0 && j == 0)
            {
                for (int k = 0; k < system->num_particles; k++)
                {
                    const double diff[3] = {
                        reference_a[k * 3 + 0] - a[k * 3 + 0],
                        reference_a[k * 3 + 1] - a[k * 3 + 1],
                        reference_a[k * 3 + 2] - a[k * 3 + 2]
                    };
                    mae += fabs(diff[0]) + fabs(diff[1]) + fabs(diff[2]);
                }
                mae /= system->num_particles;
            }
        }

        printf("Test %d:", i);
        switch (acceleration_param->method)
        {
            case ACCELERATION_METHOD_PAIRWISE:
                fputs("    Method: Pairwise\n", stdout);
                break;
            case ACCELERATION_METHOD_MASSLESS:
                fputs("    Method: Massless\n", stdout);
                break;
            case ACCELERATION_METHOD_BARNES_HUT:
                fputs("    Method: Barnes-Hut\n", stdout);
                break;
            default:
                THROW_FMT(
                    CTB_VALUE_ERROR,
                    "Unknown acceleration method. Got: %d",
                    acceleration_param->method
                );
                goto err_unknown_acceleration_method;
        }

        printf("    Number of times: %d\n", num_times);
        printf(
            "    Avg time: %.3g (+- %.3g) s\n",
            compute_mean(run_time, num_times),
            compute_std(run_time, num_times, 1) / sqrt(num_times)
        );
        printf("    MAE: %.3g\n", mae);
        printf("\n");

        free(run_time);
    }

    free(reference_a);
    free(a);

    return;

err_unknown_acceleration_method:
error:
    free(run_time);
err_malloc:
    free(reference_a);
    free(a);
    return;
}
