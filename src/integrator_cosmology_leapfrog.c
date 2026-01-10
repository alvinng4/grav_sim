/**
 * \file integrator_cosmology_leapfrog.c
 * \brief Leapfrog integrator for cosmological simulations.
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
#include "cosmology.h"
#include "integrator.h"
#include "math_functions.h"
#include "output.h"
#include "progress_bar.h"
#include "settings.h"
#include "system.h"
#include "utils.h"

void leapfrog_cosmology(
    CosmologicalSystem *restrict system,
    OutputParam *restrict output_param,
    SimulationStatus *restrict simulation_status,
    Settings *restrict settings,
    const double a_final,
    const int num_steps,
    const int pm_grid_size
)
{
    /* Declare variables */
    const int num_particles = system->num_particles;
    double *restrict x = system->x;
    double *restrict v = system->v;
    const double H0 = system->h * 100.0;
    const double omega_m = system->omega_m;
    const double omega_lambda = system->omega_lambda;

    bool is_output = (output_param->method != OUTPUT_METHOD_DISABLED);
    int *restrict output_count_ptr = &(output_param->output_count_);
    const double output_interval = output_param->output_interval;
    double next_output_time = output_interval;

    const double t0 = log(system->scale_factor);
    const double tf = log(a_final);
    double *restrict t_ptr = &(simulation_status->t);
    *t_ptr = t0;
    int64 *restrict num_steps_ptr = &(simulation_status->num_steps);

    const bool enable_progress_bar = settings->enable_progress_bar;

    double H_a;

    /* Allocate memory */
    double *restrict momentum = malloc(num_particles * 3 * sizeof(double));
    double *restrict a = malloc(num_particles * 3 * sizeof(double));

    // Check if memory allocation is successful
    if (!momentum || !a)
    {
        THROW(CTB_MEMORY_ERROR, "Failed to allocate memory for arrays");
        goto err_memory;
    }

    /* Get mean background density */
    const double G =
        6.67430e-8 *
        (system->unit_mass * system->unit_time * system->unit_time /
         (system->unit_length * system->unit_length * system->unit_length));

    /* Initial output */
    if (is_output && output_param->output_initial)
    {
        TRY_GOTO(
            output_snapshot_cosmology(
                output_param, system, simulation_status, settings
            ),
            err_initial_output
        );
    }

    /* Set periodic boundary conditions */
    TRACE(set_periodic_boundary_conditions(system));

    /* Initialize momentum */
    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            momentum[i * 3 + j] =
                (system->scale_factor) * (system->scale_factor) * v[i * 3 + j];
        }
    }

    /* Compute initial acceleration */
    TRY_GOTO(acceleration_PM(a, system, G, pm_grid_size), err_acceleration);

    /* Main Loop */
    double dt = (tf - t0) / num_steps;
    int64 total_num_steps = (int64)ceil((tf - t0) / dt);
    ProgressBarParam progress_bar_param;
    if (enable_progress_bar)
    {
        TRY_GOTO(
            start_progress_bar(&progress_bar_param, total_num_steps),
            err_start_progress_bar
        );
    }

    simulation_status->dt = dt;
    *num_steps_ptr = 0;
    while (*num_steps_ptr < total_num_steps)
    {
        /* Check dt overshoot */
        if (*t_ptr + dt > tf)
        {
            dt = tf - *t_ptr;
        }
        simulation_status->dt = dt;

        /* Kick (p_1/2) */
        H_a = TRACE(compute_H(system->scale_factor, H0, omega_m, omega_lambda));
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                momentum[i * 3 + j] -= (0.5 * dt) * a[i * 3 + j] / H_a;
            }
        }
        *t_ptr += 0.5 * dt;
        system->scale_factor = exp(*t_ptr);

        /* Drift (x_1) */
        H_a = TRACE(compute_H(system->scale_factor, H0, omega_m, omega_lambda));
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] += dt * momentum[i * 3 + j] /
                                ((system->scale_factor) * (system->scale_factor) * H_a);
            }
        }

        /* Set periodic boundary conditions */
        TRACE(set_periodic_boundary_conditions(system));

        /* Kick (p_1) */
        TRY_GOTO(acceleration_PM(a, system, G, pm_grid_size), err_acceleration);

        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                momentum[i * 3 + j] -= (0.5 * dt) * a[i * 3 + j] / H_a;
            }
        }

        (*num_steps_ptr)++;
        *t_ptr = t0 + (*num_steps_ptr) * dt;
        system->scale_factor = exp(*t_ptr);

        /* Store solution */
        if (is_output && system->scale_factor >= next_output_time)
        {
            /* Get velocity from momentum */
            for (int i = 0; i < num_particles; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    v[i * 3 + j] = momentum[i * 3 + j] /
                                   (system->scale_factor * system->scale_factor);
                }
            }
            TRY_GOTO(output_snapshot_cosmology(
                output_param, system, simulation_status, settings
            ));

            next_output_time = (*output_count_ptr) * output_interval;
        }

        if (enable_progress_bar)
        {
            TRACE(update_progress_bar(&progress_bar_param, *num_steps_ptr, false));
        }

        /* Check exit */
        if (*(settings->is_exit_ptr))
        {
            break;
        }
    }

    if (enable_progress_bar)
    {
        TRACE(update_progress_bar(&progress_bar_param, *num_steps_ptr, true));
    }

    TRACE_BLOCK(free(momentum); free(a););

    return;

err_output:
err_acceleration:
err_start_progress_bar:
err_initial_output:
err_memory:
    TRACE_BLOCK(free(momentum); free(a););

    return;
}
