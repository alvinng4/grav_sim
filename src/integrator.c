/**
 * \file integrator.c
 * \brief Function definitions for integrator-related functions and simple integrators
 *
 * This file contains the definitions for integrator-related functions and
 * simple integrators including Euler, Euler-Cromer, Runge-Kutta 4th order (RK4),
 * and Leapfrog.
 *
 * \author Ching-Yin Ng
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acceleration.h"
#include "c_traceback.h"
#include "common.h"
#include "integrator.h"
#include "output.h"
#include "progress_bar.h"
#include "settings.h"
#include "system.h"

/**
 * \brief Euler first-order integrator
 *
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 */
static void euler(
    System *restrict system,
    IntegratorParam *restrict integrator_param,
    AccelerationParam *restrict acceleration_param,
    OutputParam *restrict output_param,
    SimulationStatus *restrict simulation_status,
    Settings *restrict settings,
    const double tf
);

/**
 * \brief Euler-Cromer first-order integrator
 *
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 */
static void euler_cromer(
    System *restrict system,
    IntegratorParam *restrict integrator_param,
    AccelerationParam *restrict acceleration_param,
    OutputParam *restrict output_param,
    SimulationStatus *restrict simulation_status,
    Settings *restrict settings,
    const double tf
);

/**
 * \brief Runge-Kutta 4th order (RK4) integrator
 *
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 */
static void
rk4(System *restrict system,
    IntegratorParam *restrict integrator_param,
    AccelerationParam *restrict acceleration_param,
    OutputParam *restrict output_param,
    SimulationStatus *restrict simulation_status,
    Settings *restrict settings,
    const double tf);

/**
 * \brief Leapfrog integrator
 *
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 */
static void leapfrog(
    System *restrict system,
    IntegratorParam *restrict integrator_param,
    AccelerationParam *restrict acceleration_param,
    OutputParam *restrict output_param,
    SimulationStatus *restrict simulation_status,
    Settings *restrict settings,
    const double tf
);

IntegratorParam get_new_integrator_param(void)
{
    IntegratorParam integrator_param = {
        .integrator = -1,
        .dt = -1.0,
        .tolerance = -1.0,
        .initial_dt = -1.0,
        .whfast_remove_invalid_particles = true,
    };
    return integrator_param;
}

void finalize_integration_param(IntegratorParam *restrict integration_param)
{
    if (!integration_param)
    {
        THROW(CTB_POINTER_ERROR, "integration_param is NULL");
        return;
    }

    if (integration_param->integrator != INTEGRATOR_EULER &&
        integration_param->integrator != INTEGRATOR_EULER_CROMER &&
        integration_param->integrator != INTEGRATOR_RK4 &&
        integration_param->integrator != INTEGRATOR_LEAPFROG &&
        integration_param->integrator != INTEGRATOR_RKF45 &&
        integration_param->integrator != INTEGRATOR_DOPRI &&
        integration_param->integrator != INTEGRATOR_DVERK &&
        integration_param->integrator != INTEGRATOR_RKF78 &&
        integration_param->integrator != INTEGRATOR_IAS15 &&
        integration_param->integrator != INTEGRATOR_WHFAST)
    {
        THROW_FMT(
            CTB_VALUE_ERROR,
            "Unknown integrator. Got: %d",
            integration_param->integrator
        );
        return;
    }

    if (integration_param->integrator == INTEGRATOR_EULER ||
        integration_param->integrator == INTEGRATOR_EULER_CROMER ||
        integration_param->integrator == INTEGRATOR_RK4 ||
        integration_param->integrator == INTEGRATOR_LEAPFROG ||
        integration_param->integrator == INTEGRATOR_WHFAST)
    {
        if (integration_param->dt <= 0.0)
        {
            THROW_FMT(
                CTB_VALUE_ERROR,
                "dt must be positive. Got: %g",
                integration_param->dt
            );
            return;
        }
    }

    else
    {
        if (integration_param->tolerance <= 0.0)
        {
            THROW_FMT(
                CTB_VALUE_ERROR,
                "tolerance must be positive. Got: %g",
                integration_param->tolerance
            );
            return;
        }
    }

    return;
}

void integrator_launch_simulation(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    OutputParam *output_param,
    SimulationStatus *simulation_status,
    Settings *settings,
    const double tf
)
{
    switch (integrator_param->integrator)
    {
        case INTEGRATOR_EULER:
            TRY_GOTO(euler(
                system,
                integrator_param,
                acceleration_param,
                output_param,
                simulation_status,
                settings,
                tf
            ), error);
            break;
        case INTEGRATOR_EULER_CROMER:
            TRY_GOTO(euler_cromer(
                system,
                integrator_param,
                acceleration_param,
                output_param,
                simulation_status,
                settings,
                tf
            ), error);
            break;
        case INTEGRATOR_RK4:
            TRY_GOTO(
                rk4(system,
                    integrator_param,
                    acceleration_param,
                    output_param,
                    simulation_status,
                    settings,
                    tf), error
            );
            break;
        case INTEGRATOR_LEAPFROG:
            TRY_GOTO(leapfrog(
                system,
                integrator_param,
                acceleration_param,
                output_param,
                simulation_status,
                settings,
                tf
            ), error);
            break;
        case INTEGRATOR_RKF45:
        case INTEGRATOR_DOPRI:
        case INTEGRATOR_DVERK:
        case INTEGRATOR_RKF78:
            TRY_GOTO(rk_embedded(
                system,
                integrator_param,
                acceleration_param,
                output_param,
                simulation_status,
                settings,
                tf
            ), error);
            break;
        case INTEGRATOR_IAS15:
            TRY_GOTO(ias15(
                system,
                integrator_param,
                acceleration_param,
                output_param,
                simulation_status,
                settings,
                tf
            ), error);
            break;
        case INTEGRATOR_WHFAST:
            TRY_GOTO(whfast(
                system,
                integrator_param,
                acceleration_param,
                output_param,
                simulation_status,
                settings,
                tf
            ), error);
            break;
        default:
            THROW(CTB_VALUE_ERROR, "Invalid integrator");
            goto error;
    }

    return;

error:
    return;
}

static void euler(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    OutputParam *output_param,
    SimulationStatus *simulation_status,
    Settings *settings,
    const double tf
)
{
    /* Declare variables */

    const int num_particles = system->num_particles;
    double *restrict x = system->x;
    double *restrict v = system->v;

    double dt = integrator_param->dt;

    bool is_output = (output_param->method != OUTPUT_METHOD_DISABLED);
    int *restrict output_count_ptr = &(output_param->output_count_);
    const double output_interval = output_param->output_interval;
    double next_output_time = output_interval;

    double *restrict t_ptr = &(simulation_status->t);
    int64 *restrict num_steps_ptr = &(simulation_status->num_steps);

    const bool enable_progress_bar = settings->enable_progress_bar;

    /* Allocate memory */
    double *restrict x_0 = malloc(num_particles * 3 * sizeof(double));
    double *restrict v_0 = malloc(num_particles * 3 * sizeof(double));
    double *restrict a = malloc(num_particles * 3 * sizeof(double));

    // Compensated summation
    double *restrict x_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict v_err_comp_sum = calloc(num_particles * 3, sizeof(double));

    // Check if memory allocation is successful
    if (!x_0 || !v_0 || !a || !x_err_comp_sum || !v_err_comp_sum)
    {
        THROW(CTB_MEMORY_ERROR, "Failed to allocate memory for arrays");
        goto err_memory;
    }

    /* Initial output */
    if (is_output && output_param->output_initial)
    {
        TRY_GOTO(output_snapshot(
            output_param,
            system,
            integrator_param,
            acceleration_param,
            simulation_status,
            settings
        ), err_initial_output);
        
    }

    /* Main Loop */
    int64 total_num_steps = (int64)ceil(tf / dt);
    ProgressBarParam progress_bar_param;
    if (enable_progress_bar)
    {
        TRY_GOTO(start_progress_bar(&progress_bar_param, total_num_steps), err_start_progress_bar);
    }

    *t_ptr = 0.0;
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

        memcpy(x_0, x, num_particles * 3 * sizeof(double));
        memcpy(v_0, v, num_particles * 3 * sizeof(double));

        /* Compute acceleration */
        TRY_GOTO(acceleration(a, system, acceleration_param), err_acceleration);

        /* Update step */
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;

                x[i * 3 + j] = x_0[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                v[i * 3 + j] = v_0[i * 3 + j] + v_err_comp_sum[i * 3 + j];

                x_err_comp_sum[i * 3 + j] += x_0[i * 3 + j] - x[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += v_0[i * 3 + j] - v[i * 3 + j];
            }
        }
        (*num_steps_ptr)++;
        *t_ptr = (*num_steps_ptr) * dt;

        /* Store solution */
        if (is_output && *t_ptr >= next_output_time)
        {
            TRY_GOTO(output_snapshot(
                output_param,
                system,
                integrator_param,
                acceleration_param,
                simulation_status,
                settings
            ), err_output);

            next_output_time = (*output_count_ptr) * output_interval;
        }

        if (enable_progress_bar)
        {
            update_progress_bar(&progress_bar_param, *num_steps_ptr, false);
        }

        /* Check exit */
        if (*(settings->is_exit_ptr))
        {
            break;
        }
    }

    if (enable_progress_bar)
    {
        update_progress_bar(&progress_bar_param, *num_steps_ptr, true);
    }

    free(x_0);
    free(v_0);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    return;

err_output:
err_acceleration:
err_start_progress_bar:
err_initial_output:
err_memory:
    free(x_0);
    free(v_0);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    return;
}

static void euler_cromer(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    OutputParam *output_param,
    SimulationStatus *simulation_status,
    Settings *settings,
    const double tf
)
{
    /* Declare variables */
    

    const int num_particles = system->num_particles;
    double *restrict x = system->x;
    double *restrict v = system->v;

    double dt = integrator_param->dt;

    bool is_output = (output_param->method != OUTPUT_METHOD_DISABLED);
    int *restrict output_count_ptr = &(output_param->output_count_);
    const double output_interval = output_param->output_interval;
    double next_output_time = output_interval;

    double *restrict t_ptr = &(simulation_status->t);
    int64 *restrict num_steps_ptr = &(simulation_status->num_steps);

    const bool enable_progress_bar = settings->enable_progress_bar;

    /* Allocate memory */
    double *restrict x_0 = malloc(num_particles * 3 * sizeof(double));
    double *restrict v_0 = malloc(num_particles * 3 * sizeof(double));
    double *restrict a = malloc(num_particles * 3 * sizeof(double));

    // Compensated summation
    double *restrict x_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict v_err_comp_sum = calloc(num_particles * 3, sizeof(double));

    // Check if memory allocation is successful
    if (!x_0 || !v_0 || !a || !x_err_comp_sum || !v_err_comp_sum)
    {
        THROW(CTB_MEMORY_ERROR, "Failed to allocate memory for arrays");
        goto err_memory;
    }

    /* Initial output */
    if (is_output && output_param->output_initial)
    {
        TRY_GOTO(output_snapshot(
            output_param,
            system,
            integrator_param,
            acceleration_param,
            simulation_status,
            settings
        ), err_initial_output);
        
    }

    /* Main Loop */
    int64 total_num_steps = (int64)ceil(tf / dt);
    ProgressBarParam progress_bar_param;
    if (enable_progress_bar)
    {
        TRY_GOTO(start_progress_bar(&progress_bar_param, total_num_steps), err_start_progress_bar);
    }

    *t_ptr = 0.0;
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

        memcpy(x_0, x, num_particles * 3 * sizeof(double));
        memcpy(v_0, v, num_particles * 3 * sizeof(double));

        /* Compute acceleration */
        TRY_GOTO(acceleration(a, system, acceleration_param), err_acceleration);

        /* Update step */
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;
                v[i * 3 + j] = v_0[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += v_0[i * 3 + j] - v[i * 3 + j];

                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                x[i * 3 + j] = x_0[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += x_0[i * 3 + j] - x[i * 3 + j];
            }
        }
        (*num_steps_ptr)++;
        *t_ptr = (*num_steps_ptr) * dt;

        /* Store solution */
        if (is_output && *t_ptr >= next_output_time)
        {
            TRY_GOTO(output_snapshot(
                output_param,
                system,
                integrator_param,
                acceleration_param,
                simulation_status,
                settings
            ), err_output);
            

            next_output_time = (*output_count_ptr) * output_interval;
        }

        if (enable_progress_bar)
        {
            update_progress_bar(&progress_bar_param, *num_steps_ptr, false);
        }

        /* Check exit */
        if (*(settings->is_exit_ptr))
        {
            break;
        }
    }

    if (enable_progress_bar)
    {
        update_progress_bar(&progress_bar_param, *num_steps_ptr, true);
    }

    free(x_0);
    free(v_0);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    return;

err_output:
err_acceleration:
err_start_progress_bar:
err_initial_output:
err_memory:
    free(x_0);
    free(v_0);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    return;
}

static void
rk4(System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    OutputParam *output_param,
    SimulationStatus *simulation_status,
    Settings *settings,
    const double tf)
{
    /* Declare variables */
    

    const int num_particles = system->num_particles;
    double *restrict x = system->x;
    double *restrict v = system->v;

    double dt = integrator_param->dt;

    bool is_output = (output_param->method != OUTPUT_METHOD_DISABLED);
    int *restrict output_count_ptr = &(output_param->output_count_);
    const double output_interval = output_param->output_interval;
    double next_output_time = output_interval;

    double *restrict t_ptr = &(simulation_status->t);
    int64 *restrict num_steps_ptr = &(simulation_status->num_steps);

    const bool enable_progress_bar = settings->enable_progress_bar;

    /* Allocate memory */
    double *restrict x_0 = malloc(num_particles * 3 * sizeof(double));
    double *restrict v_0 = malloc(num_particles * 3 * sizeof(double));
    double *restrict vk1 = malloc(num_particles * 3 * sizeof(double));
    double *restrict vk2 = malloc(num_particles * 3 * sizeof(double));
    double *restrict vk3 = malloc(num_particles * 3 * sizeof(double));
    double *restrict vk4 = malloc(num_particles * 3 * sizeof(double));
    double *restrict xk1 = malloc(num_particles * 3 * sizeof(double));
    double *restrict xk2 = malloc(num_particles * 3 * sizeof(double));
    double *restrict xk3 = malloc(num_particles * 3 * sizeof(double));
    double *restrict xk4 = malloc(num_particles * 3 * sizeof(double));

    // Compensated summation
    double *restrict x_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict v_err_comp_sum = calloc(num_particles * 3, sizeof(double));

    // Check if memory allocation is successful
    if (!x_0 || !v_0 || !vk1 || !vk2 || !vk3 || !vk4 || !xk1 || !xk2 || !xk3 || !xk4 ||
        !x_err_comp_sum || !v_err_comp_sum)
    {
        THROW(CTB_MEMORY_ERROR, "Failed to allocate memory for arrays");
        goto err_memory;
    }

    /* Initial output */
    if (is_output && output_param->output_initial)
    {
        TRY_GOTO(output_snapshot(
            output_param,
            system,
            integrator_param,
            acceleration_param,
            simulation_status,
            settings
        ), err_initial_output);
        
    }

    /* Main Loop */
    int64 total_num_steps = (int64)ceil(tf / dt);
    ProgressBarParam progress_bar_param;
    if (enable_progress_bar)
    {
        TRY_GOTO(start_progress_bar(&progress_bar_param, total_num_steps), err_start_progress_bar);
    }

    *t_ptr = 0.0;
    simulation_status->dt = dt;
    *num_steps_ptr = 0;
    while (*num_steps_ptr < total_num_steps)
    {
        /* Store current state */
        memcpy(x_0, x, num_particles * 3 * sizeof(double));
        memcpy(v_0, v, num_particles * 3 * sizeof(double));

        /* Compute xk1 and vk1 */
        TRY_GOTO(acceleration(vk1, system, acceleration_param), err_acceleration);
        memcpy(xk1, v, num_particles * 3 * sizeof(double));

        /* Compute xk2 and vk2 */
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = x_0[i * 3 + j] + 0.5 * xk1[i * 3 + j] * dt;
                v[i * 3 + j] = v_0[i * 3 + j] + 0.5 * vk1[i * 3 + j] * dt;
            }
        }
        TRY_GOTO(acceleration(vk2, system, acceleration_param), err_acceleration);
        memcpy(xk2, v, num_particles * 3 * sizeof(double));

        /* Compute xk3 and vk3 */
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = x_0[i * 3 + j] + 0.5 * xk2[i * 3 + j] * dt;
                v[i * 3 + j] = v_0[i * 3 + j] + 0.5 * vk2[i * 3 + j] * dt;
            }
        }
        TRY_GOTO(acceleration(vk3, system, acceleration_param), err_acceleration);
        memcpy(xk3, v, num_particles * 3 * sizeof(double));

        /* Compute xk4 and vk4 */
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = x_0[i * 3 + j] + xk3[i * 3 + j] * dt;
                v[i * 3 + j] = v_0[i * 3 + j] + vk3[i * 3 + j] * dt;
            }
        }
        TRY_GOTO(acceleration(vk4, system, acceleration_param), err_acceleration);
        memcpy(xk4, v, num_particles * 3 * sizeof(double));

        /* Update step */
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                v_err_comp_sum[i * 3 + j] += (vk1[i * 3 + j] + 2 * vk2[i * 3 + j] +
                                              2 * vk3[i * 3 + j] + vk4[i * 3 + j]) *
                                             dt / 6.0;
                x_err_comp_sum[i * 3 + j] += (xk1[i * 3 + j] + 2 * xk2[i * 3 + j] +
                                              2 * xk3[i * 3 + j] + xk4[i * 3 + j]) *
                                             dt / 6.0;

                v[i * 3 + j] = v_0[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                x[i * 3 + j] = x_0[i * 3 + j] + x_err_comp_sum[i * 3 + j];

                v_err_comp_sum[i * 3 + j] += v_0[i * 3 + j] - v[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += x_0[i * 3 + j] - x[i * 3 + j];
            }
        }
        (*num_steps_ptr)++;
        *t_ptr = (*num_steps_ptr) * dt;

        /* Store solution */
        if (is_output && *t_ptr >= next_output_time)
        {
            TRY_GOTO(output_snapshot(
                output_param,
                system,
                integrator_param,
                acceleration_param,
                simulation_status,
                settings
            ), err_output);
            

            next_output_time = (*output_count_ptr) * output_interval;
        }

        if (enable_progress_bar)
        {
            update_progress_bar(&progress_bar_param, *num_steps_ptr, false);
        }

        /* Check exit */
        if (*(settings->is_exit_ptr))
        {
            break;
        }
    }

    if (enable_progress_bar)
    {
        update_progress_bar(&progress_bar_param, *num_steps_ptr, true);
    }

    free(x_0);
    free(v_0);
    free(vk1);
    free(vk2);
    free(vk3);
    free(vk4);
    free(xk1);
    free(xk2);
    free(xk3);
    free(xk4);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    return;

err_output:
err_acceleration:
err_start_progress_bar:
err_initial_output:
err_memory:
    free(x_0);
    free(v_0);
    free(vk1);
    free(vk2);
    free(vk3);
    free(vk4);
    free(xk1);
    free(xk2);
    free(xk3);
    free(xk4);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    return;
}

static void leapfrog(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    OutputParam *output_param,
    SimulationStatus *simulation_status,
    Settings *settings,
    const double tf
)
{
    /* Declare variables */
    

    const int num_particles = system->num_particles;
    double *restrict x = system->x;
    double *restrict v = system->v;

    double dt = integrator_param->dt;

    bool is_output = (output_param->method != OUTPUT_METHOD_DISABLED);
    int *restrict output_count_ptr = &(output_param->output_count_);
    const double output_interval = output_param->output_interval;
    double next_output_time = output_interval;

    double *restrict t_ptr = &(simulation_status->t);
    int64 *restrict num_steps_ptr = &(simulation_status->num_steps);

    const bool enable_progress_bar = settings->enable_progress_bar;

    /* Allocate memory */
    double *restrict temp_x = malloc(num_particles * 3 * sizeof(double));
    double *restrict temp_v = malloc(num_particles * 3 * sizeof(double));
    double *restrict a = malloc(num_particles * 3 * sizeof(double));

    // Compensated summation
    double *restrict x_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict v_err_comp_sum = calloc(num_particles * 3, sizeof(double));

    // Check if memory allocation is successful
    if (!temp_x || !temp_v || !a || !x_err_comp_sum || !v_err_comp_sum)
    {
        THROW(CTB_MEMORY_ERROR, "Failed to allocate memory for arrays");
        goto err_memory;
    }

    /* Initial output */
    if (is_output && output_param->output_initial)
    {
        TRY_GOTO(output_snapshot(
            output_param,
            system,
            integrator_param,
            acceleration_param,
            simulation_status,
            settings
        ), err_initial_output);
        
    }

    /* Compute initial acceleration and v_1/2 */
    TRY_GOTO(acceleration(a, system, acceleration_param), err_acceleration);

    memcpy(temp_v, v, num_particles * 3 * sizeof(double));
    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            v_err_comp_sum[i * 3 + j] += 0.5 * a[i * 3 + j] * dt;
            v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
            v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];
        }
    }

    /* Main Loop */
    int64 total_num_steps = (int64)ceil(tf / dt);
    ProgressBarParam progress_bar_param;
    if (enable_progress_bar)
    {
        TRY_GOTO(start_progress_bar(&progress_bar_param, total_num_steps), err_start_progress_bar);
    }

    *t_ptr = 0.0;
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

        /* Calculate x_1 */
        memcpy(temp_x, x, num_particles * 3 * sizeof(double));
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                x[i * 3 + j] = temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += temp_x[i * 3 + j] - x[i * 3 + j];
            }
        }

        /* Calculate v_1+1/2 */
        TRY_GOTO(acceleration(a, system, acceleration_param), err_acceleration);
        memcpy(temp_v, v, num_particles * 3 * sizeof(double));
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;
                v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];
            }
        }

        (*num_steps_ptr)++;
        *t_ptr = (*num_steps_ptr) * dt;

        /* Store solution */
        if (is_output && *t_ptr >= next_output_time)
        {
            // Get v_1 from v_1+1/2
            memcpy(temp_v, v, num_particles * 3 * sizeof(double));
            for (int i = 0; i < num_particles; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    v[i * 3 + j] -= 0.5 * a[i * 3 + j] * dt;
                }
            }
            TRY_GOTO(output_snapshot(
                output_param,
                system,
                integrator_param,
                acceleration_param,
                simulation_status,
                settings
            ), err_output);
            

            next_output_time = (*output_count_ptr) * output_interval;

            /* Restore v_1+1/2 */
            memcpy(v, temp_v, num_particles * 3 * sizeof(double));
        }

        if (enable_progress_bar)
        {
            update_progress_bar(&progress_bar_param, *num_steps_ptr, false);
        }

        /* Check exit */
        if (*(settings->is_exit_ptr))
        {
            break;
        }
    }

    /* Synchronize v_1+1/2 to v_1 */
    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            v[i * 3 + j] -= 0.5 * a[i * 3 + j] * dt;
        }
    }

    if (enable_progress_bar)
    {
        update_progress_bar(&progress_bar_param, *num_steps_ptr, true);
    }

    free(temp_x);
    free(temp_v);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    return;

err_output:
err_acceleration:
err_start_progress_bar:
err_initial_output:
err_memory:
    free(temp_x);
    free(temp_v);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    return;
}
