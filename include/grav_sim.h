#ifndef GRAV_SIM_H
#define GRAV_SIM_H

/* Acceleration */
#include "grav_sim/acceleration.h"

/* Common definitions */
#include "grav_sim/common.h"

/* Exception handling */
#include "grav_sim/error.h"

/* Integrator */
#include "grav_sim/integrator.h"

/* Output */
#include "grav_sim/output.h"

/* System */
#include "grav_sim/system.h"

/* Settings */
#include "grav_sim/settings.h"

/* Utils */
#include "grav_sim/utils.h"


/* Project version */
#ifndef VERSION_INFO
#define VERSION_INFO "unknown"
#endif

/**
 * \brief Main function to launch a simulation.
 * 
 * \param system Pointer to the system.
 * \param integrator_param Pointer to the integrator parameters.
 * \param acceleration_param Pointer to the acceleration parameters.
 * \param output_param Pointer to the output parameters.
 * \param simulation_status Pointer to the simulation status.
 * \param settings Pointer to the settings.
 * \param tf Simulation time.
 * 
 * \return Error code.
 */
int launch_simulation(
    System *restrict system,
    IntegratorParam *restrict integrator_param,
    AccelerationParam *restrict acceleration_param,
    OutputParam *restrict output_param,
    SimulationStatus *restrict simulation_status,
    Settings *restrict settings,
    const double tf
);

/**
 * \brief Main function to launch a cosmological simulation.
 * 
 * \param system Pointer to the system.
 * \param output_param Pointer to the output parameters.
 * \param simulation_status Pointer to the simulation status.
 * \param settings Pointer to the settings.
 * \param a_final Final scale factor. 
 * \param num_steps Number of steps.
 * \param pm_grid_size Particle mesh grid size.
 * 
 * \return Error code.
 */
int launch_cosmological_simulation(
    CosmologicalSystem *restrict system,
    OutputParam *restrict output_param,
    SimulationStatus *restrict simulation_status,
    Settings *restrict settings,
    const double a_final,
    const int num_steps,
    const int pm_grid_size
);

/**
 * \brief Get the logo string of grav_sim.
 * 
 * \return Pointer to the logo string.
 */
const char* get_grav_sim_logo_string(void);

/**
 * \brief Print project compilation information.
 */
void print_compilation_info(void);

#endif
