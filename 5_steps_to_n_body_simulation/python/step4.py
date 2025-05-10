import timeit

import numpy as np

import common

INITIAL_CONDITION = "solar_system"

# Default units is AU, days, and M_sun
TF = 200.0 * 365.24  # years to days
DT = 1.0
OUTPUT_INTERVAL = 0.1 * 365.24  # years to days
NUM_STEPS = int(TF / DT)


def main() -> None:
    # Get initial conditions
    system, labels, colors, legend = common.get_initial_conditions(INITIAL_CONDITION)

    # Initialize memory
    a = np.zeros((system.num_particles, 3))

    # Solution array
    sol_size = int(TF // OUTPUT_INTERVAL + 2)  # +2 for initial and final time
    sol_x = np.zeros((sol_size, system.num_particles, 3))
    sol_v = np.zeros((sol_size, system.num_particles, 3))
    sol_t = np.zeros(sol_size)
    sol_x[0] = system.x
    sol_v[0] = system.v
    sol_t[0] = 0.0
    output_count = 1

    # Launch simulation
    common.print_simulation_info_fixed_step_size(
        system, TF, DT, NUM_STEPS, OUTPUT_INTERVAL, sol_size
    )
    next_output_time = output_count * OUTPUT_INTERVAL
    start = timeit.default_timer()
    for i in range(NUM_STEPS):
        common.leapfrog(a, system, DT)

        current_time = i * DT
        if current_time >= next_output_time:
            sol_x[output_count] = system.x
            sol_v[output_count] = system.v
            sol_t[output_count] = current_time

            output_count += 1
            next_output_time = output_count * OUTPUT_INTERVAL

            print(f"Current time: {current_time:.2f} days", end="\r")

    sol_x = sol_x[:output_count]
    sol_v = sol_v[:output_count]
    sol_t = sol_t[:output_count]

    end = timeit.default_timer()

    print()
    print(f"Done! Runtime: {end - start:.3g} seconds, Solution size: {output_count}")
    common.plot_trajectory(
        sol_x=sol_x,
        labels=labels,
        colors=colors,
        legend=legend,
    )

    # Compute and plot relative energy error
    rel_energy_error = common.compute_rel_energy_error(sol_x, sol_v, system)
    print(f"Relative energy error: {rel_energy_error[-1]:.3g}")
    common.plot_rel_energy_error(rel_energy_error, sol_t / 365.24)


if __name__ == "__main__":
    main()
