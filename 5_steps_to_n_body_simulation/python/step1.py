import common

INITIAL_CONDITION = "solar_system"


def main():
    # Get initial conditions
    system, labels, colors, legend = common.get_initial_conditions(INITIAL_CONDITION)
    print("Number of particles:", system.num_particles)
    print("Initial positions (AU):\n", system.x)
    print("Initial velocities (AU/day):\n", system.v)
    print("Masses (M_sun):\n", system.m)
    print("Gravitational constant (AU^3 / day^2 / M_sun):", system.G)

    # Plot the initial conditions
    common.plot_initial_conditions(
        system=system,
        labels=labels,
        colors=colors,
        legend=legend,
    )


if __name__ == "__main__":
    main()
