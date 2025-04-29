import math
import timeit

import numpy as np
import matplotlib.pyplot as plt

# Default units is AU, days, and M_sun
# INITIAL_CONDITION = "solar_system"
# TF = 200.0 * 365.24  # 200 years to days
# TOLERANCE = 1e-8
# OUTPUT_INTERVAL = 0.01 * 365.24  # 0.01 year to days
# INITIAL_DT = 1.0  # Initial time step in days
# SOLAR_SYSTEM_COLORS = {
#     "Sun": "orange",
#     "Mercury": "slategrey",
#     "Venus": "wheat",
#     "Earth": "skyblue",
#     "Mars": "red",
#     "Jupiter": "darkgoldenrod",
#     "Saturn": "gold",
#     "Uranus": "paleturquoise",
#     "Neptune": "blue",
# }
# LABELS = list(SOLAR_SYSTEM_COLORS.keys())
# COLORS = list(SOLAR_SYSTEM_COLORS.values())
# LEGEND = True

# Pyth-3-body
INITIAL_CONDITION = "pyth-3-body"
TF = 70.0  # 70 days
TOLERANCE = 1e-13
OUTPUT_INTERVAL = 0.001  # 0.001 day
INITIAL_DT = 0.01  # Initial time step in days
LABELS = [None, None, None]
COLORS = [None, None, None]
LEGEND = False


class System:
    def __init__(
        self, num_particles: int, x: np.ndarray, v: np.ndarray, m: np.ndarray, G: float
    ) -> None:
        self.num_particles = num_particles
        self.x = x
        self.v = v
        self.m = m
        self.G = G

    def center_of_mass_correction(self) -> None:
        """Set center of mass of position and velocity to zero"""
        M = np.sum(self.m)
        x_cm = np.einsum("i,ij->j", self.m, self.x) / M
        v_cm = np.einsum("i,ij->j", self.m, self.v) / M

        self.x -= x_cm
        self.v -= v_cm


def main() -> None:
    # Get initial conditions
    system = get_initial_conditions(INITIAL_CONDITION)

    # RKF4(5) coefficients
    # fmt: off
    coeff = np.array((
        [1.0 / 4.0, 0.0, 0.0, 0.0, 0.0],
        [3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0],
        [1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0],
        [439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0],
        [-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0],
    ))
    # fmt: on
    weights = np.array(
        [25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0.0]
    )
    weights_test = np.array(
        [
            16.0 / 135.0,
            0.0,
            6656.0 / 12825.0,
            28561.0 / 56430.0,
            -9.0 / 50.0,
            2.0 / 55.0,
        ]
    )
    min_power = 4
    num_stages = len(weights)

    # Initialize memory and arrays
    a = np.zeros((system.num_particles, 3))
    temp_system = System(
        num_particles=system.num_particles,
        x=np.zeros((system.num_particles, 3)),
        v=np.zeros((system.num_particles, 3)),
        m=system.m,
        G=system.G,
    )
    x_1 = np.zeros((system.num_particles, 3))
    v_1 = np.zeros((system.num_particles, 3))
    xk = np.zeros((num_stages, system.num_particles, 3))
    vk = np.zeros((num_stages, system.num_particles, 3))
    error_estimation_delta_x = np.zeros((system.num_particles, 3))
    error_estimation_delta_v = np.zeros((system.num_particles, 3))
    tolerance_scale_x = np.zeros((system.num_particles, 3))
    tolerance_scale_v = np.zeros((system.num_particles, 3))

    # Safety factors for step-size control
    safety_fac_max = 6.0
    safety_fac_min = 0.33
    safety_fac = math.pow(0.38, 1.0 / (1.0 + float(min_power)))

    # Solution array
    sol_size = int(TF // OUTPUT_INTERVAL + 2)  # +2 for initial and final time
    sol_x = np.zeros((sol_size, system.num_particles, 3))
    sol_v = np.zeros((sol_size, system.num_particles, 3))
    sol_t = np.zeros(sol_size)
    sol_dt = np.zeros(sol_size)
    sol_x[0] = system.x
    sol_v[0] = system.v
    sol_t[0] = 0.0
    sol_dt[0] = INITIAL_DT
    output_count = 1

    # Launch simulation
    print_simulation_info(system, sol_size)
    next_output_time = output_count * OUTPUT_INTERVAL
    start = timeit.default_timer()
    dt = INITIAL_DT
    current_time = 0.0
    while current_time < TF:
        # Initial stage
        acceleration(a, system)
        xk[0] = system.v
        vk[0] = a

        # Compute the stages
        for stage in range(1, num_stages):
            # Empty temp_x and temp_v
            temp_system.x.fill(0.0)
            temp_system.v.fill(0.0)

            for i in range(stage):
                temp_system.x[:] += coeff[stage - 1, i] * xk[i]
                temp_system.v[:] += coeff[stage - 1, i] * vk[i]

            temp_system.x[:] = system.x + dt * temp_system.x
            temp_system.v[:] = system.v + dt * temp_system.v

            # Compute the acceleration
            xk[stage] = temp_system.v
            acceleration(vk[stage], temp_system)

        # Calculate x_1, v_1 and also delta x, delta v for error estimation
        x_1[:] = system.x
        v_1[:] = system.v
        error_estimation_delta_x.fill(0.0)
        error_estimation_delta_v.fill(0.0)
        for stage in range(num_stages):
            x_1[:] += dt * weights[stage] * xk[stage]
            v_1[:] += dt * weights[stage] * vk[stage]
            error_estimation_delta_x[:] += (
                dt * (weights[stage] - weights_test[stage]) * xk[stage]
            )
            error_estimation_delta_v[:] += (
                dt * (weights[stage] - weights_test[stage]) * vk[stage]
            )

        # Error estimation
        tolerance_scale_x[:] = (
            TOLERANCE + np.maximum(np.abs(system.x), np.abs(x_1)) * TOLERANCE
        )
        tolerance_scale_v[:] = (
            TOLERANCE + np.maximum(np.abs(system.v), np.abs(v_1)) * TOLERANCE
        )

        sum = np.sum(np.square(error_estimation_delta_x / tolerance_scale_x)) + np.sum(
            np.square(error_estimation_delta_v / tolerance_scale_v)
        )
        error = math.sqrt(sum / (system.num_particles * 3.0 * 2.0))

        # Advance step
        if error <= 1.0 or dt <= TF * 1e-12:
            current_time += dt
            system.x[:] = x_1
            system.v[:] = v_1

            if current_time >= next_output_time:
                sol_x[output_count] = system.x
                sol_v[output_count] = system.v
                sol_t[output_count] = current_time
                sol_dt[output_count] = dt

                output_count += 1
                next_output_time = output_count * OUTPUT_INTERVAL

                print(f"Current time: {current_time:.2f} days", end="\r")

        # Calculate dt for next step
        if error < 1e-12:
            error = 1e-12  # Prevent error from being too small

        dt_new = dt * safety_fac / math.pow(error, 1.0 / (1.0 + float(min_power)))
        if dt_new > safety_fac_max * dt:
            dt *= safety_fac_max
        elif dt_new < safety_fac_min * dt:
            dt *= safety_fac_min
        else:
            dt = dt_new

        if dt_new < TF * 1e-12:
            dt = TF * 1e-12

        # Correct overshooting
        if current_time < TF and current_time + dt > TF:
            dt = TF - current_time

    sol_x = sol_x[:output_count]
    sol_v = sol_v[:output_count]
    sol_t = sol_t[:output_count]
    sol_dt = sol_dt[:output_count]

    end = timeit.default_timer()

    print()
    print(f"Done! Runtime: {end - start:.3g} seconds, Solution size: {output_count}")
    plot_trajectory(
        sol_x=sol_x,
        labels=LABELS,
        colors=COLORS,
        legend=LEGEND,
    )

    # Compute and plot relative energy error
    rel_energy_error = compute_rel_energy_error(sol_x, sol_v, system)
    print(f"Relative energy error: {rel_energy_error[-1]:.3g}")
    plot_rel_energy_error(rel_energy_error, sol_t / 365.24)
    plot_dt(sol_dt, sol_t)


def compute_rel_energy_error(
    sol_x: np.ndarray, sol_v: np.ndarray, system: System
) -> np.ndarray:
    """
    Compute the relative energy error of the simulation.

    Parameters
    ----------
    sol_x : np.ndarray
        Solution position array with shape (N_steps, num_particles, 3).
    sol_v : np.ndarray
        Solution velocity array with shape (N_steps, num_particles, 3).
    system : System
        System object.

    Returns
    -------
    energy_error : np.ndarray
        Relative energy error of the simulation, with shape (N_steps,).
    """
    # Allocate memory and initialize arrays
    n_steps = sol_x.shape[0]
    num_particles = system.num_particles
    m = system.m
    G = system.G
    rel_energy_error = np.zeros(n_steps)

    # Compute the total energy (KE + PE)
    for count in range(n_steps):
        x = sol_x[count]
        v = sol_v[count]
        for i in range(num_particles):
            # KE
            rel_energy_error[count] += 0.5 * m[i] * np.linalg.norm(v[i]) ** 2
            # PE
            for j in range(i + 1, num_particles):
                rel_energy_error[count] -= G * m[i] * m[j] / np.linalg.norm(x[i] - x[j])

    # Compute the relative energy error
    initial_energy = rel_energy_error[0]
    rel_energy_error = (rel_energy_error - initial_energy) / initial_energy
    rel_energy_error = np.abs(rel_energy_error)

    return rel_energy_error


def plot_rel_energy_error(rel_energy_error: np.ndarray, sol_t: np.ndarray) -> None:
    """
    Plot the relative energy error.

    Parameters
    ----------
    rel_energy_error : np.ndarray
        Relative energy error of the simulation, with shape (N_steps,).
    sol_t : np.ndarray
        Solution time array with shape (N_steps,).
    """
    plt.figure()
    plt.plot(sol_t, rel_energy_error)
    plt.yscale("log")
    plt.xlabel("Time step")
    plt.ylabel("Relative Energy Error")
    plt.title("Relative Energy Error vs Time Step")
    plt.show()


def plot_dt(sol_dt: np.ndarray, sol_t: np.ndarray) -> None:
    """
    Plot the time step.

    Parameters
    ----------
    sol_dt : np.ndarray
        Time step array with shape (N_steps,).
    sol_t : np.ndarray
        Solution time array with shape (N_steps,).
    """
    plt.figure()
    plt.semilogy(sol_t, sol_dt)
    plt.xlabel("Time step")
    plt.ylabel("Time step (days)")
    plt.title("Time Step vs Time Step")
    plt.show()


def acceleration(
    a: np.ndarray,
    system: System,
) -> None:
    """
    Compute the gravitational acceleration

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array to be modified in-place,
        with shape (N, 3)
    system : System
        System object.
    """
    # Empty acceleration array
    a.fill(0.0)

    # Declare variables
    x = system.x
    m = system.m
    G = system.G

    # Compute the displacement vector
    r_ij = x[:, np.newaxis, :] - x[np.newaxis, :, :]

    # Compute the distance
    r_norm = np.linalg.norm(r_ij, axis=2)

    # Compute 1 / r^3
    with np.errstate(divide="ignore", invalid="ignore"):
        inv_r_cubed = 1.0 / (r_norm * r_norm * r_norm)

    # Set diagonal elements to 0 to avoid self-interaction
    np.fill_diagonal(inv_r_cubed, 0.0)

    # Compute the acceleration
    a[:] = G * np.einsum("ijk,ij,i->jk", r_ij, inv_r_cubed, m)


def print_simulation_info(system: System, sol_size: int) -> None:
    print("----------------------------------------------------------")
    print("Simulation Info:")
    print(f"num_particles: {system.num_particles}")
    print(f"G: {system.G}")
    print(f"tf: {TF} days")
    print(f"tolerance: {TOLERANCE}")
    print(f"Initial dt: {INITIAL_DT} days")
    print()
    print(f"Output interval: {OUTPUT_INTERVAL} days")
    print(f"Estimated solution size: {sol_size}")
    print("----------------------------------------------------------")


def plot_trajectory(
    sol_x: np.ndarray,
    labels: list,
    colors: list,
    legend: bool,
) -> None:
    """
    Plot the 2D trajectory.

    Parameters
    ----------
    sol_x : np.ndarray
        Solution position array with shape (N_steps, num_particles, 3).
    labels : list
        List of labels for the particles.
    colors : list
        List of colors for the particles.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")
    ax.set_xlabel("$x$ (AU)")
    ax.set_ylabel("$y$ (AU)")

    for i in range(sol_x.shape[1]):
        traj = ax.plot(
            sol_x[:, i, 0],
            sol_x[:, i, 1],
            color=colors[i],
        )
        # Plot the last position with marker
        ax.plot(
            sol_x[-1, i, 0],
            sol_x[-1, i, 1],
            marker="o",
            color=traj[0].get_color(),
            label=labels[i],
        )

    if legend:
        fig.legend(loc="center right", borderaxespad=0.2)
        fig.tight_layout()
    plt.show()


def get_initial_conditions(initial_condition: str) -> System:
    """
    Returns the initial conditions for solar system,
    with units AU, days, and M_sun.

    Returns
    -------
    x : np.ndarray
        Initial positions of the solar system bodies,
        with shape (N, 3).
    v : np.ndarray
        Initial velocities of the solar system bodies,
        with shape (N, 3).
    m : np.ndarray
        Masses of the solar system bodies,
        with shape (N,).
    G : float
        Gravitational constant.
    """
    # Conversion factor from km^3 s^-2 to AU^3 d^-2
    CONVERSION_FACTOR = (86400**2) / (149597870.7**3)

    # GM values (km^3 s^-2)
    # ref: https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf
    GM_KM_S = {
        "Sun": 132712440041.279419,
        "Mercury": 22031.868551,
        "Venus": 324858.592000,
        "Earth": 398600.435507,
        "Mars": 42828.375816,
        "Jupiter": 126712764.100000,
        "Saturn": 37940584.841800,
        "Uranus": 5794556.400000,
        "Neptune": 6836527.100580,
        "Moon": 4902.800118,
        "Pluto": 975.500000,
        "Ceres": 62.62890,
        "Vesta": 17.288245,
    }

    # GM values (AU^3 d^-2)
    GM_AU_DAY = {
        "Sun": 132712440041.279419 * CONVERSION_FACTOR,
        "Mercury": 22031.868551 * CONVERSION_FACTOR,
        "Venus": 324858.592000 * CONVERSION_FACTOR,
        "Earth": 398600.435507 * CONVERSION_FACTOR,
        "Mars": 42828.375816 * CONVERSION_FACTOR,
        "Jupiter": 126712764.100000 * CONVERSION_FACTOR,
        "Saturn": 37940584.841800 * CONVERSION_FACTOR,
        "Uranus": 5794556.400000 * CONVERSION_FACTOR,
        "Neptune": 6836527.100580 * CONVERSION_FACTOR,
        "Moon": 4902.800118 * CONVERSION_FACTOR,
        "Pluto": 975.500000 * CONVERSION_FACTOR,
        "Ceres": 62.62890 * CONVERSION_FACTOR,
        "Vesta": 17.288245 * CONVERSION_FACTOR,
    }

    # Solar system masses (M_sun^-1)
    SOLAR_SYSTEM_MASSES = {
        "Sun": 1.0,
        "Mercury": GM_KM_S["Mercury"] / GM_KM_S["Sun"],
        "Venus": GM_KM_S["Venus"] / GM_KM_S["Sun"],
        "Earth": GM_KM_S["Earth"] / GM_KM_S["Sun"],
        "Mars": GM_KM_S["Mars"] / GM_KM_S["Sun"],
        "Jupiter": GM_KM_S["Jupiter"] / GM_KM_S["Sun"],
        "Saturn": GM_KM_S["Saturn"] / GM_KM_S["Sun"],
        "Uranus": GM_KM_S["Uranus"] / GM_KM_S["Sun"],
        "Neptune": GM_KM_S["Neptune"] / GM_KM_S["Sun"],
        "Moon": GM_KM_S["Moon"] / GM_KM_S["Sun"],
        "Pluto": GM_KM_S["Pluto"] / GM_KM_S["Sun"],
        "Ceres": GM_KM_S["Ceres"] / GM_KM_S["Sun"],
        "Vesta": GM_KM_S["Vesta"] / GM_KM_S["Sun"],
    }

    G = GM_AU_DAY["Sun"]

    # Solar system position and velocities data
    # Units: AU-D
    # Coordinate center: Solar System Barycenter
    # Data dated on A.D. 2024-Jan-01 00:00:00.0000 TDB
    # Computational data generated by NASA JPL Horizons System https://ssd.jpl.nasa.gov/horizons/
    SOLAR_SYSTEM_POS = {
        "Sun": [-7.967955691533730e-03, -2.906227441573178e-03, 2.103054301547123e-04],
        "Mercury": [
            -2.825983269538632e-01,
            1.974559795958082e-01,
            4.177433558063677e-02,
        ],
        "Venus": [
            -7.232103701666379e-01,
            -7.948302026312400e-02,
            4.042871428174315e-02,
        ],
        "Earth": [-1.738192017257054e-01, 9.663245550235138e-01, 1.553901854897183e-04],
        "Mars": [-3.013262392582653e-01, -1.454029331393295e00, -2.300531433991428e-02],
        "Jupiter": [3.485202469657674e00, 3.552136904413157e00, -9.271035442798399e-02],
        "Saturn": [8.988104223143450e00, -3.719064854634689e00, -2.931937777323593e-01],
        "Uranus": [1.226302417897505e01, 1.529738792480545e01, -1.020549026883563e-01],
        "Neptune": [
            2.983501460984741e01,
            -1.793812957956852e00,
            -6.506401132254588e-01,
        ],
        "Moon": [-1.762788124769829e-01, 9.674377513177153e-01, 3.236901585768862e-04],
        "Pluto": [1.720200478843485e01, -3.034155683573043e01, -1.729127607100611e00],
        "Ceres": [-1.103880510367569e00, -2.533340440444230e00, 1.220283937721780e-01],
        "Vesta": [-8.092549658731499e-02, 2.558381434460076e00, -6.695836142398572e-02],
    }
    SOLAR_SYSTEM_VEL = {
        "Sun": [4.875094764261564e-06, -7.057133213976680e-06, -4.573453713094512e-08],
        "Mercury": [
            -2.232165900189702e-02,
            -2.157207103176252e-02,
            2.855193410495743e-04,
        ],
        "Venus": [
            2.034068201002341e-03,
            -2.020828626592994e-02,
            -3.945639843855159e-04,
        ],
        "Earth": [
            -1.723001232538228e-02,
            -2.967721342618870e-03,
            6.382125383116755e-07,
        ],
        "Mars": [1.424832259345280e-02, -1.579236181580905e-03, -3.823722796161561e-04],
        "Jupiter": [
            -5.470970658852281e-03,
            5.642487338479145e-03,
            9.896190602066252e-05,
        ],
        "Saturn": [
            1.822013845554067e-03,
            5.143470425888054e-03,
            -1.617235904887937e-04,
        ],
        "Uranus": [
            -3.097615358317413e-03,
            2.276781932345769e-03,
            4.860433222241686e-05,
        ],
        "Neptune": [
            1.676536611817232e-04,
            3.152098732861913e-03,
            -6.877501095688201e-05,
        ],
        "Moon": [
            -1.746667306153906e-02,
            -3.473438277358121e-03,
            -3.359028758606074e-05,
        ],
        "Pluto": [2.802810313667557e-03, 8.492056438614633e-04, -9.060790113327894e-04],
        "Ceres": [
            8.978653480111301e-03,
            -4.873256528198994e-03,
            -1.807162046049230e-03,
        ],
        "Vesta": [
            -1.017876585480054e-02,
            -5.452367109338154e-04,
            1.255870551153315e-03,
        ],
    }

    if initial_condition == "pyth-3-body":
        # Pythagorean 3-body problem
        R1 = np.array([1.0, 3.0, 0.0])
        R2 = np.array([-2.0, -1.0, 0.0])
        R3 = np.array([1.0, -1.0, 0.0])
        V1 = np.array([0.0, 0.0, 0.0])
        V2 = np.array([0.0, 0.0, 0.0])
        V3 = np.array([0.0, 0.0, 0.0])

        x = np.array([R1, R2, R3])
        v = np.array([V1, V2, V3])
        m = np.array([3.0 / G, 4.0 / G, 5.0 / G])

        system = System(
            num_particles=len(m),
            x=x,
            v=v,
            m=m,
            G=G,
        )
        system.center_of_mass_correction()

        return system

    elif initial_condition == "solar_system":
        m = np.array(
            [
                SOLAR_SYSTEM_MASSES["Sun"],
                SOLAR_SYSTEM_MASSES["Mercury"],
                SOLAR_SYSTEM_MASSES["Venus"],
                SOLAR_SYSTEM_MASSES["Earth"],
                SOLAR_SYSTEM_MASSES["Mars"],
                SOLAR_SYSTEM_MASSES["Jupiter"],
                SOLAR_SYSTEM_MASSES["Saturn"],
                SOLAR_SYSTEM_MASSES["Uranus"],
                SOLAR_SYSTEM_MASSES["Neptune"],
            ]
        )

        R1 = np.array(SOLAR_SYSTEM_POS["Sun"])
        R2 = np.array(SOLAR_SYSTEM_POS["Mercury"])
        R3 = np.array(SOLAR_SYSTEM_POS["Venus"])
        R4 = np.array(SOLAR_SYSTEM_POS["Earth"])
        R5 = np.array(SOLAR_SYSTEM_POS["Mars"])
        R6 = np.array(SOLAR_SYSTEM_POS["Jupiter"])
        R7 = np.array(SOLAR_SYSTEM_POS["Saturn"])
        R8 = np.array(SOLAR_SYSTEM_POS["Uranus"])
        R9 = np.array(SOLAR_SYSTEM_POS["Neptune"])

        V1 = np.array(SOLAR_SYSTEM_VEL["Sun"])
        V2 = np.array(SOLAR_SYSTEM_VEL["Mercury"])
        V3 = np.array(SOLAR_SYSTEM_VEL["Venus"])
        V4 = np.array(SOLAR_SYSTEM_VEL["Earth"])
        V5 = np.array(SOLAR_SYSTEM_VEL["Mars"])
        V6 = np.array(SOLAR_SYSTEM_VEL["Jupiter"])
        V7 = np.array(SOLAR_SYSTEM_VEL["Saturn"])
        V8 = np.array(SOLAR_SYSTEM_VEL["Uranus"])
        V9 = np.array(SOLAR_SYSTEM_VEL["Neptune"])

        x = np.array(
            [
                R1,
                R2,
                R3,
                R4,
                R5,
                R6,
                R7,
                R8,
                R9,
            ]
        )
        v = np.array(
            [
                V1,
                V2,
                V3,
                V4,
                V5,
                V6,
                V7,
                V8,
                V9,
            ]
        )

        system = System(
            num_particles=len(m),
            x=x,
            v=v,
            m=m,
            G=G,
        )
        system.center_of_mass_correction()

        return system

    else:
        raise ValueError(f"Initial condition not recognized: {initial_condition}.")


if __name__ == "__main__":
    main()
