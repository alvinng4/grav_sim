import numpy as np
import matplotlib.pyplot as plt

# Default units is AU, days, and M_sun
DT = 1.0  # 1.0 days
PLOT_INTERVAL = 0.01 * 365.24  # 0.1 year to days
FPS = 60

SOLAR_SYSTEM_COLORS = {
    "Sun": "orange",
    "Mercury": "slategrey",
    "Venus": "wheat",
    "Earth": "skyblue",
    "Mars": "red",
    "Jupiter": "darkgoldenrod",
    "Saturn": "gold",
    "Uranus": "paleturquoise",
    "Neptune": "blue",
}
OBJECT_NAMES = list(SOLAR_SYSTEM_COLORS.keys())
COLORS = list(SOLAR_SYSTEM_COLORS.values())
MARKER_SIZE = [6.0, 1.5, 1.5, 2.0, 1.5, 4.0, 3.5, 3.5, 3.5]

X_LIM = (-2.0, 2.0)
Y_LIM = (-2.0, 2.0)
Z_LIM = (-2.0, 2.0)


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
        r_cm = np.sum(self.m[:, np.newaxis] * self.x, axis=0) / np.sum(self.m)
        v_cm = np.sum(self.m[:, np.newaxis] * self.v, axis=0) / np.sum(self.m)

        self.x -= r_cm
        self.v -= v_cm


def main() -> None:
    # Get initial conditions
    system = get_initial_conditions()

    # Initialize memory
    a = np.zeros((system.num_particles, 3))

    # Launch simulation
    print_simulation_info(system)
    next_plot_time = 0.0
    i = 0

    plt.ion()  # turn on interactive mode
    fig = plt.figure(figsize=(5, 5), dpi=150)
    plt.style.use("dark_background")
    ax = fig.add_subplot(111, projection="3d")

    while True:
        rk4(a, system, DT)

        current_time = i * DT
        if current_time >= next_plot_time:
            next_plot_time += PLOT_INTERVAL

            print(f"Current time: {current_time:.2f} days", end="\r")

            ax.grid(False)
            ax.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
            ax.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
            ax.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.zaxis.set_visible(False)
            ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

            for j in range(system.num_particles):
                ax.scatter(
                    system.x[j, 0],
                    system.x[j, 1],
                    system.x[j, 2],
                    marker="o",
                    color=COLORS[j],
                    s=MARKER_SIZE[j],
                    label=OBJECT_NAMES[j],
                )

            # Add legend
            legend = ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
            legend.facecolor = "transparent"

            # Adjust figure for the legend
            fig.subplots_adjust(right=0.7)
            fig.tight_layout()

            ax.set_xlim3d(X_LIM)
            ax.set_ylim3d(Y_LIM)
            ax.set_zlim3d(Z_LIM)

            # Set equal aspect ratio to prevent distortion
            ax.set_aspect("equal")

            plt.pause(1.0 / FPS)

            ax.clear()
        i += 1


def euler(a: np.ndarray, system: System, dt: float) -> None:
    """
    Advance one step with the Euler's method.

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array with shape (N, 3).
    system : System
        System object.
    dt : float
        Time step.
    """
    acceleration(a, system)
    system.x += system.v * dt
    system.v += a * dt


def euler_cromer(a: np.ndarray, system: System, dt: float) -> None:
    """
    Advance one step with the Euler-Cromer method.

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array with shape (N, 3).
    system : System
        System object.
    dt : float
        Time step.
    """
    acceleration(a, system)
    system.v += a * dt
    system.x += system.v * dt


def rk4(a: np.ndarray, system: System, dt: float) -> None:
    """
    Advance one step with the RK4 method.

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array with shape (N, 3).
    system : System
        System object.
    dt : float
        Time step.
    """
    num_stages = 4

    coeff = np.array([0.0, 0.5, 0.5, 1.0])
    weights = np.array([1.0, 2.0, 2.0, 1.0]) / 6.0

    # Allocate memory and initialize arrays
    x0 = system.x.copy()
    v0 = system.v.copy()
    xk = np.zeros((num_stages, system.num_particles, 3))
    vk = np.zeros((num_stages, system.num_particles, 3))

    # Initial stage
    acceleration(a, system)
    xk[0] = v0
    vk[0] = a

    # Compute the stages
    for i in range(1, num_stages):
        # Compute acceleration
        system.x = x0 + dt * coeff[i] * xk[i - 1]
        acceleration(a, system)

        # Compute xk and vk
        xk[i] = v0 + dt * coeff[i] * vk[i - 1]
        vk[i] = a

    # Advance step
    dx = 0.0
    dv = 0.0
    for i in range(num_stages):
        dx += weights[i] * xk[i]
        dv += weights[i] * vk[i]
    system.x = x0 + dt * dx
    system.v = v0 + dt * dv


def leapfrog(a: np.ndarray, system: System, dt: float) -> None:
    """
    Advance one step with the LeapFrog method.

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array with shape (N, 3).
    system : System
        System object.
    dt : float
        Time step.
    """
    # Velocity kick (v_1/2)
    acceleration(a, system)
    system.v += a * 0.5 * dt

    # Position drift (x_1)
    system.x += system.v * dt

    # Velocity kick (v_1)
    acceleration(a, system)
    system.v += a * 0.5 * dt


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


def print_simulation_info(system: System) -> None:
    print("----------------------------------------------------------")
    print("Simulation Info:")
    print(f"num_particles: {system.num_particles}")
    print(f"G: {system.G}")
    print(f"dt: {DT} days")
    print()
    print(f"Plot interval: {PLOT_INTERVAL} days")
    print("----------------------------------------------------------")


def plot_trajectory(
    sol_x: np.ndarray,
    labels: list,
    colors: list,
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
            color=traj[0].get_color(),
            label=labels[i],
        )

    fig.legend(loc="center right", borderaxespad=0.2)
    fig.tight_layout()
    plt.show()


def get_initial_conditions() -> System:
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


if __name__ == "__main__":
    main()
