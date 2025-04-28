# Step 1: Initial setup

Welcome to the first step of "12 steps to N-body simulation"! 
This is a series of tutorials with the goal of teaching beginners
how to write fast and clean N-body gravity simulations code in Python.
In this step, we will set up the required python environment and
implement some basic functions to be used in the following steps.

First of all, make sure you have Python 3.9 or later installed. You
can check the installation by running the following command in your terminal:
```
python --version
```

Two packages are required for this tutorial: `numpy` and `matplotlib`.
They are commonly used packages in Python for scientific computing and data visualization.
You can install them using `pip`:
```
pip install numpy matplotlib
```

Now we are ready to start coding. First, let us import the required packages:
```python
import numpy as np
import matplotlib.pyplot as plt
```

## System
To keep the code clean and organized, we will create a class called `System`
to represent the N-body system. It has the following attributes:

- `num_particles (float)`: Number of particles N in the system.
- `x (np.ndarray)`: Positions of the particles in 3D space (shape: (N, 3)).
- `v (np.ndarray)`: Velocities of the particles in 3D space (shape: (N, 3)).
- `m (np.ndarray)`: Masses of the particles (shape: (N,)).
- `G (float)`: Gravitational constant.

```python
class System:
    def __init__(
        self, num_particles: int, x: np.ndarray, v: np.ndarray, m: np.ndarray, G: float
    ) -> None:
        self.num_particles = num_particles
        self.x = x
        self.v = v
        self.m = m
        self.G = G
```

We will also implement a method to set the center of mass position and velocity to zero.
This is to prevent drifting of the system during the simulation, and to set the center of mass
of the system to the origin. To do this, we subtract the center of mass position
and velocity from the particles,

$$
    \mathbf{r}_{\mathrm{com}} = \frac{1}{M} \sum_{i=1}^{N} m_i \mathbf{r}_i,
    \quad 
    \mathbf{v}_{\mathrm{com}} = \frac{1}{M} \sum_{i=1}^{N} m_i \mathbf{v}_i.
$$

where $M$ is the total mass of the system, $m_i$, $\mathbf{r}_i$, and $\mathbf{v}_i$ are the mass, position, and velocity of the $i$-th particle respectively.
To compute the center of mass terms, we first compute $m_i \mathbf{r}_i$ by
```
self.m[:, np.newaxis] * self.x
```
where `np.newaxis` is used to "broadcast" the mass array to the shape of the position array.
By adding `np.newaxis`, the mass array is reshaped from `(N,)` to `(N, 1)`, i.e. expanded along
axis 1 (column).

$$
    \begin{bmatrix}
        m_{1} \\
        m_{2} \\
        \vdots \\
        m_{N}
    \end{bmatrix}
    \to
    \begin{bmatrix}
        m_{1} \dots \\
        m_{2} \dots \\
        \vdots  \\
        m_{N} \dots
    \end{bmatrix}
$$


The shape of `self.m[:, np.newaxis]` is now `(N, 1)` and the shape of `self.x` is `(N, 3)`.
The multiplication is then done element-wise as

$$
    \begin{bmatrix}
        m_{1} r_{1,1} & m_{1} r_{1,2} & m_{1} r_{1,3} \\
        m_{2} r_{2,1} & m_{2} r_{2,2} & m_{2} r_{2,3} \\
        \vdots & \vdots & \vdots \\
        m_{N} r_{N,1} & m_{N} r_{N,2} & m_{N} r_{N,3}
    \end{bmatrix}
$$

Then, we perform the summation along the axis 0 (row) with length `N` by
```
np.sum(self.m[:, np.newaxis] * self.x, axis=0)
```
Putting it all together, we have the following code:
```python 
class System:
    ...
    def center_of_mass_correction(self) -> None:
        """Set center of mass of position and velocity to zero"""
        r_cm = np.sum(self.m[:, np.newaxis] * self.x, axis=0) / np.sum(self.m)
        v_cm = np.sum(self.m[:, np.newaxis] * self.v, axis=0) / np.sum(self.m)

        self.x -= r_cm
        self.v -= v_cm
```

## Initial conditions (Solar System)
With the System class ready, we can now implement a function to get the initial conditions.
In our tutorial, we will use the solar system as our N-body system. The initial conditions
at 1/Jan/2024 are generated using the JPL Horizons System[@Horizons]. Gravitational constant
and masses of the solar system objects are calculated using the data from R.S. Park, et. al.[@Park2021].

Also, note that in our tutorial, we will stick with the following units:

- Length: AU (Astronomical Unit), i.e. the distance from the Earth to the Sun.
- Mass: $M_\odot$ (Solar mass)
- Time: days

They are convenient units for solar system simulations.
If you want to use different units, make sure to convert all data to the same units
and be consistent.

??? note "Code (Click to expand)"
    ```python
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
    ```

## Plotting initial conditions
Finally, we will implement a function to plot the initial conditions of the solar system.
We will use the `matplotlib` package to plot the positions of the particles in 2D.
Colors and labels are optional, but they make the plot look nicer.

If `plt.show()` does not work in your environment, you may need to change it to
`plt.savefig(file_name)` to save the plot.

```python hl_lines="52 53"
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
LABELS = list(SOLAR_SYSTEM_COLORS.keys())
COLORS = list(SOLAR_SYSTEM_COLORS.values())
LEGEND = True

def plot_initial_conditions(
    system: System,
    labels: list,
    colors: list,
    legend: bool,
) -> None:
    """
    Plot the initial positions.

    Parameters
    ----------
    system : System
        System object.
    labels : list
        Labels for the particles.
    colors : list
        Colors for the particles.
    legend : bool
        Whether to show the legend.
    """
    fig, ax = plt.subplots()
    ax.set_xlabel("$x$ (AU)")
    ax.set_ylabel("$y$ (AU)")

    for i in range(system.num_particles):
        ax.plot(
            system.x[i, 0],
            system.x[i, 1],
            marker = "o",
            color=colors[i],
            label=labels[i]
        )

    if legend:
        ax.legend()

    plt.show() # Here, you may need to change to plt.savefig(file_name) if 
               # plt.show() does not work in your environment. 
```

## Main function
Now we can put everything together in the `main` function.
```python
def main():
    # Get initial conditions
    system = get_initial_conditions()
    print("Number of particles:", system.num_particles)
    print("Initial positions (AU):\n", system.x)
    print("Initial velocities (AU/day):\n", system.v)
    print("Masses (M_sun):\n", system.m)
    print("Gravitational constant (AU^3 / day^2 / M_sun):", system.G)

    # Plot the initial conditions
    plot_initial_conditions(
        system,
        labels=LABELS,
        colors=COLORS,
        legend=LEGEND,
    )

if __name__ == "__main__":
    main()
```

Now, as you run the code, you should see:
```
Number of particles: 9
Initial positions (AU):
 [[-7.96712825e-03 -2.90611166e-03  2.10213120e-04]
 [-2.82597500e-01  1.97456095e-01  4.17742433e-02]
 [-7.23209543e-01 -7.94829045e-02  4.04286220e-02]
 [-1.73818374e-01  9.66324671e-01  1.55297876e-04]
 [-3.01325412e-01 -1.45402922e+00 -2.30054066e-02]
 [ 3.48520330e+00  3.55213702e+00 -9.27104467e-02]
 [ 8.98810505e+00 -3.71906474e+00 -2.93193870e-01]
 [ 1.22630250e+01  1.52973880e+01 -1.02054995e-01]
 [ 2.98350154e+01 -1.79381284e+00 -6.50640206e-01]]
Initial velocities (AU/day):
 [[ 4.87524241e-06 -7.05716139e-06 -4.57929038e-08]
 [-2.23216589e-02 -2.15720711e-02  2.85519283e-04]
 [ 2.03406835e-03 -2.02082863e-02 -3.94564043e-04]
 [-1.72300122e-02 -2.96772137e-03  6.38154172e-07]
 [ 1.42483227e-02 -1.57923621e-03 -3.82372338e-04]
 [-5.47097051e-03  5.64248731e-03  9.89618477e-05]
 [ 1.82201399e-03  5.14347040e-03 -1.61723649e-04]
 [-3.09761521e-03  2.27678190e-03  4.86042739e-05]
 [ 1.67653809e-04  3.15209870e-03 -6.87750693e-05]]
Masses (M_sun):
 [1.00000000e+00 1.66012083e-07 2.44783829e-06 3.00348962e-06
 3.22715608e-07 9.54791910e-04 2.85885670e-04 4.36624961e-05
 5.15138377e-05]
Gravitational constant (AU^3 / day^2 / M_sun): 0.00029591220828411956
```
You should also see the following plot:

<img src="../figures/step1/initial_conditions.png" alt="Initial Condition" width="500">


## Full script
The full script is available at `6_steps_to_n_body_simulation/python/step1.py`,
or https://github.com/alvinng4/grav_sim/blob/main/6_steps_to_n_body_simulation/python/step1.py

??? note "Code (Click to expand)"
    ```python linenums="1"
    --8<-- "6_steps_to_n_body_simulation/python/step1.py"
    ```