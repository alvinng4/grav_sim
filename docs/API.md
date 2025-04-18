## Documentation
* [Quick Start](#quick-start)
    - [Prerequisite](#prerequisite)
    - [Installation](#installation)
    - [Some notes](#some-notes)
* [GravitySimulator API](#gravitysimulator-api)
* [Default systems](#default-systems)
* [Integrators](#integrators)
    - [Simple methods](#simple-methods)
    - [Embedded Runge-Kutta methods](#embdedded-runge-kutta-methods)
    - [IAS15](#IAS15)
    - [WHFast](#whfast)
* [Saving the results](#saving-the-results)
* [Output animations in .gif](#output-animations-in-gif)

## Quick Start

### Prerequisite

1. Python version 3.10 or higher. 
2. Any C compiler that supports C99

### Installation
1. Download the source files, or clone this repository by
    ```
    git clone https://github.com/alvinng4/Gravity-Simulator
    ```
2. Install the required packages by
    ```
    pip install .
    ```
    If the installation is not successful, install the following packages manually:
    ```
    matplotlib==3.8.3
    numpy==1.26.4
    rich==13.7.1
    Pillow==10.3.0
    ```
3. Compile the C library

    I have provided a compilation of the C library in the repository. If it does not run on your computer,
    you may need to recompile it.

    To compile the C library, simply go to the src folder and run
    ```
    make [CC=gcc] [USE_OPENMP=1] [USE_CUDA=1]
    ```
    `CC`: To indicate which C compiler to use.
    `USE_OPENMP=1`: To enable OpenMP acceleration.
    `USE_CUDA=1`: To enable CUDA acceleration.

    Note:
    - If the program is compiled with openmp, the program will run with openmp by default, which could be slow if
    $N$ is small. Use `export OMP_NUM_THREADS=1` to disable openmp.

### Some notes
* The default unit for this project is solar masses, AU and days, with $G = 0.00029591220828411956 \text{ M}_\odot^{-1} \text{ AU}^3 \text{ day}^{-2}$.
It is possible to change this value in the API by changing `system.G`.
* Check the `examples` folder for API tutorial and sample projects
* For WHFast, features including Barnes-Hut algorithm are not available due to implementation difficulties.

## GravitySimulator API

You may import the GravitySimulator API from `gravity_sim` to perform gravity simulation.
See `examples/tutorial.ipynb` or [Sample projects](#sample-projects) for some example usage.
If your computer cannot render jupyter notebook (files end with `.ipynb`), you can view them on github.

<!-- #### launch_simulation
launch_simulation() is the main method for launching the simulation.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `integrator` | str | Required | Name of the integrator |
| `tf` | float | Required | Simulation time (days) |
| `dt` | float | None | Time step (days) |
| `tolerance` | float | None | Tolerance for adaptive step integrators |
| `store_every_n` | int | 1 | Store results every n steps |
| `acceleration_method` | str | `pairwise` | Method for calculating accelerations |
| `storing_method` | str | `default` | Method for storing simulation results |
| `flush_path` | str | None | Path to flush intermediate results. |
| `no_progress_bar` | bool | False | If True, disables the progress bar. |
| `no_print` | bool | False | If True, disables some printing to console |
| `softening_length` | float | 0.0 | Softening length for acceleration calculation |
| `**kwargs` | dict | - | Additional keyword arguments. | -->

#### integrators 
`euler`, `euler_cromer`, `rk4`, `leapfrog`, `rkf45`, `dopri`, `dverk`, `rkf78`, `ias15`, `whfast`

#### acceleration_method
- `pairwise`
    * Brute force pairwise calculations for gravitational acceleration
    * Time complexity: $O(N^2)$
- `massless`
    * Similar to `pairwise`, but seperate the calculations for massive and massless particles
    * Time complexity: $O(M^2 + MN)$, where $M$ and $N$ are the number of massive and massless particles respectively
- `barnes_hut`
    * Calculate gravitational acceleration with Barnes-Hut algorithm
    * Time complexity: $O(N \log{N})$
    * `**kwargs`: `opening_angle`
        * Threshold for Barnes-Hut algorithm, default = 0.5

<!-- - `fast_multipole`
    * Calculate gravitational acceleration with fast multipole method (FMM)
        * Time complexity: $O(N)$
        * `**kwargs`: `softening_length` -->

#### storing_method
- `default`
    * Store solutions directly into memory
- `flush`
    * Flush intermediate results into a csv file to reduce memory pressure.
- `disabled`
    * To not store any result.

## Built-in systems
Some systems are available by default and can be loaded readily.
| System | Description |
|:-------|:------------| 
| circular_binary_orbit | A circular orbit formed by two stars |
| eccentric_binary_orbit | An eccentric orbit formed by two stars |
| 3d_helix | An upward helix consists of three stars |
| sun_earth_moon | The Sun, Earth, and Moon system |
| figure-8 | A "figure-8" orbit involving three stars  |
| pyth-3-body | Three stars arranged in a triangle with length ratios of 3, 4, and 5. It is a highly chaotic orbit with close encounters that can be used to test the difference between fixed and variable step size integrators. |
| solar_system | Solar System with the Sun and the planets |
| solar_system_plus | solar_system with the inclusion of Pluto, Ceres, and Vesta  |

## Integrators 
### Simple methods
Below are four simple fixed step size methods to simulate the system with a given step size $\text{d}t$.
| Simple methods |
|:-----------|
| Euler |
| Euler Cromer |
| Fourth Order Runge-Kutta (RK4) |
| Leapfrog |

### Embedded Runge-Kutta methods
Embedded RK methods are adaptive methods that decides the step size automatically based on the estimated error.
It can resolve close encounters but fail to conserve energy over long time scele.

| Embdedded Runge-Kutta methods | Recommended tolerance* |
|:-----------|:-------------|
| Runge–Kutta–Fehlberg 4(5) | $10^{-8}$ to $10^{-14}$ |
| Dormand–Prince method (DOPRI) 5(4) | $10^{-8}$ to $10^{-14}$ |
| Verner's method (DVERK) 6(5) | $10^{-8}$ to $10^{-14}$ |
| Runge–Kutta–Fehlberg 7(8) | $10^{-4}$ to $10^{-8}$ |

*For reference only

### IAS15
IAS15 (Implicit integrator with Adaptive time Stepping, 15th order) is a highly optimized integrator with extremely high accuracy. It is the default method for this project.

The recommended tolerance* is $10^{-9}$. Since the integrator is 15th order, changing the tolerance
results in little improvement in performance, but a huge penalty in accuracy. Therefore, it is not
recommended to change this tolerance.

*For reference only

### WHFast
WHFast is a second order symplectic method with fixed step size, which conserves energy over long integration period. This integrator cannot resolve close encounter.

#### `**kwargs` for WHFast
| Argument               | Description                                                  | Default Value |
|------------------------|--------------------------------------------------------------|---------------|
| `whfast_kepler_tol`           | Tolerance in solving the Kepler's equation                   | $10^{-12}$    |
| `whfast_kepler_max_iter`      | Maximum number of iterations in solving Kepler's equation    | 500        |
| `whfast_kepler_auto_remove`   | Integer flag to indicate whether to remove objects that failed to converge in Kepler's equation | False |
| `whfast_kepler_auto_remove_tol` | Tolerance for removing objects that failed to converge in Kepler's equation  | $10^{-8}$ |

> [!WARNING]\
> When using WHFast, the order of adding objects matters. Since WHFast use Jacobi coordinate, we must add the inner object first, followed by outer objects relative to the central star. For convenience, you may also add the objects in any order, then call `system.sort_by_distance(primary_object_name)` or `system.sort_by_distance(primary_object_index)`

## Saving the results
If you save the results, the data will be saved in the default unit (solar masses, AU and days), and follow this format:
```
time, dt, total energy, x1, y1, z1, ... vx1, vy1, vz1, ...
```
