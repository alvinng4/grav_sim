# Step 3: Your first N-body program

In this step, we will write our first N-body program! 
We have set up the initial conditions in step 1 and 
the acceleration function in step 2. Now, we will need to
solve the equations of motion for the particles.
For Newtonian machanics, we have the following coupled equations:

$$
    \frac{\mathrm{d}\textbf{r}}{\mathrm{d}t} 
    = \textbf{v}, \quad
    \frac{\mathrm{d}\textbf{v}}{\mathrm{d}t} 
    = \textbf{a}(\textbf{r}) \left( = \frac{\textbf{F}(\textbf{r})}{m} \right).
$$

They are ordinary differential equations (ODEs), and to solve them, we will
need a ODE solver. 

## Euler method
The simplest ODE solver is the Euler method, where the update
is approximated as

$$
    \Delta \mathbf{r} = \mathbf{v} \Delta t, \quad 
    \Delta \mathbf{v} = \mathbf{a}(\mathbf{r}) \Delta t,
$$

where $\Delta t$ is the time step.
By taylor expansion, we can see that this is a first-order approximation,

$$
    x(t + \Delta t) = x(t) + x'(t) \Delta t + \mathcal{O}(\Delta t^2),
$$

where $\mathcal{O}(\Delta t^2)$ is the local truncation error. 
Because the total number of steps scales with $1/\Delta t$, the global error is 

$$
    \text{global error} =
    \text{number of steps} \times \text{error per step}
    \propto \frac{1}{\Delta t} \mathcal{O}(\Delta t^2) = \mathcal{O}(\Delta t).
$$

If you are not familiar with the big-O notation, you can think of it
as the error is bounded by a polynomial in the order of $\Delta t$ if $\Delta t \to 0$.

Implementing the Euler integrator is very easy. With

$$
    \mathbf{r}_{n + 1} = \mathbf{r}_n + \mathbf{v}_n \Delta t, \quad
    \mathbf{v}_{n + 1} = \mathbf{v}_n + \mathbf{a}(\mathbf{r}_{n}) \Delta t,
$$

we have the following code:

```python
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
```

Now, we will build all the components we need for our N-body program.

## Solution output

We will need a way to store the solution. A naive way is to store the system
at every time step or every few time steps. However, this is a terrible idea
because the output size will depends on your choice of time step. A better way
is to store the solution at regular intervals.
In our simulation, we will use a output interval of 0.1 years. For a simulation of
200 years, we will have 2000 time steps.

!!! Tip
    For solar system, we only have 9 particles and it takes very little memory to store. 
    So, you don't need to worry too much about the solution size. Just be careful
    and don't set the output interval too small.

Before the simulation, we will need to set up the output array and store the initial conditions.
```python
OUTPUT_INTERVAL = 0.1 * 365.24  # years to days

def main() -> None:
    ...
    # Solution array
    sol_size = int(TF // OUTPUT_INTERVAL + 2)  # +2 for initial and final time
    sol_x = np.zeros((sol_size, system.num_particles, 3))
    sol_v = np.zeros((sol_size, system.num_particles, 3))
    sol_t = np.zeros(sol_size)
    sol_x[0] = system.x
    sol_v[0] = system.v
    sol_t[0] = 0.0
    output_count = 1
```

Also, we need to calculate the output time.
```python
def main() -> None:
    ...
    for i in range(NUM_STEPS):
        ...
        if current_time >= next_output_time:
            sol_x[output_count] = system.x
            sol_v[output_count] = system.v
            sol_t[output_count] = current_time

            output_count += 1
            next_output_time = output_count * OUTPUT_INTERVAL
```

Finally, we resize the arrays to the actual size.
```python
def main() -> None:
    ...
    sol_x = sol_x[:output_count]
    sol_v = sol_v[:output_count]
    sol_t = sol_t[:output_count]
```

## Putting it all together

Let's put everything together. We first need to setup the simulation parameters.
```python
# Default units is AU, days, and M_sun
TF = 200.0 * 365.24  # years to days
DT = 1.0 
OUTPUT_INTERVAL = 0.1 * 365.24  # years to days
NUM_STEPS = int(TF / DT)
```

Before running the simulation, it is a good idea to print
the simulation information.
```python
def print_simulation_info(system: System, sol_size: int) -> None:
    print("----------------------------------------------------------")
    print("Simulation Info:")
    print(f"num_particles: {system.num_particles}")
    print(f"G: {system.G}")
    print(f"tf: {TF} days (Actual tf = dt * num_steps = {DT * NUM_STEPS} days)")
    print(f"dt: {DT} days")
    print(f"Num_steps: {NUM_STEPS}")
    print()
    print(f"Output interval: {OUTPUT_INTERVAL} days")
    print(f"Estimated solution size: {sol_size}")
    print("----------------------------------------------------------")
```

Finally, we have the main function that runs the main simulation loop.
I have added a timer to measure the runtime, and a print statement
to show the simulation progress everytime we save a solution. The `\r` at the end
of the print statement should overwrite the previous line, but you can remove
the print statement if it is spamming your terminal.
```python hl_lines="34"
def main() -> None:
    # Get initial conditions
    system = get_initial_conditions()

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
    print_simulation_info(system, sol_size)
    next_output_time = output_count * OUTPUT_INTERVAL
    start = timeit.default_timer()
    for i in range(NUM_STEPS):
        euler(a, system, DT)

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
```

You should see the following output:

```
----------------------------------------------------------
Simulation Info:
num_particles: 9
G: 0.00029591220828411956
tf: 73048.0 days (Actual tf = dt * num_steps = 73048.0 days)
dt: 1.0 days
Num_steps: 73048

Output interval: 36.524 days
Estimated solution size: 2001
----------------------------------------------------------
Current time: 73012.00 days
Done! Runtime: 1.1 seconds, Solution size: 2000
```

## Plotting the trajectory

To visualize the trajectory, we have the following function. Notice that we 
have two `ax.plot` calls. The first one is to plot the trajectory, and the second one
is to plot the final position as a circle marker.

```python
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
    legend : bool
        Whether to show the legend.
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
```

Then, we add the function call to the end of the main function.
```python
plot_trajectory(
    sol_x=sol_x,
    labels=LABELS,
    colors=COLORS,
    legend=LEGEND,
)
```

You should see the following plot:

![Trajectory](../figures/step3/trajectory.png)

Congrats! :partying_face: You have just written your first N-body simulation program!
However, we can see that the results are not very accurate, especially for those inner planets.
(Mercury has drifted to a orbit beyond Saturn within 200 years!)

In the next step, we will implement some higher-order integrators to improve the accuracy.

## Full script
The full script is available at `6_steps_to_n_body_simulation/python/step3.py`,
or https://github.com/alvinng4/grav_sim/blob/main/6_steps_to_n_body_simulation/python/step3.py

??? note "Code (Click to expand)"
    ```python linenums="1"
    --8<-- "6_steps_to_n_body_simulation/python/step3.py"
    ```