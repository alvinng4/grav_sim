# Step 5: Adaptive time-stepping

In this step, we will implement an adaptive time-stepping integrator 
Runge-Kutta-Fehlberg (RKF) method, which belongs to the same family 
to the RK4 method we implemented in the previous step.

## RKF4(5) method

To achieve adaptive time-stepping, we will need to use two RK methods
of different orders to estimate the error. This error can then give
us a reliable estimation of the new time step. For RKF4(5), we use a fifth-order
method 

$$
    \tilde{x}_{n + 1} 
    = x_n
    + \left( \frac{16}{135} k_1 
    + \frac{6656}{12825} k_3
    + \frac{28561}{56430} k_4
    - \frac{9}{50} k_5
    + \frac{2}{55} k_6 \right) \Delta t,
$$

to estimate the local error for the fourth-order method that is used
for the actual update

$$
    x_{n + 1} 
    = x_n
    + \left( \frac{25}{216} k_1 
    + \frac{1408}{2565} k_3
    + \frac{2197}{4104} k_4
    - \frac{1}{5} k_5 \right) \Delta t.
$$

Since the $k$ between both methods mostly overlaps, we 
are able to obtain a error estimation with a very small additional cost.
They are given as follows:

$$
    \begin{aligned}
        k_1 &= f(t_n, x_n), \\
        k_2 &= f\left(t_n + \frac{1}{4} \Delta t, x_n + \frac{1}{4} k_1 \Delta t\right), \\
        k_3 &= f\left(t_n + \frac{3}{8} \Delta t, x_n + \left(\frac{3}{32} k_1 + \frac{9}{32} k_2 \right) \Delta t \right), \\
        k_4 &= f\left(t_n + \frac{12}{13} \Delta t, x_n + \left(\frac{1932}{2197} k_1 - \frac{7200}{2197} k_2 + \frac{7296}{2197} k_3\right) \Delta t\right), \\
        k_5 &= f\left(t_n + \Delta t, x_n + \left(\frac{439}{216} k_1 - 8 k_2 + \frac{3680}{513} k_3 - \frac{845}{4104} k_4\right) \Delta t\right), \\
        k_6 &= f\left(t_n + \frac{1}{2} \Delta t, x_n + \left(- \frac{8}{27} k_1 + 2 k_2 - \frac{3544}{2565} k_3 + \frac{1859}{4104} k_4 - \frac{11}{40} k_5\right) \Delta t\right).
    \end{aligned}
$$

This gives the following code:

```python
    # RKF4(5) coefficients
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
```

## Error estimation

First, let us define a tolerance parameter $\varepsilon$. This is a user-defined
parameter to control the time step (Smaller $\varepsilon$ means smaller time step).
Now, we define a $\Delta x'$ (we call it `error_estimation_delta` in our code) as the difference between the two RK methods

$$
    \Delta x' = \Delta x - \Delta \tilde{x}.
$$

Then, we compute $s$ (we call it `tolerance_scale` in our code) as follows:

$$
    s = \varepsilon + \varepsilon \times \max{(|x_n|, |x_{n +1}|)},
$$

where the maximum is taken element-wise (i.e. $s$ will have the same shape as $x$ which is `(N, 3)`).
Finally, we compute the error by taking the "norm":

$$
    \text{error} = \sqrt{\overline{\sum_i \left( \frac{\Delta x'}{s_i} \right)^2}}.
$$

The bar over the sum means that we take the average over all elements we summed over.
This gives us the code below. The denominator in the final line is `system.num_particles * 3.0 * 2.0` because we have $N$ particles, each with 3 dimensions, and we have two arrays $\mathbf{r}$ and $\mathbf{v}$.

```python
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
```

The new step is accepted if the error is less than or equal to 1.0.
Otherwise, we will reject the step and try again with a smaller time step.

## Time step estimation

With the error, we can now estimate the new time step.
We will use the following formula:

$$
    \Delta t_{n + 1} = 0.38^{1 / (1 + q)} \times \Delta t_n
    \times \text{error}^{-1 / (1 + q)},
$$

where $q$ is the lowest power of the two methods (4 in our case).
We also need to make sure that the new time step is not too small or too large.
We can do this by setting

* Safety factor `safety_fac_max = 6.0` and `safety_fac_min = 0.33` so that 
    $f_{\text{max}} \Delta t_{n + 1} \leq \Delta t_{n + 1} \leq f_{\text{min}} \Delta t_{n + 1}$.

* Lower bound of $\Delta t$ such that $\Delta t_{n + 1} \geq 10^{-12} (t_f - t_0)$.
* Lower bound of error such that $\text{error} \geq 10^{-12}$.

We have the following code:

```python
# Safety factors for step-size control
safety_fac_max = 6.0
safety_fac_min = 0.33
safety_fac = math.pow(0.38, 1.0 / (1.0 + float(min_power)))

...

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
```

## Initial dt
The final thing we need is to set the initial time step. There are methods to 
estimate the initial time step automatically. However, to keep it simple, we will
simply set it manually. Personally, I found it the best to just run
the simulation for a short time and see how the time step evolves. This could
even be more accurate than using an automatic method.

## Putting it all together
Now, we can put everything together. The final code is given below.
Note that we also store the `sol_dt` to keep track of the time step.

??? note "Code (Click to expand)"
    ```python
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
    ```

## Simulation results
Let us try to run the code again for the solar system. We use a tolerance = $10^{-8}$
and initial time step = 1.0 days. The plots are shown below, with
$\Delta t$ fluctuated quickly between $\sim 1.4 - 2.75$ days. This allows a more
flexible time stepping to reduce the computation cost.

![Relative energy error](../figures/step5/solar_system_200_yrs_rel_energy_error.png)
![dt](../figures/step5/solar_system_200_yrs_dt.png)

## Pythagorean Three-Body Problem

As solar system is a mostly stable system, we may not be able to see the benefits
of adaptive time-stepping. Here, we try the Pythagorean three-body problem, which is
a extremely chaotic system with close encounters. Below is a illustration of the
initial condition. We have three particles at rest at the vertices of a right-angled triangle
with length ratio 3:4:5. The mass of the particles are $3.0 / G, 4.0 / G$ and $5.0 / G$ respectively.

![Pythagorean three-body problem](../figures/step5/pyth_3_body_initial.png)

In Python, we have

```python
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
```

Also, we set the following simulation parameters:

```python
TF = 70.0  # 70 days
TOLERANCE = 1e-13
OUTPUT_INTERVAL = 0.001  # 0.001 day
INITIAL_DT = 0.01  # Initial time step in days
LABELS = [None, None, None]
COLORS = [None, None, None]
LEGEND = False
```

The simulation result is shown below:

![Trajectory](../figures/step5/pyth_3_body_trajectory_rkf45.png)
![Relative energy error](../figures/step5/pyth_3_body_rel_energy_error_rkf45.png)
![dt](../figures/step5/pyth_3_body_dt_rkf45.png)

As shown from the plots, $\Delta t$ fluctuated greatly between $10^{-8} - 10^{-2}$ days!
This is because of the close encounters between the particles. A video of the evolution
in real time is available on youtube:

<iframe width="560" height="315" src="https://www.youtube.com/embed/UfQAqfge_V4?si=PBgYRa_fAW6W-7ke" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

Can we simulate this with RK4? Yes, but not with the time step we used before.
Because the smallest time step for RKF4(5) is $10^{-8}$ days, we will also need to set the time
step for RK4 to $10^{-8}$ days **for the whole simulation**, which is very inefficient!
I have tested it using our `grav_sim` package written in C. The largest $\Delta t$ 
we can use is $\sim 2 \times 10^{-8}$ days, and the simulation took about 8 minutes,
while the RKF45 finished within seconds! In Python it would probably takes hours to run
the RK4 simulation.

## Summary

In this step, we have implemented the RKF4(5) method with adaptive time-stepping.
It is very efficient and allow us to save computational time, especially for
chaotic systems or systems with close encounters. In the next step, we will
first see how to add new particles to our system using orbital elements.
Then, I will show you how to make animations for the N-body simulations.

## Full script
The full script is available at `12_steps_to_n_body_simulation/python/step5.py`,
or https://github.com/alvinng4/grav_sim/blob/main/12_steps_to_n_body_simulation/python/step5.py

??? note "Code (Click to expand)"
    ```python linenums="1"
    --8<-- "12_steps_to_n_body_simulation/python/step5.py"
    ```