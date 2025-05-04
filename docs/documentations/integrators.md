## Simple methods
Below are four simple fixed step size methods to simulate the system with a given step size $\Delta t$.

| Simple methods |
|:-----------|
| Euler |
| Euler Cromer |
| Fourth Order Runge-Kutta (RK4) |
| Leapfrog |

## Embedded Runge-Kutta methods
Embedded RK methods are adaptive methods that decides the step size automatically based on the estimated error.
They can resolve close encounters but fail to conserve energy over long time scele.

| Embdedded Runge-Kutta methods | Recommended tolerance (for reference only) |
|:-----------|:-------------|
| Runge–Kutta–Fehlberg 4(5) | $10^{-8}$ to $10^{-14}$ |
| Dormand–Prince method (DOPRI) 5(4) | $10^{-8}$ to $10^{-14}$ |
| Verner's method (DVERK) 6(5) | $10^{-8}$ to $10^{-14}$ |
| Runge–Kutta–Fehlberg 7(8) | $10^{-4}$ to $10^{-8}$ |

## IAS15
IAS15 (Implicit integrator with Adaptive time Stepping, 15th order)
is a highly optimized integrator with extremely high accuracy.
The recommended tolerance is $10^{-9}$. Since the integrator is 15th order, changing the tolerance
results in little improvement in performance, but a huge penalty in accuracy. Therefore, it is not
recommended to change this tolerance.

## WHFast
WHFast is a second order symplectic method with fixed step size, which conserves 
energy over long integration period. This integrator is suitable for stable central
mass systems like the solar system.