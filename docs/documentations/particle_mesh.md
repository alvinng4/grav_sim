??? Note "Source code (Click to expand)"
    ```c linenums="1"
    --8<-- "src/acceleration_PM.c"
    ```

For cosmological simulations, we are usually interested in a large domain with
periodic boundary conditions. Interestingly, grid-based algorithm
like the Particle Mesh method naturally provides periodic boundary conditions
by solving the Poisson's equation in the Fourier Space,

\begin{equation}
    \label{eqn:poisson}
    \nabla^2 \Phi(x) = 4 \pi G \rho(\textbf{x})
    \implies - \textbf{k}^2 \Phi(\textbf{k}) = 4 \pi G \rho(\textbf{k}).
\end{equation}

The overall time complexity is \(\mathcal{O}(N + N_{\textnormal{cell} } \log N_{\textnormal{cell} })\),
with \(\mathcal{O}(N)\) for interpolation and \(\mathcal{O}(N_{\textnormal{cell} } \log N_{\textnormal{cell} })\) for FFT.
Therefore, this method is very fast. However, it is not accurate for short range force where
the length scale is smaller than the grid resolution. Thus, it is not suitable for inhomogeneous
distributions with high density regions.

In this section, we would focus on the details in implementing the Particle Mesh algorithm.
There are three steps in the algorithm: 

1. Estimating the density of the underlying grid,
2. Solving the Poisson's equation in fourier space, and
3. Computing the acceleration by taking a gradient of the potential in real space
and interpolating the acceleration to the particle.


## Density estimation
In a sense, N-body simulation is a Monte Carlo technique for large scale simulations.
Each particle is a sampling drawn from a distribution at \(t = 0\), where all \(N\)
particles together models the underlying continuous distribution[@dehnen_n-body_2011].
Therefore, when estimating the density of the underlying grid, we should not think of the particles as a
localized point with extreme high density (e.g. an enormous star), but a cloud of mass that
is represented as a particle due to limited computational capability.

In figure 1, we shows the comparison between the Nearest Grid Point (NGP)
scheme versus the Cloud-In-Cell (CIC) scheme. The latter provides a much
smoother estimation of density. In 3D, a particle at \(\textbf{x}_p = (x_p, y_p, z_p)\) could contribute density to its
parent cell \(\textbf{x}_c = (x_c, y_c, z_c)\) and seven neighboring cells. Define the weightings

\begin{equation}
    %\label{}
    \textbf{d} = \frac{\textbf{x}_p - \textbf{x}_c}{l}, \quad
    \textbf{t} = 1 - \textbf{d},
\end{equation}

where \(l\) is the cell length. We have the density contribution due to a particle[@how_to_write_pm]

\begin{align}
    \rho_{i, j, k} &= \rho_p t_x t_y t_z,       &\rho_{i + 1, j + 1, k} &= \rho_p d_x d_y t_z, \\
    \rho_{i + 1, j, k} &= \rho_p d_x t_y t_z,   &\rho_{i + 1, j, k + 1} &= \rho_p d_x t_y d_z, \nonumber \\
    \rho_{i, j + 1, k} &= \rho_p t_x d_y t_z,   &\rho_{i, j + 1, k + 1} &= \rho_p t_x d_y d_z, \nonumber \\
    \rho_{i, j, k + 1} &= \rho_p t_x t_y d_z,   &\rho_{i + 1, j + 1, k + 1} &= \rho_p d_x d_y d_z, \nonumber
\end{align}

where \(\rho_p = m_p / V_{\textnormal{unit cell}}\).  

<figure style="text-align: center;">
  <img src="../../../examples/media/cloud_in_cell.png" alt="Cloud-in-cell" width="600" style="display: block; margin: auto;" />
  <figcaption>Figure 1: Two different schemes in estimating density in a discretized particle mesh.</figcaption>
</figure>

## Solving the Poisson's equation
Now, we have a grid with assigned density value on each grid point. To transform it into
Fourier Space, we simply do a discrete fourier transform using FFT libraries
(e.g. FFTW in C or `numpy.fft` in Python).
Then, we compute the wave numbers by

\begin{equation}
    %\label{}
    \textbf{k} = \left( \frac{2 \pi n_x}{L_x}, \frac{2 \pi n_y}{L_y}, \frac{2 \pi n_z}{L_z} \right),
\end{equation}

where \(\textbf{L} = (L_x, L_y, L_z)\) is the box length, and
\(\textbf{n} = (n_x, n_y, n_z) \in [0, N)\)
is the grid indices (one may need to wrap the indices above $\frac{N}{2}$ to $-\frac{N}{2}$ so that 
\((n_x, n_y, n_z) \in [- \frac{N}{2}, \frac{N}{2})\)).
Then, we simply compute \(\textbf{k}^2\) and obtain the potential in Fourier Space

\begin{equation}
    %\label{}
    \Phi(\textbf{k}) = -\frac{4 \pi G}{\textbf{k}^2} \rho(\textbf{k}).
\end{equation} 

## Interpolating the acceleration
Now, we do an inverse discrete fourier transform to obtain the potential on the grid points.
Then, the acceleration can be obtained by computing \(\textbf{a} = - \nabla \Phi\) by a central
finite difference scheme with fourth order accuracy:

\begin{equation}
    %\label{}
    f'(x) \approx \frac{1}{l^3} \left[ \frac{1}{12} f(x - 2l) - \frac{2}{3} f(x - l) + \frac{2}{3} f(x + l) - \frac{1}{12} f(x + 2l) \right] + \mathcal{O}(l^4),
\end{equation} 

Then, we perform a linear interpolation to obtain the acceleration of the particles.
Reusing the notations from the cloud-in-cell scheme, we have

\begin{align}
    %\label{eqn:eqlabel}
    \textbf{a}_p 
    &= \textbf{a}_{i, j, k} t_x t_y t_z
    + \textbf{a}_{i + 1, j, k} d_x t_y t_z
    + \textbf{a}_{i, j + 1, k} t_x d_y t_z
    + \textbf{a}_{i, j, k + 1} t_x t_y d_z \nonumber \\
    &\quad\,+
    \textbf{a}_{i + 1, j + 1, k} d_x d_y t_z
    + \textbf{a}_{i + 1, j, k + 1} d_x t_y d_z
    + \textbf{a}_{i, j + 1, k + 1} t_x d_y d_z
    + \textbf{a}_{i + 1, j + 1, k + 1} d_x d_y d_z. 
\end{align}
