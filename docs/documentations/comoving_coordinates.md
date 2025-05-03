To account for cosmological expansion, we need to use comoving coordinates.
We first review particle motion in comoving coordinates following LSSU[@LSSU].
In comoving coordinates, we have

\begin{equation}
    %\label{}
    \textbf{r} = a(t) \textbf{x},
\end{equation}

where \(\textbf{x}\) is a vector in comoving coordinates and \(a(t)\) is the scale factor. 
We may express the field equation as

\begin{equation}
    \label{eqn:poisson_peculiar}
    \nabla^2 \phi = 4 \pi G [\rho - \overline{\rho}] a^2,
\end{equation} 

where \(\phi\) is the peculiar potential, \(\rho\) and \(\overline{\rho}\)
are the density and average background density respectively. Note that
the gradient is taken with respect to \(\textbf{x}\). 
In addition, we have the proper velocity of a particle

\begin{equation}
    %\label{}
    \textbf{u} = a \dot{\textbf{x}} + \textbf{x} \dot{a},
\end{equation}

so that the Lagrangian for the particle motion is

\begin{equation}
    %\label{}
    \mathcal{L} = \frac{1}{2} m (a \dot{\textbf{x}} + \textbf{x} \dot{a})^2 - m \Phi.
\end{equation}

The canonical transformation

\begin{equation}
    %\label{}
    \mathcal{L} \to \mathcal{L} - \frac{\mathrm{d}\psi}{\mathrm{d}t},
    \quad \psi = \frac{1}{2} m a \dot{a} x^2, 
\end{equation}

reduces the Lagrangian to

\begin{equation}
    %\label{}
    \mathcal{L} = \frac{1}{2} m a^2 \dot{x}^2 - m \phi,
    \quad \phi = \Phi + \frac{1}{2} a \ddot{a} x^2.
\end{equation}

Now, we could define the canonical momentum

\begin{equation}
    \label{eqn:canonical_momentum}
    \textbf{p} = m a^2 \dot{x},
    \quad \frac{\mathrm{d}\textbf{p}}{\mathrm{d}t} = - m \nabla \phi.
\end{equation}

With the proper peculiar velocity

\begin{equation}
    %\label{}
    \textbf{v}_{\textnormal{pec} } = a \dot{\textbf{x}},
\end{equation}

we obtain from the canonical momentum equation,

\begin{equation}
    \label{eqn:peculiar_velocity_derivative}
    \frac{\mathrm{d}\textbf{v}_\textnormal{pec} }{\mathrm{d}t}
    + \textbf{v}_\textnormal{pec}  \frac{\dot{a}}{a}
    = - \frac{\nabla \phi}{a}.
\end{equation}

In our program, since matter is the
main component that contributes to the density perturbation, we rewrite the poisson equation
as

\begin{equation}
    %\label{}
    \nabla^2 \phi = \frac{4 \pi G}{a} [\rho_{m, \textnormal{comoving} } - \overline{\rho}_{m, \textnormal{comoving} }].
\end{equation}

where

\begin{equation}
    %\label{}
    \rho_{m, \textnormal{comoving} } - \overline{\rho}_{m, \textnormal{comoving} }
    = a^3 [\rho_{m} - \overline{\rho}_{m}].
\end{equation}

In fact, we can even drop the background density term since we only care about the
acceleration \(\mathbf{a} = - \nabla \phi\).
Also, instead of the canonical momentum, we define a conjugate momentum

\begin{equation}
    %\label{}
    \textbf{p}' = a^2 \dot{x} = a \textbf{v}_{\textnormal{pec} },
    \quad \frac{\mathrm{d}\textbf{p}'}{\mathrm{d}t} = - \nabla \phi.
\end{equation}

Following Gadget-4[@gadget4], we use \(\tau = \ln a\) as the integration variable.
We have

\begin{equation}
    \frac{\mathrm{d} a}{\mathrm{d} \tau}
    = \frac{\mathrm{d} a}{\mathrm{d} \ln a}
    = \left(\frac{\mathrm{d} \ln a}{\mathrm{d} a} \right)^{-1}
    = a,
\end{equation}

\begin{equation}
    \frac{\mathrm{d} \textbf{p}'}{\mathrm{d} \tau}
    = \frac{\mathrm{d} a}{\mathrm{d} \tau} \frac{\mathrm{d}t}{\mathrm{d}a} \frac{\mathrm{d} \textbf{p}'}{\mathrm{d} t}
    = - \nabla \phi \frac{a}{\dot{a}}
    = - \frac{\nabla \phi}{H(a)},
\end{equation}

where

\begin{equation}
    H(a) = H_0 \sqrt{\frac{\Omega_{m, 0}}{a^{3} } + \frac{1 - \Omega}{a^{2}} + \Omega_{\Lambda, 0}}\,\,,
\end{equation}

assuming \(\Omega_{\textnormal{radiation} } / a^4 \sim 0\). For \(\textbf{x}\),

\begin{equation}
    %\label{}
    \frac{\partial \textbf{x}}{\partial \tau}
    = \frac{\mathrm{d} a}{\mathrm{d} \tau} \frac{\mathrm{d}t}{\mathrm{d}a}  \dot{x}
    = \frac{a}{\dot{a}} \left( \frac{\textbf{p}'}{a^2} \right)
    = \frac{\textbf{p}'}{a^2 H(a)}.
\end{equation}

With these relations, we could now do time integration using \(\textbf{x}\) and \(\textbf{p}'\) in comoving coordinates. 