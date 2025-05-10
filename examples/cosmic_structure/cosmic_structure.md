# Cosmic Structure Formation

??? Note "Source code in C (Click to expand)"
    ```python linenums="1" title="examples/cosmic_structure/c/cosmic_structure.c"
    --8<-- "examples/cosmic_structure/c/cosmic_structure.c"
    ```

??? Note "Source code in Python (Click to expand)"
    ```python linenums="1" title="examples/cosmic_structure/python/cosmic_structure.py"
    --8<-- "examples/cosmic_structure/python/cosmic_structure.py"
    ```

Related topics: 

* [Comoving coordinates](../../docs/documentations/comoving_coordinates.md)
* [Particle mesh algorithm](../../docs/documentations/particle_mesh.md)


In this example, we will simulate the formation of cosmological structures in a universe with $\Lambda$CDM model.

The initial conditions were generated with MUSIC v2[@MUSIC] at
\(z = 50\), with a resolution of \(128^3 \sim 2\) million particles in a 30 Mpc \(/ h\) box.
The MUSIC configuration file and initial conditions are provided in the repository, but with
a smaller particles count of \(64^3\) to reduce the repository size.

To run the simulation, we used the LeapFrog integrator in comoving coordinates,
with a particle mesh algorithm with mesh size = \(256^3 \sim 16\) million cells (following the suggestions
from the MUSIC initial conditions).

## Gallery
The simulation was done on Macbook Air M1. The simulation ran surprisingly fast and is
finished in less than 10 minutes. The visualization was done with gadgetviewer as it is 
directly compatible with our HDF5 output format.

<iframe width="560" height="315" src="https://www.youtube.com/embed/yof2x_0IeOA?si=8f9Ip5BYmSNwQuNJ" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

