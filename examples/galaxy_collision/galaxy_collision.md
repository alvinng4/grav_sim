??? Note "Source code in C (Click to expand)"
    ```python linenums="1" title="examples/galaxy_collision/c/galaxy_collision.c"
    --8<-- "examples/galaxy_collision/c/galaxy_collision.c"
    ```

??? Note "Source code in Python (Click to expand)"
    ```python linenums="1" title="examples/galaxy_collision/python/galaxy_collision.py"
    --8<-- "examples/galaxy_collision/python/galaxy_collision.py"
    ```

Related topics: 

* [Barnes-Hut algorithm](../../docs/documentations/barnes_hut.md)

In this example, we will simulate the collision of two galaxies with 60000 particles using the initial
conditions from Gadget-2[@gadget2]. The initial conditions is preprocessed into a
HDF5 file and is included in the repository.

To deal with the large number of particles, we used the Barnes-Hut algorithm to
compute the gravitational force. On Macbook Air M1, the simulation can be done
in 5 minutes for $\theta = 1.0$ and 30 minutes for $\theta = 0.5$.

## Gallery
The visualization is done with gadgetviewer, which is compatible with our
HDF5 output format.

<iframe width="560" height="315" src="https://www.youtube.com/embed/nXTUdjLXwtI?si=0MDsZFkRNru9G8l0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>