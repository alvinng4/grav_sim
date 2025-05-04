# Asteroid belt animation

??? Note "Source code (Click to expand)"
    ```python linenums="1" title="examples/asteroid_belt_animation/python/asteroid_belt.py"
    --8<-- "examples/asteroid_belt_animation/python/asteroid_belt.py"
    ```

In this project, we will simulate the asteroid belt in the solar system by adding
50000 particles using the `system.add_keplerian` method. We also provide 3 options:

* `0`: asteroids only
* `1`: Add a star near the asteroid belt
* `2`: Same as `1` but different position

Because the simulation only consider short time scale, 

* RK4 is chosen for higher accuracy
* Gravity due to the asteroids are ignored

To draw animation, we use matplotlib to draw the frames, then use PIL to convert
the frames to a gif file.

## Gallery

<iframe width="560" height="315" src="https://www.youtube.com/embed/C45ceYja0jE?si=WrTNKP3ht_J8wxOR" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

<iframe width="560" height="315" src="https://www.youtube.com/embed/eg7plHjP1eg?si=4dlGoy_msuW3o810" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

<iframe width="560" height="315" src="https://www.youtube.com/embed/HMv7OwqAmBY?si=0HTV_61L_3r4CEun" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>