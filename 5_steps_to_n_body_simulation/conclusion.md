# Conclusion and Final Project

Congratulations on finishing the tutorial!
Now that you have obtained the basic knowledge in N-body simulation,
it is time for you to make your own project! There are so many things
you can do with N-body simulation. Below are some possible ideas,
but you can also come up with your own project.

- **Game**: You can make a simple game using N-body simulation. (In fact,
    I made a game called [OrbitSim](https://github.com/alvinng4/OrbitSim)
    when I first started learning Python and N-body simulation.)

    <img src="../../examples/media/OrbitSim.png" alt="OrbitSim" width="400"/>

- **Rewrite the code**: You can rewrite the code in C / C++ or your
    favorite language and compare the performance. You can also try
    to use Numba or Cython to speed up your code in Python. Maybe
    write your own N-body simulation library?

- **Animation**: You can make some nice animations of
    N-body simulation using Matplotlib in Python. Simply draw the
    frames, save them as images and then combine them into a video.

- **Large-scale simulation**: So far, we have only focused on systems
    with a few objects. What about large-scale systems with thousands
    or millions of objects? Turns out it is not so easy because the 
    computation of gravity scales as $\mathcal{O}(N^2)$. Have a look
    at Barnes-Hut algorithm to see how to speed up the simulation.
    In fact, we have documentations about it on this website as well.
    Try to implement it in some low-level language like C or C++!

    <img src="../../examples/media/galaxy_collision.png" alt="Galaxy Collision" width="400"/>

- **Reproducing observations**: If you are not interested in
    writing low-level simulation code, you can still do some interesting projects
    using our `grav_sim` project written in C. One example is reproducing the
    Kirkwood gaps in our solar system.

    <img src="../../examples/media/Kirkwood_gap_visualization.png" alt="Kirkwood Gaps" width="400"/>

Good luck and have fun!
    
