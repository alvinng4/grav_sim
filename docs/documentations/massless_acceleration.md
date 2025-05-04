??? Note "Source code (Click to expand)"
    ```c linenums="1"
    --8<-- "src/acceleration.c"
    ```

Let's say you have 10 regular particles and 10 massless particles
(i.e. particles with mass so small that is negligible).
One may compute the acceleration with the brute-force algorithm

$$
\mathbf{a}_i = \sum_{j \neq i} \frac{G m_j}{r_{ij}^2} \hat{\mathbf{r}}_{ij}
\quad \text{for } i = 1, \ldots, N
$$

But this is $\mathcal{O}(N^2)$, and we don't actually need to compute
the acceleration due to the massless particles. Therefore, a more
efficient way is to separate the calculations for massive and massless particles.

We first produce two list of indices, one for the massive particles and 
one for the massless particles. Then, we compute the acceleration for the massive particles
with the brute-force algorithm. For the massless particles, we compute the acceleration
of them due to the massive particles only. This gives a time complexity: $O(M^2 + MN)$,
where $M$ and $N$ are the number of massive and massless particles respectively.