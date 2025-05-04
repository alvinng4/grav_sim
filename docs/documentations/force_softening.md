To avoid singularity for r = 0, we implemented a very simple force softening model:
```c
const double R_norm = sqrt(
    R[0] * R[0] + 
    R[1] * R[1] + 
    R[2] * R[2] +
    softening_length * softening_length
);
```            