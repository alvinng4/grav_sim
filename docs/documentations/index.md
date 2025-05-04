## Table of Contents

### Documentations
- [Integrators](integrators.md)
- [Massless acceleration](massless_acceleration.md)
- [Barnes-Hut algorithm](barnes_hut.md)
- [Particle-Mesh algorithm](particle_mesh.md)
- [Comoving coordinates](comoving_coordinates.md)
- [Output formats](output_formats.md)
- [Force softening](force_softening.md)
- [Reducing round off error](reducing_round_off_error.md)
- [Why C?](why_c.md)

### Python API
- [GravitySimulatorAPI](PythonAPI/GravitySimulatorAPI.md)
- [Parameters](PythonAPI/parameters.md)
- [Plotting](PythonAPI/plotting.md)
- [System](PythonAPI/System.md)
- [Simulator](PythonAPI/Simulator.md)

### C API
- [Grav sim](CAPI/grav_sim.md)




## Project File Structure

```
grav_sim/
├── docs/                           # Documentation
├── examples/                       # Example scripts
├── grav_sim/                       # Python wrapper
├── overrides/                      # Mkdocs overrides
├── pcg/                            # PCG random number generator
├── src/                            # C source code
├── 5_steps_to_n_body_simulation/   # 5 steps to n-body simulation
├── .gitignore
├── CMakeLists.txt
├── FindFFTW3.cmake
├── LICENSE
├── MANIFEST.in
├── mkdocs.yml
├── pyproject.toml
├── README.md            
├── requirements.txt             
├── setup.py      
└── .github/   
```