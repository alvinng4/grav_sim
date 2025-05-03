In our project, we support two output formats: CSV and HDF5.

## Regular simulations

### CSV Format

For the CSV format, we output snapshots with file names in the format of `snapshot_%05d.csv`.
The file begins with metadata. For example:

```
# num_particles: 3
# G: 0.00029591220828411951
# time: 0
# dt: 0
```

Then, the data simply follows the header:

```
particle_id,m,x,y,z,vx,vy,vz
```

### HDF5 Format

For the HDF5 format, we output snapshots with file names in the format of `snapshot_%05d.hdf5`.
It has two groups: `Header` and `PartType0`. Below is a tree structure of the HDF5 file:

```
snapshot_%05d.hdf5
├── Header
│   ├── Attribute: NumFilesPerSnapshot
│   ├── Attribute: NumPart_ThisFile
│   ├── Attribute: NumPart_Total
│   └── Attribute: Time
└── PartType0
    ├── Dataset: ParticleIDs
    ├── Dataset: Masses
    ├── Dataset: Coordinates
    └── Dataset: Velocities
```

## Cosmology simulations

In cosmology simulations, we only support the HDF5 format.

### HDF5 Format

For the HDF5 format, we output snapshots with file names in the format of `snapshot_%05d.hdf5`.
It has three groups: `Header`, `PartType0` and `Units`. Below is a tree structure of the HDF5 file:

```
HDF5 File
├── Header
│   ├── Attribute: NumFilesPerSnapshot
│   ├── Attribute: NumPart_ThisFile
│   ├── Attribute: NumPart_Total
│   ├── Attribute: Time
│   ├── Attribute: Redshift
│   ├── Attribute: Omega0
│   └── Attribute: OmegaLambda
├── PartType0
│   ├── Dataset: ParticleIDs
│   ├── Dataset: Masses (float or double)
│   ├── Dataset: Coordinates (float or double)
│   └── Dataset: Velocities (float or double)
└── Units
    ├── Attribute: Unit current in cgs (U_I)
    ├── Attribute: Unit length in cgs (U_L)
    ├── Attribute: Unit mass in cgs (U_M)
    ├── Attribute: Unit time in cgs (U_t)
    └── Attribute: Unit temperature in cgs (U_T)
```