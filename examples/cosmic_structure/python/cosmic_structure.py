from pathlib import Path

import h5py
import numpy as np
from grav_sim import GravitySimulatorAPI

INIT_CONDITION_FILE = Path(__file__).parent.parent / "ics_swift.hdf5"
A_FINAL = 1.0  # Scale factor at the end of simulation
NUM_STEPS = 1000


def main() -> None:
    gs = GravitySimulatorAPI()

    system = gs.get_new_cosmological_system()
    pm_grid_size = read_init_condition(system)

    _, _, output_param, settings = gs.get_new_parameters()
    output_param.method = "hdf5"
    output_param.output_dir = Path(__file__).parent / "snapshots/"
    output_param.output_interval = (A_FINAL - system.scale_factor) / 100
    output_param.output_initial = True
    output_param.coordinate_output_dtype = "float"
    output_param.velocity_output_dtype = "float"
    output_param.mass_output_dtype = "float"

    settings.verbose = "normal"

    gs.launch_cosmological_simulation(
        system,
        output_param,
        settings,
        A_FINAL,
        NUM_STEPS,
        pm_grid_size,
    )


def read_init_condition(system) -> int:
    with h5py.File(INIT_CONDITION_FILE, "r") as f:
        header = f["/Header"]
        units = f["/Units"]
        part_type_1 = f["/PartType1"]

        # Read header attributes
        num_particles = int(header.attrs["NumPart_ThisFile"][1])
        box_size = float(header.attrs["BoxSize"])
        pm_grid_size = int(header.attrs["suggested_pmgrid"])
        omega_m = float(header.attrs["Omega0"])
        omega_lambda = float(header.attrs["OmegaLambda"])
        redshift = float(header.attrs["Redshift"])
        h_param = float(header.attrs["HubbleParam"])

        system.box_width = box_size / 2.0
        system.box_center = np.array(
            [system.box_width, system.box_width, system.box_width], dtype=np.float64
        )
        system.scale_factor = 1.0 / (redshift + 1.0)
        system.omega_m = omega_m
        system.omega_lambda = omega_lambda
        system.h = h_param

        # Read unit attributes
        system.unit_length_in_cgs = float(units.attrs["Unit length in cgs (U_L)"])
        system.unit_mass_in_cgs = float(units.attrs["Unit mass in cgs (U_M)"])
        system.unit_time_in_cgs = float(units.attrs["Unit time in cgs (U_t)"])

        # Read particle datasets
        system.add(
            part_type_1["Coordinates"][:],
            part_type_1["Velocities"][:],
            part_type_1["Masses"][:],
        )

        return pm_grid_size


if __name__ == "__main__":
    main()
