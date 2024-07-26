# Argon Simulation Project

This project simulates the behavior of argon atoms in a box using molecular dynamics techniques. It includes a Rust program for running the simulation and a Python script for visualizing the results.

## Components

1. `main.rs`: The main simulation program written in Rust.
2. `sim.py`: A Python script for visualizing the simulation results.

## Requirements

### For the Rust simulation:
- Rust programming language
- Required Rust crates (specified in `Cargo.toml`):
  - `rand`
  - `serde`
  - `indicatif`
  - `serde_json`
  

### For the Python visualization:
- Python 3.x
- Required Python libraries:
  - `numpy`
  - `matplotlib`
  - `moviepy`

## Usage

### Running the Simulation

1. Compile and run the Rust program with the following command:

```
cargo run -- <box_length> <num_atoms> <timestep> <total_steps> <snapshot_interval>
```

Replace the placeholders with appropriate values:
- `<box_length>`: Length of the simulation box (in angstroms)
- `<num_atoms>`: Number of argon atoms to simulate
- `<timestep>`: Simulation timestep (in femtoseconds)
- `<total_steps>`: Total number of simulation steps
- `<snapshot_interval>`: Interval at which to save snapshots of the system

Example:
```
cargo run -- 10.0 100 0.001 10000 100
```

This will run the simulation and generate a `simulation_data.json` file containing the trajectory data.

### Visualizing the Results

1. After running the simulation, use the Python script to visualize the results:

```
python sim.py
```

This will create a 3D animation of the argon atoms' movement and save it as `argon_simulation.mp4`.

## Features

- Implements the Lennard-Jones potential for argon atom interactions
- Uses the Verlet integration method for updating atom positions
- Applies periodic boundary conditions
- Implements the Berendsen thermostat for temperature control
- Provides a progress bar during the simulation
- Generates a JSON output file with simulation data
- Creates a 3D visualization and animation of the simulation results

## Notes

- The simulation uses reduced units for simplicity.
- The target temperature is set to 87.3 K (adjustable in the code).
- The visualization script loads the data from `simulation_data.json`, so make sure this file is in the same directory when running `sim.py`.

## Contributing

Feel free to fork this project and submit pull requests with improvements or additional features.

## License

This project is open-source and available under the GNU GENERAL PUBLIC LICENSE.