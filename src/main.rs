use rand::Rng;
use serde::Serialize;
use std::fs::File;
use std::io::Write;
use std::env;
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Serialize)]
struct SimulationData {
    box_length: f64,
    num_atoms: usize,
    timestep: f64,
    total_steps: usize,
    snapshot_interval: usize,
    trajectory: Vec<Vec<[f64; 3]>>,
}

fn lj_potential(r: f64) -> f64 {
    let sigma = 1.0;
    let epsilon = 1.0;
    4.0 * epsilon * ((sigma / r).powi(12) - (sigma / r).powi(6))
}

fn main() {
    let args: Vec<String> = env::args().collect();
    
    if args.len() != 6 {
        eprintln!("Usage: {} <box_length> <num_atoms> <timestep> <total_steps> <snapshot_interval>", args[0]);
        std::process::exit(1);
    }

    let l: f64 = args[1].parse().expect("Invalid box length");
    let n: usize = args[2].parse().expect("Invalid number of atoms");
    let dt: f64 = args[3].parse().expect("Invalid timestep");
    let steps: usize = args[4].parse().expect("Invalid total steps");
    let snapshot_interval: usize = args[5].parse().expect("Invalid snapshot interval");

    let target_temperature: f64 = 87.3; // Target temperature
    let tau: f64 = 0.1; // Coupling constant for the Berendsen thermostat

    let mut rng = rand::thread_rng();
    let mut positions = vec![[0.0; 3]; n];
    let mut velocities = vec![[0.0; 3]; n];
    let mut positions_old = vec![[0.0; 3]; n];

    // Initialize positions randomly in the box
    for pos in positions.iter_mut() {
        for coord in pos.iter_mut() {
            *coord = rng.gen::<f64>() * l;
        }
    }

    // Initialize velocities from Maxwell-Boltzmann distribution
    let mass_argon: f64 = 39.95;
    let kb: f64 = 0.0083144621;
    let velocity_factor = (kb * target_temperature / mass_argon).sqrt();
    for vel in velocities.iter_mut() {
        for coord in vel.iter_mut() {
            *coord = rng.gen::<f64>() * velocity_factor;
        }
    }

    positions_old.copy_from_slice(&positions);

    let mut trajectory = Vec::new();

    // Set up the progress bar
    let pb = ProgressBar::new(steps as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
        .unwrap()
        .progress_chars("##-"));

    // Perform simulation
    for step in 0..steps {
        // Update progress bar
        pb.set_position(step as u64);

        // Calculate forces
        let mut forces = vec![[0.0; 3]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let mut r_ij = [0.0; 3];
                for k in 0..3 {
                    r_ij[k] = positions[i][k] - positions[j][k];
                    // Apply minimum image convention
                    r_ij[k] -= (r_ij[k] / l).round() * l;
                }
                let r = (r_ij[0].powi(2) + r_ij[1].powi(2) + r_ij[2].powi(2)).sqrt();
                let force = lj_potential(r) / r;
                for k in 0..3 {
                    forces[i][k] += force * r_ij[k];
                    forces[j][k] -= force * r_ij[k];
                }
            }
        }

        // Verlet integration
        let mut positions_new = vec![[0.0; 3]; n];
        for i in 0..n {
            for j in 0..3 {
                positions_new[i][j] = 2.0 * positions[i][j] - positions_old[i][j] + forces[i][j] * dt.powi(2);
            }
        }

        // Update velocities and handle boundary collisions
        for i in 0..n {
            for j in 0..3 {
                // Check for boundary collisions
                if positions_new[i][j] >= l {
                    positions_new[i][j] = 2.0 * l - positions_new[i][j]; // Reflect position
                    velocities[i][j] = -velocities[i][j]; // Reverse velocity
                } else if positions_new[i][j] <= 0.0 {
                    positions_new[i][j] = -positions_new[i][j]; // Reflect position
                    velocities[i][j] = -velocities[i][j]; // Reverse velocity
                }
            }
        }

        // Update velocities
        for i in 0..n {
            for j in 0..3 {
                velocities[i][j] = (positions_new[i][j] - positions_old[i][j]) / (2.0 * dt);
            }
        }

        // Calculate the current temperature
        let mut kinetic_energy = 0.0;
        for vel in &velocities {
            kinetic_energy += 0.5 * mass_argon * (vel[0].powi(2) + vel[1].powi(2) + vel[2].powi(2));
        }
        let current_temperature = (2.0 * kinetic_energy) / (3.0 * n as f64 * kb);

        // Calculate the scaling factor
        let scaling_factor = (1.0 + dt / tau * (target_temperature / current_temperature - 1.0)).sqrt();

        // Scale the velocities
        for vel in velocities.iter_mut() {
            for coord in vel.iter_mut() {
                *coord *= scaling_factor;
            }
        }

        // Store trajectory data
        if step % snapshot_interval == 0 {
            trajectory.push(positions.clone());
        }

        // Update positions
        positions_old.copy_from_slice(&positions);
        positions.copy_from_slice(&positions_new);
    }

    // Finish the progress bar
    pb.finish_with_message("Simulation complete");

    // Prepare simulation data for output
    let simulation_data = SimulationData {
        box_length: l,
        num_atoms: n,
        timestep: dt,
        total_steps: steps,
        snapshot_interval,
        trajectory,
    };

    // Save simulation data to JSON file
    let json = serde_json::to_string(&simulation_data).unwrap();
    let mut file = File::create("simulation_data.json").unwrap();
    file.write_all(json.as_bytes()).unwrap();

    println!("Simulation completed. Data saved to simulation_data.json");
}