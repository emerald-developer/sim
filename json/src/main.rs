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

fn apply_pbc(mut position: f64, l: f64) -> f64 {
    position -= (position / l).floor() * l;
    position
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
    let temperature: f64 = 87.3;
    let velocity_factor = (kb * temperature / mass_argon).sqrt();
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

        // Apply periodic boundary conditions
        for pos in positions.iter_mut() {
            for coord in pos.iter_mut() {
                *coord = apply_pbc(*coord, l);
            }
        }

        // Calculate forces
        let mut forces = vec![[0.0; 3]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let mut r_ij = [0.0; 3];
                for k in 0..3 {
                    r_ij[k] = positions[i][k] - positions[j][k];
                    r_ij[k] -= (r_ij[k] / l).round() * l; // Apply minimum image convention
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
                if positions_new[i][j] > l {
                    positions_new[i][j] -= l;
                } else if positions_new[i][j] < 0.0 {
                    positions_new[i][j] += l;
                }
            }
        }

        // Update velocities
        for i in 0..n {
            for j in 0..3 {
                velocities[i][j] = (positions_new[i][j] - positions_old[i][j]) / (2.0 * dt);
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