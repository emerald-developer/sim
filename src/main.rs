use rayon::prelude::*;
use rand::Rng;
use serde::Serialize;
use std::fs::File;
use std::io::Write;
use std::env;
use indicatif::{ProgressBar, ProgressStyle, HumanDuration};
use std::time::{Instant, Duration};

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
    let mut positions = (0..n).map(|_| {
        [rng.gen::<f64>() * l, rng.gen::<f64>() * l, rng.gen::<f64>() * l]
    }).collect::<Vec<_>>();

    let mass_argon: f64 = 39.95;
    let kb: f64 = 0.0083144621;
    let velocity_factor = (kb * target_temperature / mass_argon).sqrt();
    let mut velocities = (0..n).map(|_| {
        [
            rng.gen::<f64>() * velocity_factor,
            rng.gen::<f64>() * velocity_factor,
            rng.gen::<f64>() * velocity_factor
        ]
    }).collect::<Vec<_>>();

    let mut positions_old = positions.clone();

    let mut trajectory = Vec::new();

    let pb = ProgressBar::new(steps as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
        .unwrap()
        .progress_chars("##-"));

    let start_time = Instant::now();
    let mut last_update = start_time;
    let update_interval = Duration::from_secs(1);

    // Perform simulation
    for step in 0..steps {
        pb.set_position(step as u64);

        // Calculate forces in parallel
        let forces: Vec<_> = (0..n).into_par_iter().map(|i| {
            let mut force = [0.0; 3];
            for j in 0..n {
                if i != j {
                    let mut r_ij = [0.0; 3];
                    for k in 0..3 {
                        r_ij[k] = positions[i][k] - positions[j][k];
                        r_ij[k] -= (r_ij[k] / l).round() * l;
                    }
                    let r = (r_ij[0].powi(2) + r_ij[1].powi(2) + r_ij[2].powi(2)).sqrt();
                    let force_magnitude = lj_potential(r) / r;
                    for k in 0..3 {
                        force[k] += force_magnitude * r_ij[k];
                    }
                }
            }
            force
        }).collect();

        // Verlet integration and boundary handling in parallel
        let (positions_new, new_velocities): (Vec<_>, Vec<_>) = positions.par_iter().zip(positions_old.par_iter()).zip(forces.par_iter()).zip(velocities.par_iter())
            .map(|(((pos, pos_old), force), _vel)| {
                let mut pos_new = [0.0; 3];
                let mut vel_new = [0.0; 3];
                for j in 0..3 {
                    pos_new[j] = 2.0 * pos[j] - pos_old[j] + force[j] * dt.powi(2);
                    vel_new[j] = (pos_new[j] - pos_old[j]) / (2.0 * dt);

                    if pos_new[j] >= l {
                        pos_new[j] = 2.0 * l - pos_new[j];
                        vel_new[j] = -vel_new[j];
                    } else if pos_new[j] <= 0.0 {
                        pos_new[j] = -pos_new[j];
                        vel_new[j] = -vel_new[j];
                    }
                }
                (pos_new, vel_new)
            }).unzip();

        // Update positions and velocities
        positions_old = positions;
        positions = positions_new;
        velocities = new_velocities;

        // Calculate the current temperature
        let kinetic_energy: f64 = velocities.par_iter().map(|vel| {
            0.5 * mass_argon * (vel[0].powi(2) + vel[1].powi(2) + vel[2].powi(2))
        }).sum();
        let current_temperature = (2.0 * kinetic_energy) / (3.0 * n as f64 * kb);

        // Calculate the scaling factor and scale velocities
        let scaling_factor = (1.0 + dt / tau * (target_temperature / current_temperature - 1.0)).sqrt();
        velocities.par_iter_mut().for_each(|vel| {
            for coord in vel.iter_mut() {
                *coord *= scaling_factor;
            }
        });

        // Store trajectory data
        if step % snapshot_interval == 0 {
            trajectory.push(positions.clone());
        }

        // Update progress bar with time left and speed
        let now = Instant::now();
        if now.duration_since(last_update) >= update_interval {
            let elapsed = now.duration_since(start_time);
            let iterations_per_sec = step as f64 / elapsed.as_secs_f64();
            let estimated_total = Duration::from_secs_f64(steps as f64 / iterations_per_sec);
            let time_left = estimated_total.saturating_sub(elapsed);
            
            pb.set_message(format!(
                "Speed: {:.2} it/s | Time left: {}",
                iterations_per_sec,
                HumanDuration(time_left)
            ));
            
            last_update = now;
        }
    }

    pb.finish_with_message("Simulation complete");

    let simulation_data = SimulationData {
        box_length: l,
        num_atoms: n,
        timestep: dt,
        total_steps: steps,
        snapshot_interval,
        trajectory,
    };

    let json = serde_json::to_string(&simulation_data).unwrap();
    let mut file = File::create("simulation_data.json").unwrap();
    file.write_all(json.as_bytes()).unwrap();

    println!("Simulation completed. Data saved to simulation_data.json");
}