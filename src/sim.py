import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from moviepy.editor import ImageSequenceClip
from tqdm import tqdm
import os
from multiprocessing import Pool

def get_file_size(file_path):
    return os.path.getsize(file_path)

def process_frame(args):
    i, positions = args
    frame_path = os.path.join(temp_dir, f'frame_{i:04d}.png')
    ax.clear()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_zlim(0, L)
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c='b', marker='o', label='Argon Atoms')
    ax.legend()
    plt.title(f'Timestep {i * SNAPSHOT_INTERVAL}')
    plt.tight_layout()
    plt.savefig(frame_path)
    return frame_path

# Get file size
file_path = 'simulation_data.json'
file_size = get_file_size(file_path)
print(f"Total file size: {file_size / (1024 * 1024):.2f} MB")

# Load simulation data from JSON file
with open(file_path, 'r') as f:
    simulation_data = json.load(f)

# Extract simulation parameters
L = simulation_data['box_length']
N = simulation_data['num_atoms']
DT = simulation_data['timestep']
STEPS = simulation_data['total_steps']
SNAPSHOT_INTERVAL = simulation_data['snapshot_interval']

# Convert trajectory to numpy array
trajectory = np.array(simulation_data['trajectory'])

# Create a directory for frames
temp_dir = 'temp'
os.makedirs(temp_dir, exist_ok=True)

# Visualize the system
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create a list of arguments for multiprocessing
args = [(i, trajectory[i]) for i in range(len(trajectory))]

# Process frames in parallel
with Pool() as pool:
    frame_paths = list(tqdm(pool.imap(process_frame, args), total=len(args), desc="Processing frames"))

print("Creating video...")

# Save the movie from the temp directory
clip = ImageSequenceClip(temp_dir, fps=10)
clip.write_videofile('argon_simulation.mp4')

plt.close(fig)

# Clean up temporary files
for file in os.listdir(temp_dir):
    os.remove(os.path.join(temp_dir, file))

# Remove the temp directory
os.rmdir(temp_dir)

print(f"Visualization completed. Movie saved as argon_simulation.mp4")
