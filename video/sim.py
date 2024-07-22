import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from moviepy.editor import ImageSequenceClip

# Load simulation data from JSON file
with open('simulation_data.json', 'r') as f:
    simulation_data = json.load(f)

# Extract simulation parameters
L = simulation_data['box_length']
N = simulation_data['num_atoms']
DT = simulation_data['timestep']
STEPS = simulation_data['total_steps']
SNAPSHOT_INTERVAL = simulation_data['snapshot_interval']

# Convert trajectory to numpy array
trajectory = np.array(simulation_data['trajectory'])

# Visualize the system
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_zlim(0, L)
ax.set_xlabel('X (angstrom)')
ax.set_ylabel('Y (angstrom)')
ax.set_zlabel('Z (angstrom)')
ax.set_title('Argon System Simulation')

# Create a movie of the time evolution
frames = []
for i, positions in enumerate(trajectory):
    ax.clear()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_zlim(0, L)
    ax.set_xlabel('X (angstrom)')
    ax.set_ylabel('Y (angstrom)')
    ax.set_zlabel('Z (angstrom)')
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c='b', marker='o', label='Argon Atoms')
    ax.legend()
    plt.title(f'Timestep {i * SNAPSHOT_INTERVAL}')
    plt.tight_layout()
    fig.canvas.draw()
    frame = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    frames.append(frame)

# Save the movie
clip = ImageSequenceClip(frames, fps=10)
clip.write_videofile('argon_simulation.mp4')

plt.show()

print(f"Visualization completed. Movie saved as argon_simulation.mp4")