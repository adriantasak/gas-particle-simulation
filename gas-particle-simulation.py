import tkinter as tk
import numpy as np
from scipy.spatial.distance import pdist
from numpy.linalg import norm

# Set up Tkinter window
root = tk.Tk()
root.title("Gas Particle Collision Simulation")
w, h = 700, 700
canvas = tk.Canvas(root, width=w, height=h)
canvas.pack()

# Set simulation parameters
particle_count = int(input("Enter particle count (default is 1000):") or 1000)
radius = float(input("Enter particle radius (default is 3):") or 3)
velocity = float(input("Enter particle velocity (default is 10):") or 10)
time_step = float(input("Enter time step (default is 0.2):") or 0.2)

# Generate random initial positions and velocities for particles
particle_state = np.zeros((particle_count, 4))   # 4-vector containing the state of each particles
particle_state[:, 0] = np.random.randint(2, w - radius, particle_count)   # position on x-axis
particle_state[:, 1] = np.random.randint(2, h - radius, particle_count)   # position on y-axis
motion_angle = np.random.uniform(0, 100, particle_count)
particle_state[:, 2] = velocity * np.sin(motion_angle)   # velocity on x-axis
particle_state[:, 3] = velocity * np.cos(motion_angle)   # velocity on y-axis

# Create and display particles on the canvas
particles = []
for i in range(particle_count):
    x1, y1, x2, y2 = (
        particle_state[i, 0] - radius, 
        particle_state[i, 1] - radius, 
        particle_state[i, 0] + radius, 
        particle_state[i, 1] + radius
    )
    particle = canvas.create_oval(x1, y1, x2, y2, fill='blue')
    particles.append(particle)

# Get indices for upper triangle of a (particle_count x particle_count) matrix
upper_triangle_indices = np.triu_indices(particle_count, k=1)

def step(time_step, particle_state):
    """
    Update positions and velocities of particles.
    
    1. First, the current particle state is copied for later comparison. 
       The particle positions are then updated based on their current velocities and the time step.

    2. Next, the distances between all pairs of particles are calculated.
    
    3. Pairs of particles that are colliding (their distance is less than or equal to twice the particle radius) 
       are identified.

    4. We ensure that each pair of colliding particles only contains distinct particles. 
       In other words, we filter out pairs that contain the same particle.

    5. Then, we calculate the new velocities and positions for colliding particles. 
       For each pair of colliding particles, we calculate the difference in position and velocity. 
       We also compute a scalar product of the velocity difference and position difference. 
       Then, the velocities are updated based on these quantities. The positions are also adjusted to 
       ensure particles do not overlap after collision.
       
    6. Next, we make particles rebound off the boundaries of the canvas. 
       If a particle is going to move beyond the left or top edges of the canvas (position less than the particle radius),
       or beyond the right or bottom edges of the canvas (position greater than canvas size minus the particle radius), 
       we invert its velocity (making it move in the opposite direction) and adjust its position to ensure it stays within the canvas.

    7. Finally, we update the positions of particles on the canvas. We calculate the change in position 
       (current position minus old position) and move each particle accordingly on the canvas.
    """
    
    state_copy = particle_state.copy()
    particle_state[:, :2] += particle_state[:, 2:] * time_step

    # Find distances between all pairs of particles
    distances = pdist(particle_state[:, :2])
    
    # Find pairs of particles that are colliding
    is_colliding = distances <= 2 * radius
    colliding_1 = upper_triangle_indices[0][is_colliding]
    colliding_2 = upper_triangle_indices[1][is_colliding]

    # Make sure each pair only contains distinct particles
    unique_1, mask = np.unique(colliding_1, return_index=True)
    colliding_1, colliding_2 = colliding_1[mask], colliding_2[mask]
    unique_2, mask = np.unique(colliding_2, return_index=True)
    colliding_1, colliding_2 = colliding_1[mask], colliding_2[mask]
    distinct_pairs = np.in1d(colliding_1, colliding_2, invert=True)
    colliding_1, colliding_2 = colliding_1[distinct_pairs], colliding_2[distinct_pairs]

    # Compute new velocities and positions for colliding particles
    if colliding_1.size:
        pos_1, pos_2 = particle_state[colliding_1, :2], particle_state[colliding_2, :2]
        vel_1, vel_2 = particle_state[colliding_1, 2:], particle_state[colliding_2, 2:]
        pos_diff = pos_2 - pos_1
        vel_diff = vel_1 - vel_2
        scalar_product = np.einsum('ij,ij->i', vel_diff, pos_diff)[:, np.newaxis]
        pos_diff_norm = norm(pos_diff, axis=1)[:, np.newaxis]

        particle_state[colliding_1, 2:] = vel_1 - (scalar_product / pos_diff_norm ** 2) * pos_diff
        particle_state[colliding_2, 2:] = vel_2 + (scalar_product / pos_diff_norm ** 2) * pos_diff

        particle_state[colliding_1, :2] = pos_1 - 0.5 * ((2 * radius / pos_diff_norm) - 1) * pos_diff
        particle_state[colliding_2, :2] = pos_2 + 0.5 * ((2 * radius / pos_diff_norm) - 1) * pos_diff

    # Make particles rebound off walls
    for i in range(2):
        size = [w, h][i]
        out_of_bounds_min = particle_state[:, i] < radius
        out_of_bounds_max = particle_state[:, i] > (size - radius)

        if out_of_bounds_min.sum() > 0:
            particle_state[out_of_bounds_min, i + 2] = -particle_state[out_of_bounds_min, i + 2]
            particle_state[out_of_bounds_min, i] += radius - particle_state[out_of_bounds_min, i] + 1

        if out_of_bounds_max.sum() > 0:
            particle_state[out_of_bounds_max, i + 2] = -particle_state[out_of_bounds_max, i + 2]
            particle_state[out_of_bounds_max, i] += (size - radius) - particle_state[out_of_bounds_max, i] - 1

    # Update particles' positions on the canvas
    pos_diff = particle_state[:, :2] - state_copy[:, :2]
    for j in range(particle_count):
        canvas.move(particles[j], pos_diff[j, 0], pos_diff[j, 1])

def animate():
    """Continuously update the simulation every 40 ms using the step function."""
    step(time_step, particle_state)
    root.after(40, animate)

# Start animation
root.after(0, animate)
root.mainloop()
