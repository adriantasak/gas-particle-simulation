# Gas Particle Collision Simulation

This simulation models the dynamics of gas particles confined within a 2D box. It assumes that all particles share the same mass and undergo elastic collisions with each other and with the box's walls, thereby adhering to the principles of conservation of energy and momentum. One notable feature of the program is its ability to naturally demonstrate the [Maxwell-Boltzmann distribution of particle speeds](https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution), a core concept in statistical mechanics. This distribution emerges as a result of the described dynamics, without requiring explicit coding into the simulation.

<img src="./demo.gif" width="480" height="480" />

## Usage

Upon execution, the script prompts the user for inputs to customize the simulation:
  - particle count (default is 1000)
  - particle radius (default is 3)
  - particle velocity (default is 10)
  - time step (default is 0.2)

Once the parameters are set, the simulation will begin in a new window.

## Behind the Scenes

After initializing particles with random positions and velocities, the `step` function drives the simulation by continuously updating particle states. It calculates new positions of particles based on their velocities and the time step, identifies colliding particle pairs, calculates new velocities and positions for these colliding pairs, makes particles rebound off the walls of the box, and finally updates the particles' positions on the canvas. The simulation uses the `Tkinter` library to visualize these interactions, providing a graphical interface to represent particle movement and collisions.

## Limitations

The simulation progresses in discrete time steps, leading to a potential issue where particles can briefly "stick" to each other. This arises when particles are still overlapping after a collision due to the time step size, triggering another collision in the subsequent time step. As a result, particles continually collide and change directions until they escape the overlap, causing a "spinning" effect. The simulation attempts to address this by repositioning the particles after a collision is detected to prevent overlap. However, this issue may still arise, especially with larger time steps or high particle speeds.
