import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M = 5.972e24  # Mass of Earth (kg)

# Initial conditions for the satellite
r_satelite = 1.22e9  # Initial radial position of the satellite (m)
theta_satelite = 0.0  # Initial angular position of the satellite (rad)
v_r_satelite = 553.0 # Initial radial velocity (m/s) - increased for testing
v_theta_satelite = 753.0  # Initial angular velocity (m/s)

# Optimal radial velocity: 551.69 m/s
# Optimal angular velocity: 751.69 m/s

# Target (asteroid) position
r_asteroide = 2.24775e9  # Radial position of the asteroid (m)
theta_asteroide = 2.303834613  # Angular position of the asteroid (rad)

# Simulation parameters
dt = 10  # Time step (s)
num_steps = 1000000  # Number of steps for simulation

# Lists to store trajectory data for plotting
r_values = []
theta_values = []


# Function to calculate radial acceleration due to gravity
def radial_acceleration(r):
    return -G * M / r ** 2


# Calculate escape velocity at initial radius
v_escape_initial = np.sqrt(2 * G * M / r_satelite)
initial_total_velocity = np.sqrt(v_r_satelite ** 2 + v_theta_satelite ** 2)

# Runge-Kutta 2nd Order (RK2) simulation
for _ in range(num_steps):
    # Store current position for plotting
    r_values.append(r_satelite)
    theta_values.append(theta_satelite)

    # Calculate radial gravitational acceleration at the current position
    a_r = radial_acceleration(r_satelite)

    # RK2 midpoint estimates for radial and angular updates
    # Midpoint velocities and positions
    v_r_half = v_r_satelite + a_r * (dt / 2)
    r_half = r_satelite + v_r_satelite * (dt / 2)

    # Angular velocity at the midpoint, assuming conservation of angular momentum
    v_theta_half = v_theta_satelite * (r_satelite / r_half)
    theta_half = theta_satelite + v_theta_half * (dt / 2)

    # Full-step updates for r and theta using midpoint values
    r_satelite += v_r_half * dt
    theta_satelite += v_theta_half * dt / r_satelite  # Adjust theta based on current r

    # Update radial velocity at the end of the time step
    a_r_end = radial_acceleration(r_satelite)  # Recalculate acceleration after moving r
    v_r_satelite += a_r_end * dt

    # Ensure theta stays within [0, 2π]
    theta_satelite %= (2 * np.pi)

    # Check if the satellite is close to the asteroid
    distance_to_asteroid = np.sqrt(
        (r_satelite - r_asteroide) ** 2 + (r_satelite * theta_satelite - r_asteroide * theta_asteroide) ** 2
    )
    distance_to_earth = r_satelite
    if distance_to_asteroid < 1e7:
        print("Target reached!")
        break
    elif distance_to_earth < 0.01e9:
        print("Target crashed into Earth!")
        break

# Plot the satellite's trajectory
plt.figure(figsize=(8, 8))
plt.polar(theta_values, r_values, label="Satellite Trajectory")
plt.polar([theta_asteroide], [r_asteroide], 'ro', label="Asteroid Position")
plt.title("Satellite Trajectory to Target Asteroid")
plt.legend()
plt.show()
