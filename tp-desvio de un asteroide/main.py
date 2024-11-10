import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M = 5.972e24  # Mass of Earth (kg)

# Initial conditions for the satellite
r_satelite = 1.22e9  # Initial radial position of the satellite (m)
theta_satelite = 0.0  # Initial position of the satellite (rad)

# Target (asteroid) position
r_asteroide = 2.24775e9  # Radial position of the asteroid (m)
theta_asteroide = 2.303834613  # Angular position of the asteroid (rad)

# Simulation parameters
dt = 1  # Time step (s)
num_steps = 50000000  # Number of steps for simulation

# Lists to store trajectory data for plotting
r_values = []
theta_values = []

# Secant method values
v_r_guess1 = 500
v_r_guess2 = 700  # Small initial radial velocities
v_theta_guess1 = np.sqrt(G * M / r_satelite) * 0.9
v_theta_guess2 = np.sqrt(G * M / r_satelite) * 1.1


# Function to calculate radial acceleration due to gravity
def radial_acceleration(r):
    return -G * M / r ** 2


def rk2(r_satelite, theta_satelite, v_r_satelite, v_theta_satelite):
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
        earth_radius = 6.371e6  # Radius of Earth in meters
        if distance_to_asteroid < 5e3:
            print("Target reached!")
            break
        elif distance_to_earth <= earth_radius:
            print("Target crashed into Earth!")
            break


from scipy.optimize import minimize_scalar
def final_distance(v_r_initial, v_theta_initial, dt=0.1, max_steps=200000):
    r_sim = r_satelite
    theta_sim = theta_satelite
    v_r_sim = v_r_initial
    v_theta_sim = v_theta_initial
    min_distance = float('inf')

    # Run the simulation for a limited number of steps
    for _ in range(max_steps):
        # Calculate gravitational acceleration at current position
        a_r_sim = radial_acceleration(r_sim)

        # RK2 midpoint estimates
        v_r_half = v_r_sim + a_r_sim * (dt / 2)
        r_half = r_sim + v_r_sim * (dt / 2)
        v_theta_half = v_theta_sim * (r_sim / r_half)

        # Update the radial and angular positions
        r_sim += v_r_half * dt
        theta_sim += v_theta_half * dt / r_sim

        # Update the radial velocity for the end of the time step
        a_r_end = radial_acceleration(r_sim)
        v_r_sim += a_r_end * dt
        theta_sim %= (2 * np.pi)

        # Calculate distance between the satellite and the asteroid in polar coordinates
        delta_r = abs(r_sim - r_asteroide)
        delta_theta = abs(theta_sim - theta_asteroide)
        distance_to_asteroid = np.sqrt(delta_r ** 2 + (r_sim * delta_theta) ** 2)

        # Update minimum distance encountered
        min_distance = min(min_distance, distance_to_asteroid)

        # Check for proximity to target asteroid
        if distance_to_asteroid < 5e3:
            print("Asteroid proximity reached in simulation.")
            break

    return min_distance

def hybrid_secant_bisection_method(v_r_guess1, v_r_guess2, v_theta_guess1, v_theta_guess2, tol=1e6, max_step=25, max_iterations=150):
    f1 = final_distance(v_r_guess1, v_theta_guess1)
    f2 = final_distance(v_r_guess2, v_theta_guess2)

    for iteration in range(max_iterations):
        if abs(f2) < tol:
            print("Converged!")
            break

        # Smaller, stable range for bisection adjustments to avoid runaway velocities
        if f2 > f1:
            v_r_new = (v_r_guess1 + v_r_guess2) / 2
            v_theta_new = (v_theta_guess1 + v_theta_guess2) / 2
            print("Switching to bisection for stability.")
        else:
            v_r_new = v_r_guess2 - f2 * (v_r_guess2 - v_r_guess1) / (f2 - f1)
            v_theta_new = v_theta_guess2 - f2 * (v_theta_guess2 - v_theta_guess1) / (f2 - f1)
            v_r_new = np.clip(v_r_new, v_r_guess2 - max_step, v_r_guess2 + max_step)
            v_theta_new = np.clip(v_theta_new, v_theta_guess2 - max_step, v_theta_guess2 + max_step)

        f1, f2 = f2, final_distance(v_r_new, v_theta_new)
        v_r_guess1, v_theta_guess1 = v_r_guess2, v_theta_guess2
        v_r_guess2, v_theta_guess2 = v_r_new, v_theta_new

        print(f"Iteration {iteration+1}: Distance to asteroid = {f2:.2f} meters\n"
              f"Current initial velocity guess: Theta: {v_theta_guess2:.2f} | Radial: {v_r_guess2:.2f} m/s\n")

    return v_r_guess2, v_theta_guess2


# Ejecutar el método de la secante

print("\n\nmetodo secante:\n")
optimal_v_r, optimal_v_theta = hybrid_secant_bisection_method(v_r_guess1, v_r_guess2, v_theta_guess1, v_theta_guess2, 5000)
print()
print(f"Optimal radial velocity: {optimal_v_r:.2f} m/s")
print(f"Optimal angular velocity: {optimal_v_theta:.2f} m/s")

print("\n--------------------------------------------------------------------------------\n\n"
      "Ejecucion de rk2 con grafico:\n")

rk2(r_satelite, theta_satelite, optimal_v_r, optimal_v_theta)

# Plot the satellite's trajectory
plt.figure(figsize=(8, 8))
plt.polar(theta_values, r_values, label="Trayectoria del Satélite")
plt.polar([theta_asteroide], [r_asteroide], 'ro', label="Posición del Asteroide", markersize=1)
plt.title("Trayectoria del Satélite hacia el Asteroide Objetivo")
plt.legend()
plt.show()
