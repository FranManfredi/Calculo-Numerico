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
dt = 10  # Time step (s)
num_steps = 10000000  # Number of steps for simulation

# Lists to store trajectory data for plotting
r_values = []
theta_values = []

# Secant method values
v_r_guess1 = 500
v_r_guess2 = 600  # Un poco más alto para comenzar
v_theta_guess1 = 700
v_theta_guess2 = 800  # Un poco más alto para comenzar


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
        if distance_to_asteroid < 1e7:
            print("Target reached!")
            break
        elif distance_to_earth < 0.01e9:
            print("Target crashed into Earth!")
            break


def final_distance(v_r_initial, v_theta_initial, min_steps=600000, max_steps=700000, step_interval=10000):
    """
    Calcula la distancia mínima al asteroide en un rango de pasos.
    """
    # Condiciones iniciales
    r = r_satelite
    theta = theta_satelite
    v_r = v_r_initial
    v_theta = v_theta_initial

    min_distance = float('inf')  # Variable para almacenar la distancia mínima

    for steps in range(min_steps, max_steps, step_interval):
        # Resetear condiciones iniciales para cada rango de steps
        r = r_satelite
        theta = theta_satelite
        v_r = v_r_initial
        v_theta = v_theta_initial

        # Simulación para el número actual de steps
        for _ in range(steps):
            # Calcular aceleración radial en la posición actual
            a_r = radial_acceleration(r)

            # Método RK2 simplificado para actualizar solo el final
            v_r_half = v_r + a_r * (dt / 2)
            r_half = r + v_r * (dt / 2)
            v_theta_half = v_theta * (r / r_half)

            # Actualizar posición y velocidad
            r += v_r_half * dt
            theta += v_theta_half * dt / r

            a_r_end = radial_acceleration(r)
            v_r += a_r_end * dt

            theta %= (2 * np.pi)  # Asegurar que theta esté entre 0 y 2π

        # Calcular distancia al asteroide al final de los steps
        distance_to_asteroid = np.sqrt(
            (r - r_asteroide) ** 2 + (r * theta - r_asteroide * theta_asteroide) ** 2
        )

        # Actualizar la distancia mínima si es necesario
        if distance_to_asteroid < min_distance:
            min_distance = distance_to_asteroid

    return min_distance


# Método de la secante para encontrar velocidades óptimas
def secant_method(v_r_guess1, v_r_guess2, v_theta_guess1, v_theta_guess2, tol=1e7):
    # Evaluar la función de error para ambos conjuntos de conjeturas
    f1 = final_distance(v_r_guess1, v_theta_guess1)
    f2 = final_distance(v_r_guess2, v_theta_guess2)

    iteration = 0
    while abs(f2) > tol and iteration < 25:
        # Actualizar con las ecuaciones de la secante
        v_r_new = v_r_guess2 - f2 * (v_r_guess2 - v_r_guess1) / (f2 - f1)
        v_theta_new = v_theta_guess2 - f2 * (v_theta_guess2 - v_theta_guess1) / (f2 - f1)

        # Actualizar valores anteriores y calcular nuevo error
        v_r_guess1, v_theta_guess1 = v_r_guess2, v_theta_guess2
        f1 = f2
        v_r_guess2, v_theta_guess2 = v_r_new, v_theta_new
        f2 = final_distance(v_r_guess2, v_theta_guess2)

        iteration += 1
        print(f"Iteration {iteration}: Distance to asteroid = {f2:.2f} meters | Error: {abs(f2 - tol):.2f}\n"
              f"Initial velocity for projectile found: Theta: {v_theta_guess2:.2f} | Radial: {v_r_guess2:.2f} meters\n")

    return v_r_guess2, v_theta_guess2


# Ejecutar el método de la secante

print("\n\nmetodo secante:\n")
optimal_v_r, optimal_v_theta = secant_method(v_r_guess1, v_r_guess2, v_theta_guess1, v_theta_guess2, 9e6)
print()
print(f"Optimal radial velocity: {optimal_v_r:.2f} m/s")
print(f"Optimal angular velocity: {optimal_v_theta:.2f} m/s")

print("\n--------------------------------------------------------------------------------\n\n"
      "Ejecucion de rk2 con grafico:\n")

rk2(r_satelite, theta_satelite, optimal_v_r, optimal_v_theta)

# Plot the satellite's trajectory
plt.figure(figsize=(8, 8))
plt.polar(theta_values, r_values, label="Trayectoria del Satélite")
plt.polar([theta_asteroide], [r_asteroide], 'ro', label="Posición del Asteroide", markersize=10)
plt.title("Trayectoria del Satélite hacia el Asteroide Objetivo")
plt.legend()
plt.show()

