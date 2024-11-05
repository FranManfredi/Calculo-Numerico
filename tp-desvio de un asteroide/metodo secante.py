import numpy as np

# Constantes
G = 6.67430e-11  # Constante gravitacional (m^3 kg^-1 s^-2)
M = 5.972e24  # Masa de la Tierra (kg)

# Condiciones iniciales del satélite
r_satelite = 1.22e9  # Posición radial inicial del satélite (m)
theta_satelite = 0.0  # Posición angular inicial del satélite (rad)

# Posición del asteroide
r_asteroide = 2.24775e9  # Posición radial del asteroide (m)
theta_asteroide = 2.303834613  # Posición angular del asteroide (rad)

# Parámetros de simulación (optimizados)
dt = 10  # Aumentado el paso de tiempo (s)
num_steps = 604035  # Reducido el número de pasos para la simulación

# Función de aceleración radial debido a la gravedad
def radial_acceleration(r):
    return -G * M / r ** 2

# Función para evaluar solo el resultado final de la simulación
def final_distance(v_r_initial, v_theta_initial):
    # Condiciones iniciales
    r = r_satelite
    theta = theta_satelite
    v_r = v_r_initial
    v_theta = v_theta_initial

    # Iterar sin almacenar datos intermedios, solo calculando el estado final
    for _ in range(num_steps):
        # Calcular aceleración radial en la posición actual
        a_r = radial_acceleration(r)

        # Método RK2 simplificado para actualizar solo el final
        v_r_half = v_r + a_r * (dt / 2)
        r_half = r + v_r * (dt / 2)

        v_theta_half = v_theta * (r / r_half)
        theta_half = theta + v_theta_half * (dt / 2)

        r += v_r_half * dt
        theta += v_theta_half * dt / r

        a_r_end = radial_acceleration(r)
        v_r += a_r_end * dt

        theta %= (2 * np.pi)

    # Calcular la distancia final al asteroide
    distance_to_asteroid = np.sqrt(
        (r - r_asteroide) ** 2 + (r * theta - r_asteroide * theta_asteroide) ** 2
    )
    return distance_to_asteroid

# Método de la secante para encontrar velocidades óptimas
def secant_method(v_r_guess1, v_r_guess2, v_theta_guess1, v_theta_guess2, tol=0.9e7):
    # Evaluar la función de error para ambos conjuntos de conjeturas
    f1 = final_distance(v_r_guess1, v_theta_guess1)
    f2 = final_distance(v_r_guess2, v_theta_guess2)

    iteration = 0
    while abs(f2) > tol and iteration < 200:
        # Actualizar con las ecuaciones de la secante
        v_r_new = v_r_guess2 - f2 * (v_r_guess2 - v_r_guess1) / (f2 - f1)
        v_theta_new = v_theta_guess2 - f2 * (v_theta_guess2 - v_theta_guess1) / (f2 - f1)

        # Actualizar valores anteriores y calcular nuevo error
        v_r_guess1, v_theta_guess1 = v_r_guess2, v_theta_guess2
        f1 = f2
        v_r_guess2, v_theta_guess2 = v_r_new, v_theta_new
        f2 = final_distance(v_r_guess2, v_theta_guess2)

        iteration += 1
        print(f"Iteration {iteration}: Distance to asteroid = {f2:.2f} meters")

    return v_r_guess2, v_theta_guess2

# Conjeturas iniciales para el método de la secante
v_r_guess1 = 500
v_r_guess2 = 600  # Un poco más alto para comenzar
v_theta_guess1 = 700
v_theta_guess2 = 800  # Un poco más alto para comenzar

# Ejecutar el método de la secante
optimal_v_r, optimal_v_theta = secant_method(v_r_guess1, v_r_guess2, v_theta_guess1, v_theta_guess2)
print(f"Optimal radial velocity: {optimal_v_r:.2f} m/s")
print(f"Optimal angular velocity: {optimal_v_theta:.2f} m/s")
