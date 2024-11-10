import numpy as np
import matplotlib.pyplot as plt

# Constantes
G = 6.67430e-11  # Constante gravitacional
M = 5.972e24  # Masa de la Tierra
h = 0.01  # Paso de tiempo
max_secant_iter = 4500 # Número máximo de iteraciones secante antes de cambiar a bisección


def get_initial_velocity():
    # Conjeturas iniciales para la velocidad radial y angular
    vr1, vr2 = 800, 100
    vtheta1, vtheta2 = 750, 800
    i = 0

    while abs(vr2 - vr1) > 1e-4 or abs(vtheta2 - vtheta1) > 1e-7:
        dist1 = get_trajectory_distance(vr1, vtheta1)
        dist2 = get_trajectory_distance(vr2, vtheta2)

        # Si el método de la secante no converge, cambiar a bisección
        if i >= max_secant_iter:
            vr3 = (vr1 + vr2) / 2
            vtheta3 = (vtheta1 + vtheta2) / 2
        else:
            vr3 = secante(vr1, vr2, dist1, dist2)
            vtheta3 = secante(vtheta1, vtheta2, dist1, dist2)

        # Limitar el tamaño del salto en vr3 y vtheta3 para evitar cambios grandes
        step_limit = 1e5
        if abs(vr3 - vr2) > step_limit:
            vr3 = vr2 + np.sign(vr3 - vr2) * step_limit
        if abs(vtheta3 - vtheta2) > step_limit:
            vtheta3 = vtheta2 + np.sign(vtheta3 - vtheta2) * step_limit

        vr1, vr2 = vr2, vr3
        vtheta1, vtheta2 = vtheta2, vtheta3
        i += 1
        print(f"Iteration {i}: vr = {vr3}, vtheta = {vtheta3}, g(vr, vtheta) = {dist2}")

        # Comprobación de convergencia en la distancia
        if dist2 < 1e-4:
            break

    return vr2, vtheta2


def get_trajectory_distance(vr_initial, vtheta_initial):
    r = 1.22e9  # Distancia radial inicial
    theta = 0
    vr = vr_initial
    vtheta = vtheta_initial
    max_time = 10
    target_r, target_theta = 2.24775e9, 2.303834613

    time = 0
    while time < max_time:
        # Método RK2 para actualizar velocidad radial y posición
        ar = -G * M / r ** 2
        vr_mid = vr + 0.5 * ar * h
        r_mid = r + 0.5 * vr * h

        # RK2 para la componente angular
        vtheta_mid = vtheta
        theta_mid = theta + 0.5 * vtheta * h / r

        # Actualización de posiciones y velocidades
        vr += (-G * M / r_mid ** 2) * h
        vtheta = vtheta_mid  # No cambia la velocidad angular
        r += vr * h
        theta += vtheta * h / r if r != 0 else 0

        time += h

    # Calcular la distancia al objetivo al final de la trayectoria
    return np.sqrt((r * np.cos(theta) - target_r * np.cos(target_theta)) ** 2 +
                   (r * np.sin(theta) - target_r * np.sin(target_theta)) ** 2)


# Método de la secante para actualizar la velocidad radial o angular
def secante(v1, v2, g1, g2):
    return v2 - g2 * (v2 - v1) / (g2 - g1)


# Calcular las velocidades radial y angular iniciales requeridas
required_vr, required_vtheta = get_initial_velocity()
print(f"La velocidad radial inicial requerida es aproximadamente: {required_vr} m/s")
print(f"La velocidad angular inicial requerida es aproximadamente: {required_vtheta} m/s")


def plot_trajectory(vr_initial, vtheta_initial, target_r, target_theta):
    r = 1.22e9  # Distancia radial inicial (m)
    theta = 0  # Ángulo inicial (rad)
    vr = vr_initial  # Velocidad radial inicial
    vtheta = vtheta_initial  # Velocidad angular inicial
    max_time = 1000  # Tiempo máximo de simulación (s)

    time = 0
    trajectory_r = [r]  # Lista para almacenar la distancia radial a lo largo del tiempo
    trajectory_theta = [theta]  # Lista para almacenar el ángulo a lo largo del tiempo

    # Simulación de la trayectoria utilizando el método RK2
    while time < max_time:
        ar = -G * M / r ** 2
        vr_mid = vr + 0.5 * ar * h
        r_mid = r + 0.5 * vr * h

        # RK2 para la componente angular
        vtheta_mid = vtheta
        theta_mid = theta + 0.5 * vtheta * h / r

        # Actualización de posiciones y velocidades
        vr += (-G * M / r_mid ** 2) * h
        vtheta = vtheta_mid  # La velocidad angular no cambia
        r += vr * h
        theta += vtheta * h / r if r != 0 else 0

        # Agregar los valores actuales a la lista de trayectoria
        trajectory_r.append(r)
        trajectory_theta.append(theta)

        time += h

    # Convertir la trayectoria a coordenadas cartesianas para la gráfica
    x_traj = np.array(trajectory_r) * np.cos(trajectory_theta)
    y_traj = np.array(trajectory_r) * np.sin(trajectory_theta)

    # Posición del objetivo
    target_x = target_r * np.cos(target_theta)
    target_y = target_r * np.sin(target_theta)

    # Crear la gráfica
    plt.figure(figsize=(8, 8))
    plt.plot(x_traj, y_traj, label="Projectile Trajectory")
    plt.plot(target_x, target_y, 'ro', label="Target Position (Asteroid)")
    plt.xlabel("X Distance (m)")
    plt.ylabel("Y Distance (m)")
    plt.title("Projectile Trajectory towards Target")
    plt.legend()
    plt.grid()
    plt.axis('equal')  # Para que las escalas de los ejes X e Y sean iguales
    plt.show()


target_r, target_theta = 2.24775e9, 2.303834613

plot_trajectory(required_vr, required_vtheta, target_r, target_theta)