import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Constants
G = 6.67430e-11  # Gravitational constant
M = 5.972e24  # Mass of the Earth
R_earth = 6.371e6  # Radius of the Earth (in meters)
h = 1  # Time step (increased for speed)
target_r, target_theta = 2.24775e9, 2.303834613
max_time = 24642  # Maximum simulation time

initial_theta_guess = 280000
initial_r_guess = initial_theta_guess/1.22e9

def get_trajectory_distance(v_initial):
    """Calculate final distance from the target given initial velocities (vr, vtheta)."""
    vr, vtheta = v_initial
    r = 1.22e9  # Initial radial distance
    theta = 0  # Initial angle
    time = 0

    while time < max_time:
        # Calculate the acceleration and update velocities and position (RK-like)
        ar = -G * M / r ** 2
        vr_mid = vr + 0.5 * ar * h
        r_mid = r + 0.5 * vr * h
        theta_mid = theta + 0.5 * vtheta * h / r if r != 0 else 0

        vr += (-G * M / r_mid ** 2) * h
        r += vr_mid * h
        theta += vtheta * h / r if r != 0 else 0

        # Compute distance to target
        distance_to_target = np.sqrt((r * np.cos(theta) - target_r * np.cos(target_theta)) ** 2 +
                                     (r * np.sin(theta) - target_r * np.sin(target_theta)) ** 2)

        # Check if within impact distance
        if distance_to_target <= 1000:
            return 0  # Successful impact

        time += h

    # Return final distance if no impact occurred
    return distance_to_target


# Optimize initial velocities
initial_guess = [initial_r_guess, initial_theta_guess]  # Initial guesses for vr and vtheta
result = minimize(get_trajectory_distance, initial_guess, method='Nelder-Mead')
required_vr, required_vtheta = result.x

print(f"The optimized initial radial velocity is approximately: {required_vr:.2f} m/s")
print(f"The optimized initial angular velocity is approximately: {required_vtheta:.2f} m/s")


def plot_trajectory_polar(vr_initial, vtheta_initial):
    """Simulate and plot the trajectory in polar coordinates."""
    r = 1.22e9  # Initial radial distance
    theta = 0  # Initial angle
    vr = vr_initial
    vtheta = vtheta_initial
    time = 0
    trajectory_r = [r]
    trajectory_theta = [theta]

    while time < max_time:
        # RK-like update
        ar = -G * M / r ** 2
        vr_mid = vr + 0.5 * ar * h
        r_mid = r + 0.5 * vr * h
        theta_mid = theta + 0.5 * vtheta * h / r if r != 0 else 0

        vr += (-G * M / r_mid ** 2) * h
        r += vr_mid * h
        theta += vtheta * h / r if r != 0 else 0

        trajectory_r.append(r)
        trajectory_theta.append(theta)

        # Calculate distance to target and check for impact
        distance_to_target = np.sqrt((r * np.cos(theta) - target_r * np.cos(target_theta)) ** 2 +
                                     (r * np.sin(theta) - target_r * np.sin(target_theta)) ** 2)

        if distance_to_target <= 1000:
            print(f"Impact on the asteroid at {time:.2f} seconds.")
            break
        elif r <= R_earth:
            print(f"Impact on Earth at {time:.2f} seconds.")
            break
        time += h

    # Plot in polar coordinates
    plt.figure(figsize=(8, 8))
    plt.polar(trajectory_theta, trajectory_r, label="Projectile Trajectory", color='purple', linewidth=1)
    plt.polar([target_theta], [target_r], 'ro', label="Target Position (Asteroid)")

    # Draw Earth as a circle at the center
    earth_circle = plt.Circle((0, 0), R_earth, transform=plt.gca().transData._b, color='blue', alpha=0.3, label="Earth")
    plt.gca().add_artist(earth_circle)

    plt.title("Projectile Trajectory towards Target (Polar Coordinates)")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper right')
    plt.show()


# Plot the trajectory with optimized initial velocities
plot_trajectory_polar(required_vr, required_vtheta)