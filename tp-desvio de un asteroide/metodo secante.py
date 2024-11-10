import numpy as np


def trajectory(initial_vr, initial_vtheta, initial_r, initial_theta=0, target_r=None, target_theta=0.1, time_step=0.01,
               max_time=10):
    """
    Simulates the trajectory of a projectile given initial radial and tangential velocities
    and returns the final position and distance to target.

    Parameters:
    initial_vr (float): Initial radial velocity of the projectile.
    initial_vtheta (float): Initial tangential (angular) velocity of the projectile.
    initial_r (float): Initial radial distance of the projectile from Earth.
    initial_theta (float): Initial angle of the projectile (default is 0).
    target_r (float): Radial distance of the target.
    target_theta (float): Angle of the target (meteor).
    time_step (float): The time increment for the simulation.
    max_time (float): Maximum simulation time to avoid infinite loops.

    Returns:
    tuple: Final (r, theta) position of the projectile and distance to the target.
    """

    # Constants
    G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
    M = 5.972e24  # Mass of Earth (kg)

    # Initial conditions
    r = initial_r
    theta = initial_theta
    vr = initial_vr  # Initial radial velocity
    vtheta = initial_vtheta  # Initial angular velocity

    # Integration loop
    time = 0
    while time < max_time:
        # Calculate accelerations
        ar = -G * M / r ** 2  # Radial acceleration due to gravity
        atheta = 0  # No external torque is applied, so angular acceleration is zero

        # Update radial and angular velocities using RK2
        vr_mid = vr + 0.5 * ar * time_step
        vtheta_mid = vtheta + 0.5 * atheta * time_step

        r_mid = r + 0.5 * vr * time_step
        theta_mid = theta + 0.5 * vtheta * time_step / r if r != 0 else 0

        ar_mid = -G * M / r_mid ** 2
        vr += ar_mid * time_step
        vtheta += atheta * time_step  # No change in angular acceleration

        # Update positions
        r += vr * time_step
        theta += vtheta * time_step / r if r != 0 else 0

        time += time_step

    # Calculate the Euclidean distance to the target
    if target_r is not None:
        distance_to_target = np.sqrt((r * np.cos(theta) - target_r * np.cos(target_theta)) ** 2 +
                                     (r * np.sin(theta) - target_r * np.sin(target_theta)) ** 2)
    else:
        distance_to_target = None

    return r, theta, distance_to_target


def secant_velocity(target_r, target_theta, initial_guess_vr1, initial_guess_vr2, initial_guess_vtheta, initial_r,
                    tolerance=1e-4, max_iterations=10000, step_limit=1e6):
    """
    Finds the radial velocity required for a projectile to hit a meteor at (target_r, target_theta)
    using the secant method and an initial angular component, with a step limit to control velocity changes.

    Parameters:
    target_r (float): The radial distance of the meteor from Earth.
    target_theta (float): The angular position of the meteor from Earth in radians.
    initial_guess_vr1 (float): First initial guess for the radial velocity.
    initial_guess_vr2 (float): Second initial guess for the radial velocity.
    initial_guess_vtheta (float): Initial tangential (angular) velocity.
    initial_r (float): Initial radial distance of the projectile from Earth.
    tolerance (float): The tolerance for convergence.
    max_iterations (int): Maximum number of iterations.
    step_limit (float): Maximum allowable change in radial velocity per iteration.

    Returns:
    float: The radial velocity required to hit the meteor.
    """

    def distance_to_target(vr):
        # Use the trajectory function to get the distance to the target directly
        final_r, final_theta, distance_to_target = trajectory(vr, initial_guess_vtheta, initial_r, initial_theta=0,
                                                              target_r=target_r, target_theta=target_theta)
        print(
            f"Testing vr={vr}: final position (r, theta) = ({final_r}, {final_theta}), distance to target = {distance_to_target}")
        return distance_to_target

    # Initial guesses for radial velocity and corresponding distances
    vr1, vr2 = initial_guess_vr1, initial_guess_vr2
    f1, f2 = distance_to_target(vr1), distance_to_target(vr2)

    # Secant loop
    for iteration in range(max_iterations):
        if abs(f2 - f1) < 1e-10:  # Prevent division by zero in the secant formula
            print("Secant method failed: Division by zero in secant update.")
            return None

        # Secant method formula
        vr_next = vr2 - f2 * (vr2 - vr1) / (f2 - f1)

        # Limit the step size of vr_next to prevent large jumps
        if abs(vr_next - vr2) > step_limit:
            vr_next = vr2 + np.sign(vr_next - vr2) * step_limit

        # Check for convergence
        if abs(vr_next - vr2) < tolerance:
            print(f"Converged at iteration {iteration}: Required radial velocity = {vr_next}")
            return vr_next

        # Update guesses for the next iteration
        vr1, vr2 = vr2, vr_next
        f1, f2 = f2, distance_to_target(vr2)

        # Debug output to monitor the progression of values
        print(f"Iteration {iteration}: vr1={vr1}, vr2={vr2}, f1={f1}, f2={f2}, vr_next={vr_next}")

    print("Secant method did not converge within the maximum number of iterations.")
    return None


# Define target parameters for the meteor
target_r = 2.24775e9  # Target radial distance (e.g., 10,000 km from Earth)
target_theta = 2.303834613  # Target angular position in radians (e.g., about 5.7 degrees)

# Define initial conditions for the projectile
initial_r = 1.22e9  # Earth's radius in meters (assuming the projectile starts from Earth's surface)
initial_guess_vr1 = 0  # First guess for the initial radial velocity (m/s)
initial_guess_vr2 = 1000  # Second guess for the initial radial velocity (m/s)
initial_guess_vtheta = 790  # Initial guess for the tangential (angular) velocity (m/s)

# Run the secant_velocity function to find the required initial radial velocity
required_radial_velocity = secant_velocity(
    target_r=target_r,
    target_theta=target_theta,
    initial_guess_vr1=initial_guess_vr1,
    initial_guess_vr2=initial_guess_vr2,
    initial_guess_vtheta=initial_guess_vtheta,
    initial_r=initial_r,
    tolerance=1e-4,
    max_iterations=50,
    step_limit=1e6  # Adjust step limit as necessary
)

print(f"The required initial radial velocity to hit the meteor is approximately: {required_radial_velocity} m/s")
print(f"The initial angular (tangential) velocity used was: {initial_guess_vtheta} m/s")
