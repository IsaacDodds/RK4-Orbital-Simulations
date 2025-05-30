import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M = 1.989e30  # Mass of the Sun (kg)


def f(s, t):
    """Defines the system of differential equations for Halley's Comet orbit."""
    x, v_x, y, v_y = s
    r = np.sqrt(x**2 + y**2)  # Radial distance from the Sun
    a_x = -G * M * x / r**3  # Acceleration in x direction
    a_y = -G * M * y / r**3  # Acceleration in y direction
    return np.array([v_x, a_x, v_y, a_y])

def rk4_step(s, t, h, f):
    """Performs a single step of the 4th-order Runge-Kutta method."""
    k1 = h * f(s, t)
    k2 = h * f(s + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(s + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(s + k3, t + h)
    return s + (k1 + 2*k2 + 2*k3 + k4) / 6

def Solver(t_end, t0, N, initial, f):
    """Solves the differential equations using the RK4 method over the given time span."""
    h = (t_end - t0) / (N - 1)  # Time step size
    t_values = [t0]
    s_values = [initial]
    
    while t_values[-1] < t_end:
        t_current = t_values[-1]
        s_current = s_values[-1]    
        s_next = rk4_step(s_current, t_current, h, f)
        s_values.append(s_next)
        t_next = t_current + h 
        t_values.append(t_next)
    
    return np.array(t_values), np.array(s_values).T

N = 50000  # Number of time steps
initial = np.array([5.2e12, 0, 0, 880])  # Initial conditions [x, v_x, y, v_y]
t0 = 0  # Initial time (s)
t_end = 76 * 365.25 * 24 * 60 * 60  # Total simulation time (seconds)
# Solver was run
t_arr, S_arr = Solver(t_end, t0, N, initial, f)

# Plot of orbit
plt.figure(figsize=(8, 8))
plt.plot(S_arr[0, :] / 1e9, S_arr[2, :] / 1e9, 'bo', markersize=0.2, label="Halley's Comet")
plt.plot(0, 0, 'yo', markersize=10, label="Sun")
plt.xlabel("x position (billion m)", fontsize=18)
plt.ylabel("y position (billion m)", fontsize=18)
plt.title("Orbit of Halley's Comet (RK4 Method)", fontsize=18)
plt.legend( fontsize=16)
plt.axis('equal')
plt.grid()
plt.show()


