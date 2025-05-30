import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Gravitational constant (m³/kg·s²)
M = 1.989e30     # Mass of Sun (kg)
AU = 1.496e11    # Astronomical unit (m)

# Function to calculate center of mass position
def rcom(s, m):
    N = len(m)
    r = np.zeros(3)
    total_mass = sum(m)
    for i in range(N):
        r += m[i] * np.array([s[i, 0], s[i, 2], s[i, 4]])
    return r / total_mass

# Function to calculate center of mass velocity
def vcom(s, m):
    N = len(m)
    v = np.zeros(3)
    total_mass = sum(m)
    for i in range(N):
        v += m[i] * np.array([s[i, 1], s[i, 3], s[i, 5]])
    return v / total_mass

# Function defining the system of ODEs
def f(s, t, m):
    N = len(m)
    ode = np.zeros((N, 6))

    for i in range(N):
        # Set velocity derivatives
        ode[i, 0] = s[i, 1]
        ode[i, 2] = s[i, 3]
        ode[i, 4] = s[i, 5]

        # Compute gravitational forces
        for j in range(N):
            if i != j:
                diff = np.array([s[i, 0] - s[j, 0], s[i, 2] - s[j, 2], s[i, 4] - s[j, 4]])
                radius = np.linalg.norm(diff)

                # Acceleration due to gravity
                ode[i, 1] += -G * m[j] * diff[0] / radius**3
                ode[i, 3] += -G * m[j] * diff[1] / radius**3
                ode[i, 5] += -G * m[j] * diff[2] / radius**3

    return ode

# RK4 Stepper
def rk4_step(s, t, h, f):
    k1 = h * f(s, t, m)
    k2 = h * f(s + 0.5 * k1, t + 0.5 * h, m)
    k3 = h * f(s + 0.5 * k2, t + 0.5 * h, m)
    k4 = h * f(s + k3, t + h, m)
    return s + (k1 + 2 * k2 + 2 * k3 + k4) / 6

# Multi-Variable ODE Solver
def MVODE(t_end, t0, num_steps, initial, f, m):
    h = (t_end - t0) / (num_steps - 1)
    t_values = [t0]
    s_values = [initial]

    while t_values[-1] < t_end:
        t_current = t_values[-1]
        s_current = s_values[-1]

        s_next = rk4_step(s_current, t_current, h, f)
        s_values.append(s_next)
        t_next = t_current + h
        t_values.append(t_next)

    return np.array(t_values), np.array(s_values)

# Simulation parameters
t_end = 30 * 365.25 * 24 * 3600  # 30 years in seconds
num_steps = 10000

# Initial conditions (semi-major axis and masses)
a1 = 5.204 * AU  # Jupiter
m1 = 1.898e27
v1 = np.sqrt(G * M / a1)

a2 = 9.583 * AU  # Saturn
m2 = 5.683e26
v2 = np.sqrt(G * M / a2)

m = [m1, m2, M]  # Masses

# Initial positions and velocities
initial = np.array([
    [a1, 0, 0, v1, 0, 0],  # Jupiter
    [a2, 0, 0, v2, 0, 0],  # Saturn
    [0, 0, 0, 0, 0, 0]     # Sun
])

# Compute initial center of mass position and velocity
r_cm_initial = rcom(initial, m)
v_cm_initial = vcom(initial, m)

# Shift positions and velocities to center-of-mass frame
for i in range(len(m)):
    initial[i, 0] -= r_cm_initial[0]
    initial[i, 2] -= r_cm_initial[1]
    initial[i, 4] -= r_cm_initial[2]

    initial[i, 1] -= v_cm_initial[0]
    initial[i, 3] -= v_cm_initial[1]
    initial[i, 5] -= v_cm_initial[2]

# Run the simulation
t_arr, S_arr = MVODE(t_end, 0, num_steps, initial, f, m)

# Plot results
plt.figure(figsize=(8, 8))
planet_labels = ["Jupiter", "Saturn", "Sun"]

for i in range(len(m)):
    plt.plot(S_arr[:, i, 0] / 1e9, S_arr[:, i, 2] / 1e9, label=planet_labels[i])

plt.xlabel("x position (billion m)", fontsize=18)
plt.ylabel("y position (billion m)", fontsize=18)
plt.title("Jupiter and Saturn Orbiting the Sun (Center-of-Mass RK4 Method)", fontsize=18)
plt.legend(fontsize=16, loc='upper right')
plt.axis('equal')
plt.grid()
plt.show()




