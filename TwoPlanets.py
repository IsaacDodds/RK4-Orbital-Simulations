import numpy as np
import matplotlib.pyplot as plt

G = 6.67430e-11  # Gravitational constant
M = 1.989e30  # Mass of the Sun
AU = 1.496e11  # Astronomical unit


def f(s, t, m):
    N = len(m)
    ode = np.zeros((N, 6))

    for i in range(N):
        # Set velocity derivatives
        ode[i, 0] = s[i, 1]
        ode[i, 2] = s[i, 3]
        ode[i, 4] = s[i, 5]

        r = np.sqrt(s[i, 0]**2 + s[i, 2]**2 + s[i, 4]**2)
        ode[i, 1] = -G * M * s[i, 0] / r**3
        ode[i, 3] = -G * M * s[i, 2] / r**3
        ode[i, 5] = -G * M * s[i, 4] / r**3

        # Compute gravitational forces between bodies
        for j in range(N):
            if i != j:
                diff = np.array([s[i, 0] - s[j, 0], s[i, 2] - s[j, 2], s[i, 4] - s[j, 4]])
                radius=np.sqrt((s[i,0]-s[j,0])**2+(s[i,2]-s[j,2])**2+(s[i,4]-s[j,4])**2)
                
                factor = -G * m[j] / radius**3
                ode[i, 1] += factor * diff[0]
                ode[i, 3] += factor * diff[1]
                ode[i, 5] += factor * diff[2]
    
    return ode

def rk4_step(s, t, h, f):
    k1 = h * f(s, t, m)
    k2 = h * f(s + 0.5 * k1, t + 0.5 * h, m)
    k3 = h * f(s + 0.5 * k2, t + 0.5 * h, m)
    k4 = h * f(s + k3, t + h, m)
    return s + (k1 + 2*k2 + 2*k3 + k4) / 6

def Solver(t_end, t0, N, initial, f, m):
    h = (t_end - t0) / (N - 1)
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

a1=2.52*AU
m1=1e-3*M
v1=np.sqrt(G*M/a1)

a2=5.24*AU
m2=4e-2*M
v2=np.sqrt(G*M/a2)

initial = np.array([
    [a1, 0, 0, v1, 0, 0],  # Jupiter
    [a2, 0, 0, v2, 0, 0]   # Saturn
])
m = [m1, m2]
t_arr, S_arr = Solver(5e8, 0, 10000, initial, f, m)


plt.figure(figsize=(8, 8))
for i in range(len(m)):
    plt.plot(S_arr[:, i, 0] / 1e9, S_arr[:, i, 2] / 1e9, label=f'Body {i}')

plt.plot(0, 0, 'yo', markersize=10, label="Sun")
plt.xlabel("x position (billion m)", fontsize=18)
plt.ylabel("y position (billion m)", fontsize=18)
plt.title("2 planets orbiting the sun (RK4 Method)", fontsize=18)
plt.legend(fontsize=16, loc='upper right')
plt.axis('equal')
plt.grid()
plt.show()

a1 = 5.204 * AU
m1 = 1.898e27
v1 = np.sqrt(G * M / a1)

a2 = 9.583 * AU
m2 = 5.683e26
v2 = np.sqrt(G * M / a2)

# Initial positions and velocities (x, vx, y, vy, z, vz)
initial = np.array([
    [a1, 0, 0, v1, 0, 0],  # Jupiter
    [a2, 0, 0, v2, 0, 0]   # Saturn
])
m = [m1, m2]
t_arr, S_arr = Solver(30 * 365.25 * 24 * 3600, 0, 10000, initial, f, m)

plt.figure(figsize=(8, 8))
planet_labels = ["Jupiter", "Saturn"]

for i in range(len(m)):
    plt.plot(S_arr[:, i, 0] / 1e9, S_arr[:, i, 2] / 1e9, label=planet_labels[i])

plt.plot(0, 0, 'yo', markersize=10, label="Sun")
plt.xlabel("x position (billion m)", fontsize=18)
plt.ylabel("y position (billion m)", fontsize=18)
plt.title("Jupiter and Saturn Orbiting the Sun (RK4 Method)", fontsize=18)
plt.legend(fontsize=16, loc='upper right')
plt.axis('equal')
plt.grid()
plt.show()
