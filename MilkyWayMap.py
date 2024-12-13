import numpy as np
import matplotlib.pyplot as plt

R0 = 8.5  # kpc
V0 = 220  # km/s
theta_values = np.linspace(0, 2 * np.pi, 1000)  # Angles in radians

def calculate_R(Vr, l):
    l_rad = np.radians(l)
    return R0 * V0 * np.sin(l_rad) / (V0 * np.sin(l_rad) + Vr)

def calculate_r(R, l):
    l_rad = np.radians(l)
    term = R**2 - R0**2 * np.sin(l_rad)**2
    if term < 0:
        return None, None  # Discard unphysical solutions
    sqrt_term = np.sqrt(term)
    r_plus = R0 * np.cos(l_rad) + sqrt_term
    r_minus = R0 * np.cos(l_rad) - sqrt_term
    return r_plus, r_minus

# Function to calculate Cartesian coordinates
def calculate_xy(r, l):
    theta = np.radians(l - 90)
    x = r * np.cos(theta)
    y = r * np.sin(theta) + R0  # Adjust y to center on Galactic center
    return x, y

observed_data = [
    (39.82, 63.74), (39.82, 30.93), (39.82, 6.45), (39.82, -45.51),  # First observation
    (42.36, 59.60), (42.36, 26.84), (42.36, 5.55), (42.36, -47.19),  # Second observation
    (44.95, 55.97), (44.95, 22.62), (44.95, 5.52), (44.95, -49.57),  # Third observation
    (47.45, 50.63), (47.45, 12.67), (47.45, -50.89),                # Fourth observation
    (49.78, 46.96), (49.78, 11.08), (49.78, -51.87),                # Fifth observation
    (52.46, 37.52), (52.46, 7.28), (52.46, -54.15),                 # Sixth observation
    (54.94, 18.56), (54.94, -60.20),                               # Seventh observation
    (57.37, 15.87), (57.37, -62.65),                               # Eighth observation
    (59.15, 14.17), (59.15, -64.58),                               # Ninth observation
    (62.42, 11.81), (62.42, -66.23),                               # Tenth observation
    (64.89, 9.25), (64.89, -66.77),                                # Eleventh observation
    (67.38, 6.72), (67.38, -66.87),                                # Twelfth observation
    (69.81, 4.65), (69.81, -67.90),                                # Thirteenth observation
    (72.41, 3.07), (72.41, -70.34),                                # Fourteenth observation
    (74.97, 2.47), (74.97, -73.37),                                # Fifteenth observation
    (77.34, 3.03), (77.34, -62.42),                                # Sixteenth observation
    (79.79, 2.34), (79.79, -40.67), (79.79, -73.86),               # Seventeenth observation
    (82.33, 0.66), (82.33, -42.58), (82.33, -76.09),               # Eighteenth observation
    (84.81, -1.10), (84.81, -43.35), (84.81, -77.37),              # Nineteenth observation
    (87.37, -0.32), (87.37, -32.44), (87.37, -77.97),              # Twentieth observation
    (89.75, 4.26), (89.75, -7.76), (89.75, -54.11)                 # Twenty-first observation
]

x_coords = []
y_coords = []

for l, Vr in observed_data:
    R = calculate_R(Vr, l)
    if R is None:
        continue
    r_plus, r_minus = calculate_r(R, l)
    if r_plus:
        x, y = calculate_xy(r_plus, l)
        x_coords.append(x)
        y_coords.append(y)
    if r_minus:
        x, y = calculate_xy(r_minus, l)
        x_coords.append(x)
        y_coords.append(y)

plt.style.use('classic')
# Spiral arm function
def spiral_arm(theta, R0, k, theta0):
    return R0 * np.exp(k * (-theta - theta0))

# Parameters for the arms
k_orion = -0.1
theta0_orion = np.radians(20)
k_perseus = -0.23
theta0_perseus = np.radians(40)

# Compute x, y for the arms
R_orion = spiral_arm(theta_values, R0, k_orion, theta0_orion)
x_orion = R_orion * np.cos(theta_values)
y_orion = R_orion * np.sin(theta_values)

R_perseus = spiral_arm(theta_values, R0, k_perseus, theta0_perseus)
x_perseus = R_perseus * np.cos(theta_values)
y_perseus = R_perseus * np.sin(theta_values)

# Function to find the cutoff for the third quadrant (x < 0 and y < 0)
def find_cutoff(x, y):
    # Find the index where both x and y become negative
    cutoff_idx = np.where((x < 0) & (y < 0))[0]
    if len(cutoff_idx) > 0:
        return cutoff_idx[0]  # return the first index where this happens
    return len(x)  # If never negative, return full length

# Find the cutoff indices for both arms
cutoff_idx_orion = find_cutoff(x_orion, y_orion)
cutoff_idx_perseus = find_cutoff(x_perseus, y_perseus)

# Limit the data to stop plotting when in the third quadrant
x_orion_limited = x_orion[:cutoff_idx_orion]
y_orion_limited = y_orion[:cutoff_idx_orion]
x_perseus_limited = x_perseus[:cutoff_idx_perseus]
y_perseus_limited = y_perseus[:cutoff_idx_perseus]

# Plotting the map with arms
plt.figure(figsize=(10, 10))
plt.scatter(x_coords, y_coords, c='blue', label='HI clouds')  # Observed points
plt.plot(x_orion_limited, y_orion_limited, color='orange', label='Orion-Cygnus Arm')
plt.plot(x_perseus_limited, y_perseus_limited, color='green', label='Perseus Arm')

# Add Sun and Galactic center
plt.scatter(0, 0, color='red', label='Galactic Center', s=100)
plt.scatter(0, R0, color='yellow', label='Sun', s=100)

# Adjust plot
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('X (kpc)')
plt.ylabel('Y (kpc)')
plt.xlim(-15, 15)
plt.ylim(-15, 15)
plt.title('Map of the Milky Way with Spiral Arms')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()