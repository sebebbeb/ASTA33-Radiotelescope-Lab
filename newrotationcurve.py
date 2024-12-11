import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Given data
l_values = [59.8, 69.8, 79.8, 89.9, 99.8, 110.0, 119.8, 129.8, 139.8, 149.8, 159.8, 169.9, 179.8, 189.8, 199.9, 209.8, 219.8, 229.9, 239.9]
V_r_max_values = [54.6, 36.8, 35.1, 23.3, 20.3, 16.9, 17.3, 14.3, 16.4, 10.0, 35.0, 26.3, 28.0, 36.0, 61.0, 59.3, 91.4, 88.1, 113.2]

# Constants
R_0 = 8.2  # kpc
V_0 = 220  # km/s

# Convert l to radians
l_radians = np.radians(l_values)

# Initialize lists for results
R_values = []
V_values = []

# Loop through data points to calculate R and V
for l, V_r in zip(l_radians, V_r_max_values):
    if 0 <= np.degrees(l) <= 180:
        # Quadrants I and II
        R = R_0 * np.sin(l)
        V = V_r + V_0 * np.sin(l)
    else:
        # Quadrants III and IV
        R = R_0 * np.sin(l)
        V = V_r - V_0 * np.sin(l)

    R_values.append(abs(R))
    V_values.append(V)

# Create a DataFrame for the results
results = pd.DataFrame({
    'l (degrees)': l_values,
    'V_r,max': V_r_max_values,
    'R (kpc)': R_values,
    'V (km/s)': V_values
})

print(results)

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(R_values, V_values, marker='o', linestyle='-', color='b', label='V vs R')
plt.xlabel('R (kpc)')
plt.ylabel('V (km/s)')
plt.title('Plot of V vs R for Galactic Longitudes')
plt.legend()
plt.grid()
plt.show()
