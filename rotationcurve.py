import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

'''
OLD VALUES
l_values = [59.8, 69.8, 79.8, 89.9, 99.8, 110.0, 119.8, 129.8, 139.8, 149.8, 159.8, 169.9, 179.8, 189.8, 199.9, 209.8, 219.8, 229.9, 239.9]
V_r_max_values = [54.6, 36.8, 35.1, 23.3, 20.3, 16.9, 17.3, 14.3, 16.4, 10.0, 35.0, 26.3, 28.0, 36.0, 61.0, 59.3, 91.4, 88.1, 113.2]
l_values = [59.8, 69.8, 79.8, 89.9] # only <90 degree values
V_r_max_values = [54.6, 36.8, 35.1, 23.3] # only <90 degree values
OLD VALUES
'''

l_values = [39.8, 42.4, 45.0, 47.6, 49.7, 52.5, 55.0, 57.4, 59.1, 62.4, 64.9, 67.4, 69.8, 72.4, 75.0, 77.3, 79.8, 82.3, 84.8, 87.4, 89.8]
V_r_max_values = [91.9, 89.9, 87.4, 90.1, 86.6, 84.0, 80.9, 71.4, 62.1, 54.1, 46.5, 42.0, 37.1, 35.1, 33.3, 34.6, 36.4, 35.8, 33.3, 27.3, 16.1]
V_r_errors = [1.830, 1.958, 2.030, 2.518, 2.52, 2.947, 3.930, 3.508, 3.039, 2.691, 2.379, 2.261, 2.088, 2.063, 1.989, 2.036, 2.184, 2.251, 2.206, 1.800, 0.893] 

R_0 = 8.5  # kpc
V_0 = 220  # km/s
G = 4.300e-6  # Gravitational constant in kpc km^2 s^-2 M_sun^-1

l_radians = np.radians(l_values)

R_values = R_0 * np.sin(l_radians)
V_values = [V_r + V_0 * np.sin(l) for V_r, l in zip(V_r_max_values, l_radians)]
M_enclosed = [(V**2 * R) / G for V, R in zip(V_values, R_values)]

M_central = 5e10  # Solar masses

V_keplerian_values = [np.sqrt(G * M_central / R) for R in R_values]


sigma_l_deg = 0.2  
sigma_l_rad = np.radians(sigma_l_deg)  

R_errors = R_0 * np.cos(l_radians) * sigma_l_rad
V_errors = [
    np.sqrt(V_r_err**2 + (V_0 * np.cos(l) * sigma_l_rad)**2)
    for V_r_err, l in zip(V_r_errors, l_radians)
]

M_enclosed_errors = [
    M_enclosed[i] * np.sqrt((2 * V_errors[i] / V_values[i])**2 + (R_errors[i] / R_values[i])**2)
    for i in range(len(V_values))
]

results = pd.DataFrame({
    'l (degrees)': l_values,
    'V_r,max (km/s)': V_r_max_values,
    'V_r,max error (km/s)': V_r_errors,
    'R (kpc)': R_values,
    'R error (kpc)': R_errors,
    'V (km/s)': V_values,
    'V error (km/s)': V_errors,
    'V Keplerian (km/s)': V_keplerian_values,
    'M_enclosed (M_sun)': M_enclosed,
    'M_enclosed error (M_sun)': M_enclosed_errors
})

print(results)

# Rotation curve
plt.style.use('classic')
plt.figure(figsize=(8, 6))
plt.errorbar(R_values, V_values, xerr=R_errors, yerr=V_errors, fmt='o', color='b', label='V vs R')
plt.plot(R_values, V_values, linestyle='-', color='b')
plt.plot(R_values, V_keplerian_values, linestyle='--', color='r', label='V (Keplerian) vs R')
plt.xlabel('R (kpc)')
plt.ylabel('V (km/s)')
plt.ylim(0, 300)
plt.xlim(5.3, 8.6)
plt.title('Plot of V vs R')
plt.legend(loc='best')
plt.grid()
plt.show()

# Enclosed Mass
plt.figure(figsize=(8, 6))
plt.errorbar(R_values, M_enclosed, xerr=R_errors, yerr=M_enclosed_errors, fmt='o', color='g', label='Enclosed Mass (Observed)')
plt.xlabel('R (kpc)')
plt.ylabel('Enclosed Mass ($M_{\odot}$)')
plt.title('Enclosed Mass as a Function of Radius')
plt.legend(loc='best')
plt.grid()
plt.show()