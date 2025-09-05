import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# === Load data from file ===
data = np.loadtxt('pyramid.dat', delimiter=',')

# === Generate axes ===
time_steps = np.arange(data.shape[0]) * 5  # "time (5 day time steps)"
ages = np.arange(data.shape[1])            # "age (years)"

T, A = np.meshgrid(time_steps, ages, indexing='ij')  # Create meshgrid

# === Create 3D plot ===
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# === Plot surface ===
surf = ax.plot_surface(T, A, data, cmap='viridis', edgecolor='none')

# === Labels and title ===
ax.set_xlabel('Time (5 day time steps)')
ax.set_ylabel('Age (years)')
ax.set_zlabel('Number of humans')
ax.set_title('3D Plot of Human Count by Time and Age')

# === Color bar ===
fig.colorbar(surf, shrink=0.5, aspect=10, label='Number of humans')

plt.tight_layout()
plt.show()

