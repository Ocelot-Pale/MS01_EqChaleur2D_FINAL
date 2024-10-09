import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime

current_dateTime = datetime.now()

print("Import des donnees")
solution_simulee = np.genfromtxt("solution_sortie.txt", delimiter=';')
#source = np.genfromtxt("Source.txt", delimiter=';')

print("Affichage de solution_simulee...")

Ly = solution_simulee.shape[1]
Lx = solution_simulee.shape[0]

x = np.linspace(0, 1, Ly)
y = np.linspace(0, 1, Lx)
X, Y = np.meshgrid(x, y)

solution_exacte = np.sin(np.pi/2 * X) * np.cos(np.pi * Y)

# Calculate L2 error
L2_error = np.sqrt(np.nansum((solution_simulee - solution_exacte)**2))/np.sqrt(np.nansum((solution_exacte)**2))
#print(solution_simulee)
# Determine the color scale limits based on the min and max values from both datasets
vmax = max(np.nanmax(solution_simulee) , np.nanmax(solution_exacte))
vmin = min(np.nanmin(solution_simulee) , np.nanmin(solution_exacte))
print(vmin)
print(vmax)
fig, axs = plt.subplots(1, 2, figsize=(10, 5))

contour1 = axs[0].contourf(X, Y,  solution_simulee,np.linspace(-1,1,20) , cmap='hot', vmin=vmin, vmax=vmax)
axs[0].set_xlim(0, 1)
axs[0].set_ylim(0, 1)
axs[0].set_title("Solution Simulee")
fig.colorbar(contour1, ax=axs[0])

contour2 = axs[1].contourf(X, Y, solution_exacte,np.linspace(-1,1,20) , cmap='hot', vmin=vmin, vmax=vmax)
axs[1].set_xlim(0, 1)
axs[1].set_ylim(0, 1)
axs[1].set_title("Solution Exacte")
fig.colorbar(contour2, ax=axs[1])

# Add L2 error text on the figure

print("Sauvegarde Figure...")
#plt.savefig("solution_simulee_et_exactee"+str(current_dateTime)+".png", bbox_inches='tight')
print(L2_error)
print("FIN")
