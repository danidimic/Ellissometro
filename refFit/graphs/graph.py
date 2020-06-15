import math
import numpy as np
import matplotlib.pyplot as plt

#singola riflessione
#filePsi = "Psi_1interfaccia.txt"
#fileDelta = "Delta_1interfaccia.txt"

#layer + substrato
filePsi = "Psi_2interfacce.txt"
fileDelta = "Delta_2interfacce.txt"

theta, psi = np.loadtxt(filePsi, unpack=True)
theta, delta = np.loadtxt(fileDelta, unpack=True)

plt.plot(theta, psi, label='$\Psi$')
plt.plot(theta, delta, label='$\Delta$')
plt.xlabel("Angolo incidente [°]")
plt.ylabel("Grandezze ellissometriche [°]")
#plt.title("RefFit: singola interfaccia")
plt.title("Reffit: layer e substrato")
plt.grid(True)
plt.legend()
plt.show()
