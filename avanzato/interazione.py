import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from mueller import *
from stokes import *

n = 100
pi=math.pi
theta = np.linspace(0, pi/2, n)

v = stokes_vector()
v.linear_polarization(1, pi/4)

Psi = []
Delta = []

for i in range(n):
	m = MM(1+0j, 1.5+0.00000002j, theta[i], 0.5, 0.00002)
	mat_ref = m.mueller_reflection()
	mat_tra = m.mueller_transmission()
	mat_lay = m.mueller_layer()
	
	v.layer_interaction(mat_ref, mat_tra, mat_lay)

	psi, delta = v.ellipsometric_parameters()

	Psi.append(psi)
	Delta.append(delta)	


plt.plot(theta, Psi, label="$\Psi$")
plt.plot(theta, Delta, label="$\Delta$")
plt.legend()
plt.grid(True)
plt.show()
