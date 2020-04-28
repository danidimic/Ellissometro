import math
import numpy as np
import matplotlib.pyplot as plt
from grandell import ellissometro

passo = 0.01
#print(passo)
nvalues = int(math.pi/2/passo)

r_pi = []
r_sigma = []
delta1 = []
psi1 = []

for i in range(nvalues):
	theta = passo*i  #angoli theta iniziali
	E = ellissometro(0.0001, theta)
	res = E.grandell() #calcolo delle grandezze ellissometriche

	r_pi.append   (res[0])
	r_sigma.append(res[1])
	delta1.append (res[2])
	psi1.append   (res[3])
	
	#r_pi[i] = 180*r_pi[i]/math.pi
	#r_sigma[i] = 180*r_sigma[i]/math.pi
	delta1[i] = 180*delta1[i]/math.pi
	psi1[i] = 180*psi1[i]/math.pi
	
	r_pi[i] = abs(r_pi[i])**2
	r_sigma[i] = abs(r_sigma[i])**2
	


theta = np.linspace(0, 90, nvalues)

plt.xlabel("angolo incidente [°]")
plt.ylabel("$\Delta$ | $\Psi$")
plt.plot(theta, delta1, label='$\Delta$') #grafico delta
plt.plot(theta, psi1, label='$\Psi$')	#grafico psi
plt.legend()
plt.grid(True)
plt.show()

plt.xlabel("angolo incidente [°]")
plt.ylabel("$R_{\pi}$ | $R_{\sigma}$")
plt.plot(theta, r_pi, label='$R_{\pi}$')
plt.plot(theta, r_sigma, label='$R_{\sigma}$')
plt.legend()
plt.grid(True)
plt.show()
