import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from grandell import ellissometro

passo = 0.01
nvalues = int(math.pi/2/passo)

r_pi = []
r_sigma = []
delta1 = []
psi1 = []
Theta = []
E = ellissometro(0.0001)

I = []
Q = []
U = []
V = []

for i in range(nvalues):
	theta = passo*i  #angoli theta iniziali
	res = E.grandell(theta) #calcolo delle grandezze ellissometriche

	r_pi.append   (res[0])
	r_sigma.append(res[1])
	delta1.append (res[2])
	psi1.append   (res[3])
	Theta.append  (theta )
	
	#conversione da radianti in gradi
	delta1[i] = 180*delta1[i]/math.pi
	psi1[i] = 180*psi1[i]/math.pi
	#calcolo dei moduli quadri
	#r_pi[i] = abs(r_pi[i])**2
	#r_sigma[i] = abs(r_sigma[i])**2
	r_pi[i] = r_pi[i]*math.tan(E.sorg_val.psi_0)*cmath.exp(1j*E.sorg_val.delta_0)
    
	I.append((abs(r_sigma[i]))**2+(abs(r_pi[i]))**2)
	Q.append((abs(r_sigma[i]))**2-(abs(r_pi[i]))**2)
	U.append(2*abs(r_sigma[i])*abs(r_pi[i])*np.cos(res[2]))
	V.append(2*abs(r_sigma[i])*abs(r_pi[i])*np.sin(res[2]))
    
Theta = np.dot(Theta, 180/math.pi)    

#Grafico di Delta e Psi
plt.xlabel("angolo incidente [°]")
plt.ylabel("$\Delta$ | $\Psi$ [°]")
plt.plot(Theta, delta1, label='$\Delta$') #grafico delta
plt.plot(Theta, psi1, label='$\Psi$')	#grafico psi
plt.legend()
plt.grid(True)
plt.show()

#Grafico di r_pi e r_sigma
plt.xlabel("angolo incidente [°]")
plt.ylabel("$R_{\pi}$ | $R_{\sigma}$ [°]")
plt.plot(Theta, r_pi, label='$R_{\pi}$')
plt.plot(Theta, r_sigma, label='$R_{\sigma}$')
plt.legend()
plt.grid(True)
plt.show()

#Grafico di I, Q, U, V
plt.title("Metodo standard")
plt.xlabel("angolo incidente [°]")
plt.ylabel("$I$ | $Q$ [°]")
plt.plot(Theta, I, label='$I$')
plt.plot(Theta, Q, label='$Q$')
plt.plot(Theta, U, label='$U$')
plt.plot(Theta, V, label='$V$')
plt.legend()
plt.grid(True)
plt.show()