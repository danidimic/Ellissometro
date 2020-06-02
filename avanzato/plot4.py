import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from sympy import I
from sympy.algebras.quaternion import Quaternion
from stokes import stokes_vector
from sorgente import sorgente
from campione import campione
from interazione import *


#indici rifrzione
n0 = 1
n1 = 1.5+0.00000002j
#n1 = 2.1+0.1j
#Brewster angle
br = np.arctan(n1/n0)
pi = math.pi

#definisco il campione
camp = campione(n0, [n1])

#definisco la sorgente
sorg = sorgente(1.95954488, 1, 0.78539163) #per ora non uso i due argomenti a dx
riniz = stokes_vector()	#inizializzo oggetto vettore di stokes

riniz.generic_polarization(1, 1)	#vettore di Stokes incidente, polarizzazione lineare

#Calcolo quaternione raggio inziale  
s = Quaternion( riniz.I(), riniz.Q()*I, riniz.U()*I, riniz.V()*I )	#quaternione corrispondente al vett Stoke
 		
inter = interazione(0.0001, camp, sorg, riniz)  #inizializzo oggetto interazione

#Liste per risultati quaternioni
Psiq = []
Deltaq = []
Sq = []
Qq = []
Uq = []
Vq = []


nvalues = 100
theta = np.linspace(0, pi/2, nvalues)

for i in range(nvalues):
   
    #calcolo quaternioni del materiale e propago
	inter.materials_to_jones(theta[i])
    
	print('ciao1')
	while inter.nraggi != 0:
		inter.propagazione()
	print('ciao2')

	h_fin = inter.interference()
    
	psi, delta, rfin = grandell([h_fin], s, svfinal=True)

	#Salva i risultati con quaternioni
	Sq.append( rfin.I().real )
	Qq.append( rfin.Q().real )
	Uq.append( rfin.U().real )
	Vq.append( rfin.V().real )
	#Parametri ellissometrici
	Psiq.append( psi )
	Deltaq.append( delta )

#Grafico dei risultati con quaternioni
plt.title("Calcolo tramite quaternioni")
plt.xlabel("Angolo incidente [rad]")
plt.plot(theta, Sq , label='I')
plt.plot(theta, Qq , label='Q')
plt.plot(theta, Uq , label='U')
plt.plot(theta, Vq , label='V')
plt.grid(True)
plt.legend()
plt.show()

Psi = np.dot(Psi, 180/pi)
Psiq = np.dot(Psiq, 180/pi)
Deltaq = np.dot(Deltaq, 180/pi)

#con i quaternioni
plt.title("Calcolo tramite quaternioni")
plt.xlabel("Angolo incidente [rad]")
plt.plot(theta, Psi, label='$\Psi$ Mueller', lw=2.5)
plt.plot(theta, Psiq, label='$\Psi$ quaternioni', linestyle='dashed')
plt.plot(theta, Deltaq, label='$\Delta$ quaternioni')
plt.plot(theta, Delta, label='$\Delta$ Mueller')
plt.grid(True)
plt.legend()
plt.show()



