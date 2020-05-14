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

riniz = stokes_vector()			#inizializzo oggetto vettore di stokes
inter = interazione(0.0001, camp, sorg)  #inizializzo oggetto interazione

#Liste per risultati matrici Mueller
S = []
Q = []
U = []
V = []
Psi = []
Delta = []

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

	inter.materials_to_mueller(theta[i])	#calcolo le matrici di Mueller
	MMrif = inter.muellers.loc[0, 'M_rif']	#matrice di Mueller per la riflessione
	#print("matrice di Mueller: ", MMrif)
	#print()

	#Vettori di Stokes con matrici di Mueller
	riniz.generic_polarization(1, 1)	#vettore di Stokes incidente, polarizzazione lineare
	riniz.mueller_product(MMrif)			#vettore di Stokes riflesso
    
	#Salva i risultati con Matrici Mueller
	S.append( riniz.I() )
	Q.append( riniz.Q() )
	U.append( riniz.U() )
	V.append( riniz.V() )
	Psi.append( riniz.ellipsometric_Psi() )
	Delta.append( riniz.ellipsometric_Delta() )



	#Calcolo con formalismo quaternioni
	riniz.generic_polarization(1, 1)	#vettore di Stokes incidente, polarizzazione lineare  
	s = Quaternion( riniz.I(), riniz.Q()*I, riniz.U()*I, riniz.V()*I )	#quaternione corrispondente al vett Stoke
	hrif = inter.biquaternions.loc[0, 'h_rif']*np.sqrt(MMrif[0,0])				#quaternione corrispondente alla matrice di Mueller	
	
	psi, delta, rfin = grandell([hrif], s, svfinal=True)

	#Salva i risultati con quaternioni
	Sq.append( rfin.I().real )
	Qq.append( rfin.Q().real )
	Uq.append( rfin.U().real )
	Vq.append( rfin.V().real )
	#Parametri ellissometrici
	#Psiq.append( rfin.ellipsometric_Psi().real )
	#Deltaq.append( cmath.phase(shs).real )
	Psiq.append( psi )
	Deltaq.append( delta )


#Grafico dei risultati con matrici di Mueller
plt.title("Calcolo tramite Mueller")
plt.xlabel("Angolo incidente [rad]")
plt.plot(theta, S , label='I')
plt.plot(theta, Q , label='Q')
plt.plot(theta, U , label='U')
plt.plot(theta, V , label='V')
plt.grid(True)
plt.legend()
plt.show()

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
#plt.plot(theta, Delta, label='$\Delta$ Mueller')
plt.grid(True)
plt.legend()
plt.show()



