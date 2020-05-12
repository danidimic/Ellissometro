import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from sympy import I
from sympy.algebras.quaternion import Quaternion
from stokes import stokes_vector
from sorgente import sorgente
from campione import campione
from interazione import interazione

#prodotto scalare fra biquaterioni, per ora inutile
def scalar_prod(q1, q2):
	return q1.a*q2.a + q1.b*q2.b + q1.c*q2.c + q1.d*q2.d


#coniugato di un quaternione
def conjugate(h):
	tau   =  np.conj( h.a )
	alfa  =  np.conj(-h.b*I)
	beta  =  np.conj(-h.c*I)
	gamma =  np.conj(-h.d*I)
	
	return Quaternion( tau, alfa*I, beta*I, gamma*I )


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
	h = inter.biquaternions.loc[0, 'h_rif']				#quaternione corrispondente alla matrice di Mueller	
	
	hdaga = conjugate(h)
	shdaga = s.mul(hdaga)
	sfin	= h.mul(shdaga)
	rfin = stokes_vector( complex(sfin.a), complex(sfin.b*(-1j)), complex(sfin.c*(-1j)), complex(sfin.d*(-1j)) )  #vettore di Stokes finale

	hs = h.mul(s)			#prodotto hs tra quaternioni
	shs = scalar_prod(h, hs)*MMrif[0,0]	#prodotto scalare s.hs

	#Salva i risultati con quaternioni
	Sq.append( rfin.I() )
	Qq.append( rfin.Q() )
	Uq.append( rfin.U() )
	Vq.append( rfin.V() )
	#Parametri ellissometrici
	Psiq.append( rfin.ellipsometric_Psi() )
	Deltaq.append( cmath.phase(shs) )


#Grafico dei risultati con matrici di Mueller
plt.title("Calcolo tramite Mueller")
plt.xlabel("Angolo incidente [rad]")
plt.plot(theta, S , label='I')
plt.plot(theta, Q , label='Q')
plt.plot(theta, U , label='U')
plt.plot(theta, V , label='V')
#plt.plot(theta, Psi, label='$\Psi$', lw=2.5)
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
plt.plot(theta, Psiq, label='$\Psi$ quaternioni', lw=2.5)
plt.plot(theta, Psiq, label='$\Psi$ Mueller', linestyle='dashed')
plt.plot(theta, Deltaq, label='$\Delta$ quaternioni')
#plt.plot(theta, Delta, label='$\Delta$ Mueller')
plt.grid(True)
plt.legend()
plt.show()




