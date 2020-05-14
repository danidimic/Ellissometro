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
sp1 = 0.00002
n2 = 1+0.00000005j
pi = math.pi

#definisco il campione
camp = campione(n0, [n1, sp1, n2])

#definisco la sorgente
sorg = sorgente(1.95954488, 1, 0.78539163) #per ora non uso i due argomenti a dx

r = stokes_vector()			#inizializzo oggetto vettore di stokes
inter = interazione(0.0001, camp, sorg)  #inizializzo oggetto interazione

S = []
Q = []
U = []
V = []
Psi = []
Delta = []


nvalues = 100
theta = np.linspace(0, pi/2, nvalues)

for i in range(nvalues):
	print("iterazione: ", i+1)
	inter.materials_to_mueller(theta[i])	#calcolo le matrici di Mueller

	#Matrice Mueller raggio 1
	MMrif1 = inter.muellers.loc[0, 'M_rif']	
	#Matrici Mueller raggio 2
	MMdw   = inter.muellers.loc[0, 'M_tra_dw']
	MMmat  = inter.muellers.loc[1, 'M_mat']
	MMrif2 = inter.muellers.loc[1, 'M_rif']
	MMup   = inter.muellers.loc[0, 'M_tra_up']
	
	#Raggio 1
	r.generic_polarization(1, 1)	#vettore di Stokes incidente1, polarizzazione lineare 
	s = Quaternion( r.I(), r.Q()*I, r.U()*I, r.V()*I )	#quaternione corrispondente al vett Stoke

	h1rif = inter.biquaternions.loc[0, 'h_rif']*np.sqrt(MMrif1[0,0])				#quaternione corrispondente alla matrice di Mueller	

	h2dw  = inter.biquaternions.loc[0, 'h_tra_dw']*np.sqrt(MMdw[0,0])				#quaternione corrispondente alla matrice di Mueller	
	h2mat = inter.biquaternions.loc[1, 'h_mat']*np.sqrt(MMmat[0,0])					#quaternione corrispondente alla matrice di Mueller	
	h2rif = inter.biquaternions.loc[1, 'h_rif']*np.sqrt(MMrif2[0,0])					#quaternione corrispondente alla matrice di Mueller	
	h2up  = inter.biquaternions.loc[0, 'h_tra_up']*np.sqrt(MMup[0,0])				#quaternione corrispondente alla matrice di Mueller	

	#quaternioni finali corrispondenti ai raggi
	h1, h1daga = multiplication([h1rif])
	h2, h2daga = multiplication([h2up, h2mat, h2rif, h2mat, h2dw])
	htot = h1 + h2

	psi, delta, rfin = grandell([htot], s, svfinal=True)

	S.append( rfin.I() )
	Q.append( rfin.Q() )
	U.append( rfin.U() )
	V.append( rfin.V() )
	Psi.append( psi )
	Delta.append( delta )


Psi = np.dot(Psi, 180/pi)
Delta = np.dot(Delta, 180/pi)

plt.title("Calcolo tramite quaternioni")
plt.xlabel("Angolo incidente [rad]")
plt.plot(theta, S , label='I')
plt.plot(theta, Q , label='Q')
plt.plot(theta, U , label='U')
plt.plot(theta, V , label='V')
plt.grid(True)
plt.legend()
plt.show()

plt.title("Calcolo tramite quaternioni")
plt.xlabel("Angolo incidente [rad]")
plt.plot(theta, Psi, label='$\Psi$ quaternioni', lw=2.5)
plt.plot(theta, Delta, label='$\Delta$ quaternioni')
plt.grid(True)
plt.legend()
plt.show()

