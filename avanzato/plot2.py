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

	#Matrice Mueller interfaccia 1
	MMrif1dw = inter.muellers.loc[0, 'M_rif_dw']
	MMrif1up = inter.muellers.loc[0, 'M_rif_up']
	MMtra1up   = inter.muellers.loc[0, 'M_tra_up']
	MMtra1dw   = inter.muellers.loc[0, 'M_tra_dw']
	#Matrici Mueller interfaccia 2
	MMmat  = inter.muellers.loc[1, 'M_mat']
	MMrif2dw = inter.muellers.loc[1, 'M_rif_dw']
	MMrif2up = inter.muellers.loc[1, 'M_rif_up']

	#Quaternioni interfaccia 1
	hrif1dw = inter.muellers.loc[0, 'h_rif_dw']
	hrif1up = inter.muellers.loc[0, 'h_rif_up']
	htra1up   = inter.muellers.loc[0, 'h_tra_up']
	htra1dw   = inter.muellers.loc[0, 'h_tra_dw']
	#Quaternioni interfaccia 2
	hmat  = inter.muellers.loc[1, 'h_mat']
	hrif2dw = inter.muellers.loc[1, 'h_rif_dw']
	hrif2up = inter.muellers.loc[1, 'h_rif_up']
	
	'''
	h1rifdw = inter.biquaternions.loc[0, 'h_rif_dw']*np.sqrt(MMrif1[0,0])				#quaternione corrispondente alla matrice di Mueller	

	h2dw  = inter.biquaternions.loc[0, 'h_tra_dw']*np.sqrt(MMdw[0,0])				#quaternione corrispondente alla matrice di Mueller	
	h2mat = inter.biquaternions.loc[1, 'h_mat']*np.sqrt(MMmat[0,0])					#quaternione corrispondente alla matrice di Mueller	
	h2rif = inter.biquaternions.loc[1, 'h_rif']*np.sqrt(MMrif2[0,0])					#quaternione corrispondente alla matrice di Mueller	
	h2up  = inter.biquaternions.loc[0, 'h_tra_up']*np.sqrt(MMup[0,0])				#quaternione corrispondente alla matrice di Mueller	
	'''

	#Raggio 1
	r.generic_polarization(1, 1)	#vettore di Stokes incidente1, polarizzazione lineare 
	s = Quaternion( r.I(), r.Q()*I, r.U()*I, r.V()*I )	#quaternione corrispondente al vett Stoke

	#quaternioni finali corrispondenti ai raggi
	h1, h1daga = multiplication([hrif1dw])
	h2, h2daga = multiplication([htra1up, hmat, hrif2dw, hmat, htra1dw])
	h3, h3daga = multiplication([htra1up, hmat, hrif2dw, hmat, hrif1up, hmat, hrif2dw, hmat, htra1dw])

	htot = h1 + h2 + h3

	#psi, delta, rfin = grandell([h1rif], s, svfinal=True)
	#psi, delta, rfin = grandell([h2up, h2mat, h2rif, h2mat, h2dw], s, svfinal=True)
	psi, delta, rfin = grandell([htot], s, svfinal=True)

	S.append( rfin.I().real )
	Q.append( rfin.Q().real )
	U.append( rfin.U().real )
	V.append( rfin.V().real )
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

