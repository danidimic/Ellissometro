from interazione.py import interazione
import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from stokes import stokes_vector
from sorgente import sorgente
from campione import campione
from matrix import *
from sympy.algebras.quaternion import Quaternion
from sympy import I

my_sorgente = sorgente(1.95954488, 1, 0.78539163)

#indici rifrzione
n0 = 1
n1 = 2.1+1j

#Brewster angle
br = np.arctan(n1/n0)

#definisco il raggio incidente
raggio_iniz = stokes_vector()
raggio_iniz.generic_polarization(1,1.5, 0.12, 0.34)

#definisco il campione
my_campione = campione(n0, [n1])

#prodotto scalare fra biquaterioni, per ora inutile
def scalar_prod(q1, q2):
	return q1.a*q2.a + q1.b*q2.b + q1.c*q2.c + q1.d*q2.d

passo = 0.1
nvalues = int(math.pi/2/passo)

raggio_finale_M = []
raggio_finale_h = []

for i in range(nvalues):
    
	theta = passo*i  #angoli theta iniziali
	my_interazione = interazione(0.0001, my_campione,my_sorgente, theta)
    
    #con matrici di Mueller
	my_interazione.materials_to_mueller()
	raggio_iniz.mueller_product(my_interazione.muellers.loc[0, 'M_rif'])
	raggio_finale_M.append(raggio_iniz)
	
    #con biquaternioni
	biq_iniz = Quaternion(raggio_iniz.parameters[0], raggio_iniz.parameters[1]*I, raggio_iniz.parameters[2]*I, raggio_iniz.parameters[3]*I)
	h_riflessione = my_interazione.biquaternions.loc[0, 'h_rif']
	biq_fin = h_riflessione.mul(biq_iniz) #giusto il verso della moltipl?
	raggio_finale = [biq_fin.a, biq_fin.b*(-1j), biq_fin.c*(-1j), biq_fin.d*(-1j)]
	raggio_finale_h.append(raggio_finale)
    
    #per ora inutile
	prodotto_scalare = scalar_prod(biq_iniz, biq_fin)
	total_phase = cmath.phase(prodotto_scalare)
	result = my_interazione.muellers.loc[0, 'M_rif'][0,0]*prodotto_scalare
    
raggio_finale_h_I = []   
raggio_finale_h_Q = [] 

raggio_finale_M_I = [] 
raggio_finale_M_Q = [] 
    
for i in range(nvalues):
    
    raggio_finale_h_I.append(raggio_finale_h[i][0])    
    raggio_finale_h_Q.append(raggio_finale_h[i][1])
    
    raggio_finale_M_I.append(raggio_finale_M[i].parameters[0]) 
    raggio_finale_M_Q.append(raggio_finale_M[i].parameters[1])
    
theta = np.linspace(0, 90, nvalues)

plt.title("con Mueller")
plt.xlabel("angolo incidente [°]")
plt.ylabel("I | Q")

plt.plot(theta, raggio_finale_M_I, label='I') 
plt.plot(theta, raggio_finale_M_Q, label='Q')	
plt.legend()
plt.grid(True)
plt.show()

plt.title("con biquaternioni")
plt.xlabel("angolo incidente [°]")
plt.ylabel("I | Q")
plt.plot(theta, raggio_finale_h_I, label='I') 
plt.plot(theta, raggio_finale_h_Q, label='Q')	
plt.legend()
plt.grid(True)
plt.show()