import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from campione import campione
from interazione import *
from stokes import stokes_vector
from sorgente import sorgente
from sympy import I
from sympy.algebras.quaternion import Quaternion


#indici rifrzione
n0 = 1
n1 = 1.5+0.00000002j
sp1 = 0.00002
n2 = 1+0.00000005j
pi = math.pi

#definisco il campione
camp = campione(n0, [n1, sp1, n2])

#definisco la sorgente
sorg = sorgente(1.95954488, 1, 0.78539163)

#definisco il raggio iniziale
riniz = stokes_vector()
riniz.generic_polarization(1, 1)

#definisco il quaternione del raggio inziale  
s = Quaternion( riniz.I(), riniz.Q()*I, riniz.U()*I, riniz.V()*I )	#quaternione corrispondente al vett Stoke

#Liste per risultati quaternioni
Psi = []
Delta = []
S = []
Q = []
U = []
V = []

nvalues = 100
theta = np.linspace(0, pi/2, nvalues)

for i in range(nvalues):
    
	print()
	print('***********************************************')
	print("Angolo di incidenza ", round(180/pi * theta[i], 2), "°")
	print()
	
    #calcolo le interazioni di tutti i possibili raggi con le interfacce
	inter = interazione(0.0001, camp, sorg, riniz)  #inizializzo oggetto interazione
	inter.materials_to_jones(theta[i])
    
	n=0
	while inter.nraggi != 0:
		inter.propagazione()
		print("iterazione = ", n)
		n += 1

	#determino l'interferenza tra tutti i raggi ottenuti
	hfin = inter.interference()
 
	#grandezze ellissometriche e vettore di Stokes finale
	psi, delta, rfin = grandell([hfin], s, svfinal=True)

	#Salva i risultati con quaternioni
	S.append( rfin.I().real )
	Q.append( rfin.Q().real )
	U.append( rfin.U().real )
	V.append( rfin.V().real )
	#Parametri ellissometrici
	Psi.append( psi )
	Delta.append( delta )


#Grafico dei risultati con quaternioni
plt.title("Calcolo tramite quaternioni")
plt.xlabel("Angolo incidente [°]")
plt.ylabel("Componenti vettori ellissometriche")
plt.plot(theta, S , label='I')
plt.plot(theta, Q , label='Q')
plt.plot(theta, U , label='U')
plt.plot(theta, V , label='V')
plt.grid(True)
plt.legend()
plt.show()

Psi = np.dot(Psi, 180/pi)
Delta = np.dot(Delta, 180/pi)
theta = np.dot(theta, 180/pi)

#con i quaternioni
plt.title("Calcolo tramite quaternioni")
plt.xlabel("Angolo incidente [°]")
plt.ylabel("Grandezze ellissometriche")
plt.plot(theta, Psi, label='$\Psi$ quaternioni', linestyle='dashed')
plt.plot(theta, Delta, label='$\Delta$ quaternioni')
plt.grid(True)
plt.legend()
plt.show()



