import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from sorgente import sorgente
from biquaternion import *
from campione import *
from interazione import *
from stokes import stokes_vector
from loadBar.bar import IncrementalBar
pi = math.pi

#definisco il campione
camp = loadCampione()

#definisco la sorgente
sorg = sorgente(lenght=630.)

#definisco il raggio iniziale
riniz = stokes_vector()
riniz.generic_polarization(1, 1)

#definisco il quaternione del raggio inziale 
s = StokesToBiq(riniz)

#Liste per risultati quaternioni
Psi = []
Delta = []
S = []
Q = []
U = []
V = []

nvalues = 200
theta = np.linspace(0, pi/2, nvalues)#[0, 1.46]#

print("Calcolo delle grandezze ellissometriche al variare dell'angolo di incidenza")
print()

print("Caratteristiche della sorgente:")
print("Energia =", round(sorg.energia, 3), "eV")
print("Lunghezza d'onda =", round(sorg.lunghezza, 2), "nm")
print()

print("Caratteristiche del campione:")
print("Numero di strati =", camp.strati)
print("--- Strato iniziale ---")
print("Indice di rifrazione iniziale =", camp.nc[0])

for i in range(camp.strati):
	print("--- Strato", i+1, "---")
	print("Indice di rifrazione =", camp.nc[i+1])
	print("Spessore =", round(camp.spessori[i+1]*10**9, 2), "nm")

print("--- Substrato ---")
print("Indice di rifrazione del substrato =", camp.nc[camp.strati+1])
print()

bar = IncrementalBar('Progresso:', suffix='%(percent)d%%', max=nvalues)
for i in range(nvalues):
    
	bar.next()
	consoleOut = "   Angolo di incidenza: " + str(round(180/pi * theta[i], 2)) + "°"
	print(consoleOut, end="\r", flush=True)

    #calcolo le interazioni di tutti i possibili raggi con le interfacce
	inter = interazione(0.01, camp, sorg, riniz)  #inizializzo oggetto interazione
	inter.materials_to_jones(theta[i])
	#print('############# theta = ', theta[i], '#############')
	if camp.strati != 0:
		while inter.nraggi != 0:
			inter.propagazione()

		#determino l'interferenza tra tutti i raggi ottenuti
		hfin = inter.interference()
	else:
		hfin = inter.biquaternions.loc[0, 'h_rif_dw']

	#grandezze ellissometriche e vettore di Stokes finale
	psi, delta, rfin = grandell([hfin], s, svfinal=True)


	#Salva i risultati con quaternioni
	S.append( rfin.I().real )
	Q.append( rfin.Q().real )
	U.append( rfin.U().real )
	V.append( rfin.V().real )

	#Parametri ellissometrici
	Psi.append( psi.real )
	Delta.append( delta.real )

bar.finish()



#conversione in gradi
Psi = np.dot(Psi, 180/pi)
Delta = np.dot(Delta, 180/pi)
theta = np.dot(theta, 180/pi)

#write data in 2 file
file = open("Delta.txt","w")

for i in range(nvalues):
    line = str(theta[i]) + "\t" + str(Delta[i]) + "\n"
    file.write(line)
    
file.close()

file = open("Psi.txt","w")

for i in range(nvalues):
    line = str(theta[i]) + "\t" + str(Psi[i]) + "\n"
    file.write(line)
    
file.close()



#Grafico delle componenti del vettore di Stokes
plt.title("Componenti vettori di Stokes")
plt.xlabel("Angolo incidente [°]")
plt.ylabel("[unità arbitrarie]")
plt.plot(theta, S , label='I')
plt.plot(theta, Q , label='Q')
plt.plot(theta, U , label='U')
plt.plot(theta, V , label='V')
plt.grid(True)
plt.legend()
plt.show()

#Grafico delle grandezze ellissometriche
plt.title("Grandezze ellissometriche")
plt.xlabel("Angolo incidente [°]")
plt.ylabel("[°]")
plt.plot(theta, Psi, label='$\Psi$')
plt.plot(theta, Delta, label='$\Delta$')
plt.grid(True)
plt.legend()
plt.show()

