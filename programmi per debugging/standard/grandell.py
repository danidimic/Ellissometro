import numpy as np
import array as arr
import math
import cmath
from campione import campione
from sorgente import sorgente
from interazione import interazione

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)

class ellissometro:
		
	def __init__(self, precisione):
		self.precisione = precisione
		self.sorg_val = 0

	def grandell(self, theta0, write=False):
		#lettura dei dati del campione
		f = open("campione.txt", 'r')
		cont = f.readlines()
		for i in range(len(cont)):
			cont[i] = complex(cont[i])
		#inizializzo il campione
		camp = campione(cont[0], cont[1:] )

		#lettura dei dati della sorgente
		f = open("sorgente.txt", 'r')
		cont = f.readlines()
		for i in range(len(cont)):
			cont[i] = float(cont[i])
		#inizializzo la sorgente
		sorg = sorgente(cont[0], cont[1], cont[2])
		self.sorg_val = sorg

		#stampo le caratteristiche di campione e sorgente
		if write==True:
			print("Analisi delle grandezze elissometriche del seguente campione")
			print("numero di strati: ", camp.strati)
			print("indice di rifrazione iniziale: ", camp.nc[0])
			for i in range(camp.strati):
				print("indice di rifrazione dello strato", i+1, ": ", camp.nc[i+1])
				print("spessore dello strato", i+1, ": ", camp.spessori[i+1])
			print()
			print("Caratteristiche della sorgente usata per l'analisi")
			print("lunghezza d'onda: ", h*c/sorg.energia, " m")
			print("Psi_0: ", sorg.psi_0)
			print("delta_0: ", sorg.delta_0)
	
			
		inter = interazione(0.0001, camp, sorg, theta0)
		inter.inizializza()

		while (len(inter.ii)!=0):
			inter.propagazione()

		r_pi = inter.somma_pi/(math.tan(sorg.psi_0)*cmath.exp(1j*sorg.delta_0))
		r_sigma = inter.somma_sigma

		psi1 = math.atan( abs(r_pi/r_sigma) );
		delta1 = -cmath.phase( r_pi/r_sigma );

		if delta1<0:
			delta1= delta1 + 2*math.pi

		results = [r_pi, r_sigma, delta1, psi1]
		return results
