import numpy as np
import array as arr
import math
import cmath
from campione import campione
from sorgente import sorgente
from interazione import interazione

theta0 = 45

class ellissometro:
		
	def __init__(self, precisione):
		self.precisione = precisione

	def grandell(self):

		f = open("campione.txt", 'r')
		cont = f.readlines()
		for i in range(len(cont)):
			cont[i] = complex(cont[i])

		camp = campione(cont[0], cont[1:] ) 

		f = open("sorgente.txt", 'r')
		cont = f.readlines()
		for i in range(len(cont)):
			cont[i] = float(cont[i])
			#print(cont[i])

		sorg = sorgente(cont[0], cont[1], cont[2])

		inter = interazione(0.0001, camp, sorg)
		inter.inizializza(0.2)
		#print('funzia')

		while (len(inter.ii)!=0):
			#print('funzia')
			inter.propagazione()

		r_pi = inter.somma_pi/(math.tan(sorg.psi_0)*cmath.exp(1j*sorg.delta_0))
		r_sigma = inter.somma_sigma
		#print(inter.somma_sigma)

		psi1 = math.atan( abs(r_pi/r_sigma) );
		delta1 = -cmath.phase( r_pi/r_sigma );
        
		if delta1<0:
			delta1= delta1 + 2*math.pi

		#psi1=math.atan(abs(r_sigma./r_pi));
		#delta1=angle(r_sigma./r_pi);v#delta1=2*pi*(delta1<0)+delta1;
		print('r_pi = ', r_pi)
		print('r_sigma = ', r_sigma)
		print('delta1 = ', delta1)
		print('psi1 = ', psi1)
        
E = ellissometro(0.0001)
E.grandell()


