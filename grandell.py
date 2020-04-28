import numpy as np
import array as arr
import math
import cmath
from campione import campione
from sorgente import sorgente
from interazione import interazione


class ellissometro:
		
	def __init__(self, precisione, theta0):
		self.precisione = precisione
		self.theta0 = theta0

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

		sorg = sorgente(cont[0], cont[1], cont[2])

		inter = interazione(0.1, camp, sorg, self.theta0)
		inter.inizializza()

		while (len(inter.ii)!=0):
			inter.propagazione()

		r_pi = inter.somma_pi/(math.tan(sorg.psi_0)*cmath.exp(1j*sorg.delta_0))
		r_sigma = inter.somma_sigma

		psi1 = math.atan( abs(r_pi/r_sigma) );
		delta1 = -cmath.phase( r_pi/r_sigma );

		#psi1=math.atan(abs(r_sigma./r_pi));
		#delta1=angle(r_sigma./r_pi);
		#delta1=2*pi*(delta1<0)+delta1;
        
		if delta1<0:
			delta1= delta1 + 2*math.pi

		#print()
		print('r_pi = ', r_pi)
		print('r_sigma = ', r_sigma)
		print('delta1 = ', delta1)
		print('psi1 = ', psi1)

		results = [r_pi, r_sigma, delta1, psi1]
		return results

 
E = ellissometro(0.1, 0.2)
E.grandell()
