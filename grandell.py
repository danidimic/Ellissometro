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

		sorg = sorgente(cont[0], cont[1], cont[2])

		inter = interazione(0.1, camp, sorg)
		inter.inizializza(math.pi/4)

		while (len(inter.ii)==0):
			inter.propagazione()

		r_pi = inter.somma_pi/(math.tan(sorg.psi_0)*cmath.exp(1j*sorg.delta_0))
		r_sigma = inter.somma_sigma
		print(inter.somma_sigma)

		psi1 = math.atan( abs(r_pi/r_sigma) );
		delta1 = -angle( r_pi/r_sigma );

		# psi1=atan(abs(r_sigma./r_pi));
		# delta1=angle(r_sigma./r_pi);
		#delta1=2*pi*(delta1<0)+delta1;

E = ellissometro(1)
E.grandell()
