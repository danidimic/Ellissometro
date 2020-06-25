import math

pi = math.pi
c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)


class sorgente:

	def __init__(self, lenght, energy=0, frequency=0):
		
		#lunghezza d'onda della sorgente in nm
		self.lunghezza = lenght
		l = lenght*1e-9

		#energia della sorgente in eV
		self.energia = h*c / l

		#frequenza della sorgente in Hz
		self.frequenza = c / l

		#frequenza angolare in Hz normalizzata a c
		self.wsuc = 2*pi*self.frequenza/c

