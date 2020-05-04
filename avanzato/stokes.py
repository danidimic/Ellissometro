import math
import cmath
import numpy as np

pi=math.pi

class stokes_vector:

	''' 
	A partire dal campo elettrico lungo x e lungo y
	determino i parametri di Stokes dell'onda polarizzata-
	parzialmente o totalmente- o non polarizzata 
	I campi elettrici sono descrivibili mediante le equazioni seguenti:

	Ex = E0x cos( omega*t+deltax)
	Ey = E0y cos( omega*t+deltay)
	
	E0x, E0y sono le ampiezze istantanee
	deltax, deltay sono le fasi istantanee 
	omega è la frequenza angolare 
	'''

	def __init__(self):
		self.parameters = []

	#generica polarizzazione
	def generic_polarization(self, E0x, E0y, deltax=0, deltay=0):
		delta = deltay - deltax
		I = E0x**2 + E0y**2 			#intensità totale
		Q = E0x**2 - E0y**2 			#intensità polarizzata orizzontalmente o verticalmente
		U = 2*E0x*E0y*np.cos(delta)		#intensità polarizzata a 45 gradi 
		V = 2*E0x*E0y*np.sin(delta)		#intensità polarizzata circolarmente

		self.parameters = [I, Q, U, V]

	#polarizzazione lineare con angolo theta rispetto all'asse x
	def linear_polarization(self, E, theta):
		E0x = E*np.cos(theta)
		E0y = E*np.sin(theta)
		self.generic_polarization(E0x, E0y)

	#polarizzazione circolare
	def circular_polarization(self, E, clockwise=False):
		if clockwise == False:
			delta = pi/2
		else:
			delta = -pi/2
	
		E0x = E/np.sqrt(2)
		E0y = E/np.sqrt(2)
		self.generic_polarization(E0x, E0y, deltay=delta)
	
	#luce non polarizzata
	def unpolarized(self, I):
		self.parameters = [I, 0, 0, 0]

	
	#grado di polarizzazione
	def polarization_degree(self):
		I = self.parameters[0]
		Q = self.parameters[1]	
		U = self.parameters[2]		
		V = self.parameters[3]			
		
		return np.sqrt( Q*Q+U*U+V*V ) / I
	
	#Calcolo di Psi, Chi dell'ellisse di polarizzazione
	def polarization_elipse(self):
		I = self.parameters[0]
		Q = self.parameters[1]	
		U = self.parameters[2]		
		V = self.parameters[3]	

		psi = 0.5*np.arctan(U/Q)
		chi = 0.5*np.arctan(V/np.sqrt(Q*Q+U*U))
		return [psi, chi]


v = stokes_vector()
v.unpolarized(1)
print( v.polarization_degree() )
