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

	def __init__(self, I=0, Q=0, U=0, V=0):
		I = np.complex128(I)
		Q = np.complex128(Q)
		U = np.complex128(U)
		V = np.complex128(V)
		self.parameters = [I, Q, U, V]

	def __str__(self):
		rep = "Vettore di Stokes di componenti: \n"
		rep += "I = " + str(self.parameters[0]) + "\n"
		rep += "Q = " + str(self.parameters[1]) + "\n"
		rep += "U = " + str(self.parameters[2]) + "\n"
		rep += "V = " + str(self.parameters[3])
		return rep

	def I(self):
		return self.parameters[0]

	def Q(self):
		return self.parameters[1]

	def U(self):
		return self.parameters[2]

	def V(self):
		return self.parameters[3]



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
		S1 = self.parameters[1]	
		S2 = self.parameters[2]		
		S3 = self.parameters[3]	

		psi = 0.5*np.arctan(S2/S1)
		chi = 0.5*np.arctan(S3/np.sqrt(S1*S1+S2*S2))
		return [psi, chi]

	#Componenti del campo elettrico
	def electric_components(self):
		S0 = self.I()
		S1 = self.Q()

		Ex = np.sqrt( (S0+S1)/2 )
		Ey = np.sqrt( (S0-S1)/2 )
		return [Ex, Ey]



	#Calcolo del parametro ellissometrico Psi
	def ellipsometric_Psi(self):
		Ex, Ey = self.electric_components()
		return np.arctan(Ey/Ex)

	#Calcolo del parametro ellissometrico Delta
	def ellipsometric_Delta(self):
		Ex, Ey = self.electric_components()
		return (self.U()-1j*self.V())/ (2*Ey*Ex)
      
	def alfa(self):
		return 0.5*abs(np.arccos(self.Q()/self.I()))
    
	def delta(self):
		delta = np.arctan(self.V()/self.U())
		if delta < 0:
			delta += 2*math.pi
		return delta



	#Prodotto per una matrice mueller
	def mueller_product(self, mat):
		product = np.dot(mat, self.parameters)
		self.parameters = product



