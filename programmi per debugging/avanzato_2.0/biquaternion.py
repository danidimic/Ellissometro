import math
import numpy as np
from stokes import stokes_vector

class Biquaternion:

	#costruttore classe Biquaternion
	def __init__(self, R=0, I=0, J=0, K=0):
		self.a = np.complex128(R)
		self.b = np.complex128(I)
		self.c = np.complex128(J)
		self.d = np.complex128(K)

	#Overload funzione print su Biquaternion
	def __str__(self):
		rep = str(self.a) + "  " + str(self.b)+"I" + "  " + str(self.c)+"J" + "  " + str(self.d)+"K" 
		return rep

	#Overload operatore somma + per Biquaternion 
	def __add__(self, o):
		A = self.a + o.a
		B = self.b + o.b
		C = self.c + o.c
		D = self.d + o.d
		return Biquaternion(A, B, C, D) 

	#Overload operatore somma += per Biquaternion
	def __iadd__(self, o):
		self.a += o.a
		self.b += o.b
		self.c += o.c	
		self.d += o.d
		return self

	#Overload operatore differenza - per Biquaternion 
	def __sub__(self, o):
		A = self.a - o.a
		B = self.b - o.b
		C = self.c - o.c
		D = self.d - o.d
		return Biquaternion(A, B, C, D)

	#Overload operatore moltiplicazione * per Biquaternion 
	def __mul__(self, o):
		prod = self.mul(o)
		return prod


	#Metodo per moltiplicazione di quaternioni
	def mul(self, q):
		a1 = self.a
		b1 = self.b
		c1 = self.c
		d1 = self.d

		a2 = q.a
		b2 = q.b
		c2 = q.c
		d2 = q.d
	
		R = a1*a2 - b1*b2 - c1*c2 - d1*d2
		I = a1*b2 + b1*a2 + c1*d2 - d1*c2
		J = a1*c2 + c1*a2 + d1*b2 - b1*d2
		K = a1*d2 + d1*a2 + b1*c2 - c1*b2 
		return Biquaternion(R, I, J, K)

	#Metodo per prodotto scalare tra quaternioni
	def scalar_prod(self, other):
		return self.R*other.R + self.I*other.I + self.J*other.J + self.K*other.K



#Biquaternion to Stokes vector
def BiqToStokes(quaternion):
	a = quaternion.a
	b = quaternion.b*-1j
	c = quaternion.c*-1j
	d = quaternion.d*-1j
	return stokes_vector(a, b, c, d)

#Stokes vector to Biquaternion
def StokesToBiq(stokesVector):
	a = stokesVector.I()
	b = stokesVector.Q()*1j
	c = stokesVector.U()*1j
	d = stokesVector.V()*1j
	return Biquaternion(a, b, c, d)


#Prodotto scalare fra biquaterioni
def scalar_prod(q1, q2):
	return q1.a*q2.a + q1.b*q2.b + q1.c*q2.c + q1.d*q2.d


#Coniugato di un quaternione
def conjugate(h):
	tau   =  np.conj(h.a)
	alfa  =  np.conj(h.b)*1j
	beta  =  np.conj(h.c)*1j
	gamma =  np.conj(h.d)*1j
	return Biquaternion( tau, alfa*1j, beta*1j, gamma*1j )

