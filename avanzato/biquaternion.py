import math
import cmath
import numpy as np

class Biquaternion:

	def __init__(self, R, I, J, K):
		self.R = complex(R)
		self.I = complex(I)
		self.J = complex(J)
		self.K = complex(K)


	def __str__(self):
		rep = str(self.R) + "  " + str(self.I)+"I" + "  " + str(self.J)+"J" + "  " + str(self.K)+"K" 
		return rep

	def scalar_prod(self, other):
		return self.R*other.R + self.I*other.I + self.J*other.J + self.K*other.K
	
	def prod(self, other):
		r1 = self.R*other.R + self.I*other.I + self.J*other.J + self.K*other.K
		i1 = 1j * ( self.R*other.I + self.I*other.R + 1j*self.J*other.K - 1j*self.K*other.J )
		j1 = 1j * ( self.R*other.J + self.J*other.R - 1j*self.I*other.K + 1j*self.K*other.I )
		k1 = 1j * ( self.R*other.K + self.K*other.R + 1j*self.I*other.J - 1j*self.J*other.I )
		return Biquaternion(r1,i1,j1,k1)
            
	def exp(self):
        # exp(q) = e^a(cos||v|| + v/||v||*sin||v||)
		q = self
		vector_norm = np.sqrt(q.I**2 + q.J**2 + q.K**2)
		r = np.exp(q.R) * np.cos(vector_norm)
		i = np.exp(q.R) * np.sin(vector_norm) * q.I / vector_norm
		j = np.exp(q.R) * np.sin(vector_norm) * q.J / vector_norm
		q = np.exp(q.R) * np.sin(vector_norm) * q.K / vector_norm
        
		return Biquaternion(r, i, j, q)

    
q = Biquaternion(1, 3+1j, 2j, 0)
k = Biquaternion(7, 1j, 1j, 2)
scalar = q.scalar_prod(k)
exp_= k.exp()
print(exp_)

		

