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
		return self.R*other.R - self.I*other.I - self.J*other.J - self.K*other.K
            

q = Biquaternion(1, 3+1j, 2j, 0)
k = Biquaternion(0, 1, 1j, 2)
scalar = q.scalar_prod(k)
print(scalar)

		

