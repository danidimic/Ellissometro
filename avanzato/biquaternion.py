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



q = Biquaternion(1, 3+1j, 2j, 0)
print(q)

		

