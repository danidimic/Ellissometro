import math
import numpy as np

nval = 100
pi = math.pi
theta = np.linspace(0, 180, nval)
names = []

for i in range(nval):
	#t = str(round(theta[i], 3)) + "\n"
	t = str(round(theta[i], 3))
	names.append(t)	

#scrittura su file LOOP.DAT
#out = open("LOOP.DAT", 'w')
#out.writelines(names)
#out.close()

resE1 = []
resE2 = []

for i in range(nval):
	name1 = "results/E1_theta" + names[i] + ".txt"
	name2 = "results/E2_theta" + names[i] + ".txt"
	l1, E1 = np.loadtxt(name1, unpack=True)
	l2, E2 = np.loadtxt(name2, unpack=True)

	mean = (E1[0] + E1[1]) / 2.
	resE1.append(mean)
	mean = (E2[0] + E2[1]) / 2.
	resE2.append(mean)


np.savetxt("E1.txt", resE1, delimiter='  ')
np.savetxt("E2.txt", resE2, delimiter='  ')
