import math
import cmath
import numpy as np
import scipy.linalg as la
#from pyquaternion import Quaternion

pi=math.pi

#Definisco le matrici di Pauli sigma_i
sigma0 = np.matrix([ [1+0j, 0+0j], [0+0j, 1+0j] ])
sigma1 = np.matrix([ [1+0j, 0+0j], [0+0j,-1+0j] ])
sigma2 = np.matrix([ [0+0j, 1+0j], [1+0j, 0+0j] ])
sigma3 = np.matrix([ [0+0j, 0-1j], [0+1j, 0+0j] ])
#Costruisco un vettore di matrici di Pauli
PauliMat = [ sigma0, sigma1, sigma2, sigma3 ]

#Definisco la matrice A
A = np.matrix([ [1,0,0,1],[1,0,0,-1],[0,1,1,0],[0,1j,-1j,0] ])
#Definisco la matrice inversa di A
Ainv = np.linalg.inv( np.matrix(A) )

#Calcolo la matrice covarianza data la matrice di Mueller
def covariance_matrix(mat):
	H = np.zeros( (4,4), dtype=np.complex64 )

	for i in range(4):
		for j in range(4):
			Pij = A * np.kron( PauliMat[i], np.conj(PauliMat[j]) ) * Ainv
			H += 1./4. * mat[i][j] * Pij
	return H


#Calcolo il vettore covariante data la matrice covarianza
def covariance_vector(H):

	h = np.zeros(4, dtype=np.complex64)
	rank = np.linalg.matrix_rank(H)

    #Se la matrice è di rango=1 c'è un solo autovalore non nullo
	if rank==1 :
		eigvals, eigvecs = la.eig(H)

		for i in range(4):
			if abs(eigvals[i]) > 1e-5 :
				indexeig = i

        #l'autovettore è scritto come vettore colonna
		h[0] = eigvecs[0][indexeig]
		h[1] = eigvecs[1][indexeig]
		h[2] = eigvecs[2][indexeig]
		h[3] = eigvecs[3][indexeig]
            
		norm = la.norm(h) #inutile? (la.eig li restituisce già normalizzati)
		h /= np.sqrt(norm)
            
	return h

#Calcolo la matrice Z dal vettore covariante
def Zeta_matrix(h):
	t = h[0]
	a = h[1]
	b = h[2]
	c = h[3]

	Z = np.matrix([ [t, a, b, c],
					[a, t, -c*1j, b*1j],
					[b, c*1j, t, -a*1j],
					[c, -b*1j, a*1j, t] ])

	return Z

#PER ORA TOLTO, USIAMO BIQUATERNIONI
#Definisco h tramite quaternioni
#def h_quaternion(h):
#	t = h[0]
#	a = h[1]
#	b = h[2]
#	c = h[3]

#	qh = Quaternion(t, 1j*a, 1j*b, 1j*c)
#	return qh






"""
#########################################
#prova di calcolo con risultato a pag 85
M = [[30., 8., 4., 12.],
	[8., 26., 8., 0.],
	[12., 0., 12., 26.],
	[4., 8., -22., 12.]]

M = np.dot((1./30.),M)

H = covariance_matrix(M)
h = covariance_vector(H)
Z = Zeta_matrix(h)
#qh = h_quaternion(h)
#qh_conj = qh.conjugate

print("vettore h:")
print(h)
print()
print("matrice Z: ")
print(np.asmatrix(Z))
print()
print("verifica: M = ZZ*")
print(Z*np.conj(Z))
print()
print(np.conj(Z)*Z)
print()
print("matrice M")
print(M)

print()
#print(qh)
#print(qh_conj)
"""
