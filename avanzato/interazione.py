import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from mueller import MM
from stokes import stokes_vector
import pandas as pd
from sorgente import sorgente
from campione import campione
from matrix import *
from sympy.algebras.quaternion import Quaternion

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)

'''
n = 100
pi=math.pi
theta = np.linspace(0, pi/2, n)

v = stokes_vector()
v.linear_polarization(1, pi/4)

Psi = []
Delta = []

for i in range(n):
	m = MM(1+0j, 1.5+0.00000002j, theta[i], 0.5, 0.00002)
	mat_ref = m.mueller_reflection()
	mat_tra = m.mueller_transmission()
	mat_lay = m.mueller_layer()
	
	v.layer_interaction(mat_ref, mat_tra, mat_lay)

	psi, delta = v.ellipsometric_parameters()

	Psi.append(psi)
	Delta.append(delta)	


plt.plot(theta, Psi, label="$\Psi$")
plt.plot(theta, Delta, label="$\Delta$")
plt.legend()
plt.grid(True)
plt.show()
'''

class interazione:
    
    def __init__(self, precisione, campione, sorgente, theta0):
        
        self.campione = campione
        self.sorgente = sorgente
        self.precisione = precisione
        
        self.theta0 = theta0
        
        #inizializzo un Data Frame dove metterò le matrici di Mueller
        righe = np.arange(0, self.campione.strati + 1) #il vuoto è lo strato 0
        col = ['M_mat', 'M_rif', 'M_tra_dw', 'M_tra_up']
        self.muellers = pd.DataFrame(0, righe, col)
        self.muellers = self.muellers.astype(object)
        
        #inizializzo un Data Frame dove metterò i biquaternioni
        col = ['h_mat', 'h_rif', 'h_tra_dw', 'h_tra_up']
        self.biquaternions = pd.DataFrame(0, righe, col)
        self.biquaternions = self.biquaternions.astype(object)
        
        
    def materials_to_mueller(self):
        
        '''
        Raccoglie nel Data Frame, con indici di riga riferiti agli strati (0 
        è il vuoto): matrice di Mueller relativa alla propagazione nel materiale 
        (M_mat), matrice di Mueller relativa alla riflessione sull'interfaccia (in
        basso) (M_rif), matrice di Mueller relativa alla trasmissione attraverso 
        l'interfaccia (in basso) (M_tra). NOTA BENE: in questa versione del programma
        considero il mezzo isotropo, per cui gli effetti della propagazione nel mezzo
        non dipendono dal segno dell'angolo, né dal fatto che la luce propaghi dal basso
        verso l'alto o viceversa.
        '''
        omega=2*math.pi*(self.sorgente.energia)/h;
        wsuc=omega/c;
        
        #Il primo theta è quello con cui la luce parte dalla sorgente
        theta_start = self.theta0
        
        for i in range(self.campione.strati+1):
            
            #inizializzo al layer a cui sono arrivato
            n0 = self.campione.nc[i] # questo deve essere coerente con...
            n1 = self.campione.nc[i+1] #...
            spessore = self.campione.spessori[i] # ...questo!
            mueller_class = MM(n0, n1, theta_start, wsuc, spessore) #definisco la classe ad ogni passaggio, meglio cambiare?
            
            #assegno le matrici di Mueller al Data Frame:
            self.muellers.loc[i, 'M_mat'] = np.array(mueller_class.mueller_layer())
            self.muellers.loc[i, 'M_rif'] = np.array(mueller_class.mueller_reflection())
            self.muellers.loc[i, 'M_tra_dw'] = np.array(mueller_class.mueller_transmission())
            
            #la matrice di Mueller per la trasmissione per ragggi diretti dal basso verso l'alto è
            #diversa da quella per raggi diretti dall'alto verso il basso
            #quindi:
            #aggiorno theta
            theta_start = mueller_class.theta1
            #e me la ricalcolo
            mueller_class = MM(n1, n0, theta_start, wsuc, spessore)
            self.muellers.loc[i, 'M_tra_up'] = np.array(mueller_class.mueller_transmission())
            
            #print('rif alto a basso: ', mueller_class.mueller_reflection())
            #print('tras alto a basso: ', mueller_class.mueller_transmission())
            #print('rif basso a alto: ', mueller_class.mueller_reflection())
            #print('tras basso a alto: ', mueller_class.mueller_transmission())
            #print()
           
        for i in range(self.campione.strati+1):
            #calcolo matrici di covarianza
            H_mat = covariance_matrix(self.muellers.loc[i, 'M_mat'])
            H_rif = covariance_matrix(self.muellers.loc[i, 'M_rif'])
            H_tra_dw = covariance_matrix(self.muellers.loc[i, 'M_tra_dw'])
            H_tra_up = covariance_matrix(self.muellers.loc[i, 'M_tra_up'])
            
            #calcolo i vettori di covarianza
            h_mat = covariance_vector(H_mat)
            h_rif = covariance_vector(H_rif)
            h_tra_dw = covariance_vector(H_tra_dw)
            h_tra_up = covariance_vector(H_tra_up)
            
            #print(H_rif) #NOTA: h_mat è autovettore di una matrice multipla dell'identità
            print()
            
            self.biquaternions.loc[i, 'h_mat'] = Quaternion(h_mat[0], h_mat[1]*1j, h_mat[2]*1j, h_mat[3]*1j)
            self.biquaternions.loc[i, 'h_rif'] = Quaternion(h_rif[0], h_rif[1]*1j, h_rif[2]*1j, h_rif[3]*1j)
            self.biquaternions.loc[i, 'h_tra_dw'] = Quaternion(h_tra_dw[0], h_tra_dw[1]*1j, h_tra_dw[2]*1j, h_tra_dw[3]*1j)
            self.biquaternions.loc[i, 'h_tra_up'] = Quaternion(h_tra_up[0], h_tra_up[1]*1j, h_tra_up[2]*1j, h_tra_up[3]*1j)
'''prova'''


my_sorgente = sorgente(1.95954488, 1, 0.78539163)

#indici rifrzione
n0 = 1
n1 = 2.1+1j

#Brewster angle
br = np.arctan(n1/n0)

raggio_iniz = stokes_vector()
raggio_iniz.generic_polarization(1,1.5, 0.12, 0.34)

#print('raggio iniziale: ',raggio_iniz)

my_campione = campione(n0, [n1])

my_interazione = interazione(0.0001, my_campione,my_sorgente, br)
my_interazione.materials_to_mueller()

raggio_iniz.mueller_product(my_interazione.muellers.loc[0, 'M_rif'])
#print()
#print(my_interazione.muellers.loc[0, 'M_rif'])
#print()

raggio_finale = raggio_iniz
#print('raggio finale: ',raggio_finale)
#print()
#formalismo biquaternioni
biq_iniz = Quaternion(raggio_iniz.parameters[0], raggio_iniz.parameters[1]*1j, raggio_iniz.parameters[2]*1j, raggio_iniz.parameters[3]*1j)
h_riflessione = my_interazione.biquaternions.loc[0, 'h_rif']

print('biquaternione riflessione: ',h_riflessione)
biq_fin = h_riflessione.mul(biq_iniz)

def scalar_prod(q1, q2):
	return q1.a*q2.a + q1.b*q2.b + q1.c*q2.c + q1.d*q2.d

prodotto_scalare = scalar_prod(biq_iniz, biq_fin)
total_phase = cmath.phase(prodotto_scalare)

print('hs: ', biq_fin)
print('prodotto scalare: ',prodotto_scalare)

print(total_phase)

