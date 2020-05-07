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
from biquaternion import Biquaternion

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
        self.precisione = precisione
        
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
            
            print(h_mat) #NOTA: h_mat è autovettore di una matrice multipla dell'identità
            print()
            
            biq_mat = Biquaternion(h_mat[0], h_mat[1], h_mat[2], h_mat[3])
            biq_rif = Biquaternion(h_rif[0], h_rif[1], h_rif[2], h_rif[3])
            biq_tra_dw = Biquaternion(h_tra_dw[0], h_tra_dw[1], h_tra_dw[2], h_tra_dw[3])
            biq_tra_up = Biquaternion(h_tra_up[0], h_tra_up[1], h_tra_up[2], h_tra_up[3])
'''prova'''

my_campione = campione(1.1, [2.1+0.00000001j,0.02, 4.1+0.00000002j,0.04, 5.1+0.00000003j,0.05, 9.1])
my_sorgente = sorgente(0.3, 0.75, 0.43)

my_interazione = interazione(0.0001, my_campione,my_sorgente, 0.15)
my_interazione.materials_to_mueller()