import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
import pandas as pd
from mueller import MM
from stokes import stokes_vector
from sorgente import sorgente
from campione import campione
from matrix import *
from sympy import I
from sympy.algebras.quaternion import Quaternion


c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)

class interazione:
    
    def __init__(self, precisione, campione, sorgente):
        
        self.campione = campione
        self.sorgente = sorgente
        self.precisione = precisione
        
        #inizializzo un Data Frame dove metterò le matrici di Mueller
        righe = np.arange(0, self.campione.strati + 1) #il vuoto è lo strato 0
        col = ['M_mat', 'M_rif_up', 'M_rif_dw', 'M_tra_dw', 'M_tra_up']
        self.muellers = pd.DataFrame(0, righe, col)
        self.muellers = self.muellers.astype(object)
        
        #inizializzo un Data Frame dove metterò i biquaternioni
        col = ['h_mat', 'h_rif_up', 'h_rif_dw', 'h_tra_dw', 'h_tra_up']
        self.biquaternions = pd.DataFrame(0, righe, col)
        self.biquaternions = self.biquaternions.astype(object)
        
        
    def materials_to_mueller(self, theta0):
        
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
        theta_start = theta0
        
        for i in range(self.campione.strati+1):
            
            #inizializzo al layer a cui sono arrivato
            n0 = self.campione.nc[i] # questo deve essere coerente con...
            n1 = self.campione.nc[i+1] #...
            spessore = self.campione.spessori[i] # ...questo!
            mueller_class = MM(n0, n1, theta_start, wsuc, spessore) #definisco la classe ad ogni passaggio, meglio cambiare?
            
            #assegno le matrici di Mueller al Data Frame:
            self.muellers.loc[i, 'M_mat'] = np.array(mueller_class.mueller_layer())
            self.muellers.loc[i, 'M_rif_dw'] = np.array(mueller_class.mueller_reflection())
            self.muellers.loc[i, 'M_tra_dw'] = np.array(mueller_class.mueller_transmission())
            
            #la matrice di Mueller per la trasmissione per ragggi diretti dal basso verso l'alto è
            #diversa da quella per raggi diretti dall'alto verso il basso
            #quindi:
            #aggiorno theta
            theta_start = mueller_class.theta1
            #e me la ricalcolo
            mueller_class = MM(n1, n0, theta_start, wsuc, spessore)
            self.muellers.loc[i, 'M_tra_up'] = np.array(mueller_class.mueller_transmission())
            self.muellers.loc[i, 'M_rif_up'] = np.array(mueller_class.mueller_reflection())
            
            #print('rif alto a basso: ', mueller_class.mueller_reflection())
            #print('tras alto a basso: ', mueller_class.mueller_transmission())
            #print('rif basso a alto: ', mueller_class.mueller_reflection())
            #print('tras basso a alto: ', mueller_class.mueller_transmission())
            #print()
           
        for i in range(self.campione.strati+1):
            
            #print('elemento 00: ',self.muellers.loc[i, 'M_rif'][0,0])
            #normalizzo le matrici di Mueller
            M_mat_nr = normalize(self.muellers.loc[i, 'M_mat'])
            M_rif_dw_nr = normalize(self.muellers.loc[i, 'M_rif_dw'])
            M_rif_up_nr = normalize(self.muellers.loc[i, 'M_rif_up'])
            M_tra_dw_nr = normalize(self.muellers.loc[i, 'M_tra_dw'])
            M_tra_up_nr = normalize(self.muellers.loc[i, 'M_tra_up'])
            
            #calcolo matrici di covarianza
            H_mat = covariance_matrix(M_mat_nr)
            H_rif_dw = covariance_matrix(M_rif_dw_nr)
            H_rif_up = covariance_matrix(M_rif_up_nr)
            H_tra_dw = covariance_matrix(M_tra_dw_nr)
            H_tra_up = covariance_matrix(M_tra_up_nr)
            
            #calcolo i vettori di covarianza
            h_mat = covariance_vector(H_mat)
            h_rif_dw = covariance_vector(H_rif_dw)
            h_rif_up = covariance_vector(H_rif_up)
            h_tra_dw = covariance_vector(H_tra_dw)
            h_tra_up = covariance_vector(H_tra_up)
            
            #print(H_rif) #NOTA: h_mat è autovettore di una matrice multipla dell'identità
            #print()
            
            self.biquaternions.loc[i, 'h_mat'] = Quaternion(h_mat[0], h_mat[1]*I, h_mat[2]*I, h_mat[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_rif_dw'] = Quaternion(h_rif_dw[0], h_rif_dw[1]*I, h_rif_dw[2]*I, h_rif_dw[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_rif_up'] = Quaternion(h_rif_up[0], h_rif_up[1]*I, h_rif_up[2]*I, h_rif_up[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_tra_dw'] = Quaternion(h_tra_dw[0], h_tra_dw[1]*I, h_tra_dw[2]*I, h_tra_dw[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_tra_up'] = Quaternion(h_tra_up[0], h_tra_up[1]*I, h_tra_up[2]*I, h_tra_up[3]*I, real_field = False)


#Normalizzazione della matrice di Mueller       
def normalize(Mueller_mat):
        M_00 = Mueller_mat[0,0]
        normalized = Mueller_mat*(1./M_00)
        return normalized


#Prodotto scalare fra biquaterioni
def scalar_prod(q1, q2):
	return q1.a*q2.a + q1.b*q2.b + q1.c*q2.c + q1.d*q2.d


#Coniugato di un quaternione
def conjugate(h):

	tau   =  np.conj( h.a )
	alfa  =  np.conj(-h.b*I)
	beta  =  np.conj(-h.c*I)
	gamma =  np.conj(-h.d*I)
	
	return Quaternion( tau, alfa*I, beta*I, gamma*I )

     
#Moltiplicazione di quaternioni per ottenere h, hdaga
def multiplication(arr):
	h = arr[0]
	
	if len(arr)==1:
		return [ h, conjugate(h) ]

	else:
		for i in range(1, len(arr)):
			h.mul(arr[i])
		return [ h, conjugate(h) ]


#Calcolo le grandezze ellissometriche Psi, delta
#a partire dai quaternioni arr di un singolo raggio 
def grandell(arr, s, svfinal=False, quatfinal=False):

	h, hdaga = multiplication(arr)

	shdaga = s.mul(hdaga)
	sfin	= h.mul(shdaga)
	rfin = stokes_vector( complex(sfin.a), complex(sfin.b*(-1j)), complex(sfin.c*(-1j)), complex(sfin.d*(-1j)) )  #vettore di Stokes finale
	
	hs = h.mul(s)			#prodotto hs tra quaternioni
	shs = scalar_prod(h, hs)	#prodotto scalare s.hs

	psi = rfin.ellipsometric_Psi().real
	delta = cmath.phase(shs).real

	if svfinal == False and quatfinal == False:
		return [psi, delta]

	elif svfinal == True and quatfinal == False:
		return [psi, delta, rfin]

	elif svfinal == False and quatfinal == True:
		return [psi, delta, sfin]

	else:
		return [psi, delta, rfin, sfin]


