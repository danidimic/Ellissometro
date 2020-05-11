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
from sympy import I

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)

class interazione:
    
    def __init__(self, precisione, campione, sorgente):
        
        self.campione = campione
        self.sorgente = sorgente
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
            
            #normalizzo le matrici di Mueller
            M_mat_nr = normalize(self.muellers.loc[i, 'M_mat'])
            M_rif_nr = normalize(self.muellers.loc[i, 'M_rif'])
            M_tra_dw_nr = normalize(self.muellers.loc[i, 'M_tra_dw'])
            M_tra_up_nr = normalize(self.muellers.loc[i, 'M_tra_up'])
            
            #calcolo matrici di covarianza
            H_mat = covariance_matrix(M_mat_nr)
            H_rif = covariance_matrix(M_rif_nr)
            H_tra_dw = covariance_matrix(M_tra_dw_nr)
            H_tra_up = covariance_matrix(M_tra_up_nr)
            
            #calcolo i vettori di covarianza
            h_mat = covariance_vector(H_mat)
            h_rif = covariance_vector(H_rif)
            h_tra_dw = covariance_vector(H_tra_dw)
            h_tra_up = covariance_vector(H_tra_up)
            
            #print(H_rif) #NOTA: h_mat è autovettore di una matrice multipla dell'identità
            #print()
            
            self.biquaternions.loc[i, 'h_mat'] = Quaternion(h_mat[0], h_mat[1]*I, h_mat[2]*I, h_mat[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_rif'] = Quaternion(h_rif[0], h_rif[1]*I, h_rif[2]*I, h_rif[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_tra_dw'] = Quaternion(h_tra_dw[0], h_tra_dw[1]*I, h_tra_dw[2]*I, h_tra_dw[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_tra_up'] = Quaternion(h_tra_up[0], h_tra_up[1]*I, h_tra_up[2]*I, h_tra_up[3]*I, real_field = False)

            
def normalize(Mueller_mat):
        M_00 = Mueller_mat[0,0]
        normalized = Mueller_mat*(1./M_00)
        return normalized
    
