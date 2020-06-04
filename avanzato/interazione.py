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
from jones import jones


c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)

class interazione:
    
    def __init__(self, precisione, campione, sorgente, stokes_iniz):
        
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
        
        
        #lista con i quaternioni che dovrò usare per calcolare la differenze di fase prima di
        #effettuare la combinazione lineare dei cammini (è GIUSTA L'IDEA?) e dei quaternioni 
        #dei rispettivi cammini
        righe = []
        col = ['s_fin', 'h_fin']
        self.s_h_fin = pd.DataFrame(0, righe, col)
        self.s_h_fin = self.s_h_fin.astype(object)        
        
        
        #inizializzo il vettore con i quaternioni dei raggi durante la propagazione,
        #al SOLO RAGGIO INIZIALE (DA FARE)
        righe = []
        col = ['s', 'provenienza', 'arrivo', 'h_parziale']
        self.s_quaternions = pd.DataFrame(0, righe, col)
        self.s_quaternions = self.s_quaternions.astype(object)
        
        #inizializzazione del quaterione del raggio iniziale
        s_iniz = Quaternion(stokes_iniz.I(), stokes_iniz.Q()*I, stokes_iniz.U()*I, stokes_iniz.V()*I)
        self.s_quaternions = self.s_quaternions.append({'s': s_iniz, "provenienza": 0, "arrivo": 1, "h_partial": 1}, ignore_index=True)

        self.nraggi = 1
        


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
            phase_aquired, M = mueller_class.mueller_layer()
            self.muellers.loc[i, 'M_mat'] = np.array(M)
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
            
        for i in range(self.campione.strati+1):
            
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
            h_mat = covariance_vector(H_mat)*phase_aquired
            h_rif_dw = covariance_vector(H_rif_dw)
            h_rif_up = covariance_vector(H_rif_up)
            h_tra_dw = covariance_vector(H_tra_dw)
            h_tra_up = covariance_vector(H_tra_up)
            
            self.biquaternions.loc[i, 'h_mat']    = Quaternion(h_mat[0], h_mat[1]*I, h_mat[2]*I, h_mat[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_rif_dw'] = Quaternion(h_rif_dw[0], h_rif_dw[1]*I, h_rif_dw[2]*I, h_rif_dw[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_rif_up'] = Quaternion(h_rif_up[0], h_rif_up[1]*I, h_rif_up[2]*I, h_rif_up[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_tra_dw'] = Quaternion(h_tra_dw[0], h_tra_dw[1]*I, h_tra_dw[2]*I, h_tra_dw[3]*I, real_field = False)
            self.biquaternions.loc[i, 'h_tra_up'] = Quaternion(h_tra_up[0], h_tra_up[1]*I, h_tra_up[2]*I, h_tra_up[3]*I, real_field = False)
            #Normalizzo subito tutti i quaternioni per non doverlo fare dopo
            self.biquaternions.loc[i, 'h_mat']    *= np.sqrt( self.muellers.loc[i, 'M_mat'][0,0] )
            self.biquaternions.loc[i, 'h_rif_dw'] *= np.sqrt( self.muellers.loc[i, 'M_rif_dw'][0,0] )
            self.biquaternions.loc[i, 'h_rif_up'] *= np.sqrt( self.muellers.loc[i, 'M_rif_up'][0,0] )
            self.biquaternions.loc[i, 'h_tra_dw'] *= np.sqrt( self.muellers.loc[i, 'M_tra_dw'][0,0] )
            self.biquaternions.loc[i, 'h_tra_up'] *= np.sqrt( self.muellers.loc[i, 'M_tra_up'][0,0] )
            
            
    def materials_to_jones(self, theta0):
        
        '''
        Raccoglie nel Data Frame, con indici di riga riferiti agli strati (0 
        è il vuoto): quaternione ottenuto da matrice di Jones relativa alla propagazione 
        nel materiale (h_mat), quaternioni ottenuti da matrici di Jones relative alla 
        riflessione sull'interfaccia (in basso) (h_rif_up/_dw), quaternioni ottenuti da 
        matrici di Jones relative alla trasmissione attraverso l'interfaccia (in basso) 
        (h_tra_up/dw). NOTA BENE: in questa versione del programma considero il mezzo 
        isotropo, per cui gli effetti della propagazione nel mezzo non dipendono dal 
        segno dell'angolo, né dal fatto che la luce propaghi dal basso verso l'alto o viceversa.
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

            jones_class = jones(n0, n1, theta_start, wsuc, spessore)
        
            tau_mat, alfa_mat, beta_mat, gamma_mat = jones_class.jones_propagation()
            self.biquaternions.loc[i, 'h_mat'] = Quaternion(tau_mat, alfa_mat*I, beta_mat*I, gamma_mat*I)
            
            tau_rif_dw, alfa_rif_dw, beta_rif_dw, gamma_rif_dw = jones_class.jones_reflection()
            self.biquaternions.loc[i, 'h_rif_dw'] = Quaternion(tau_rif_dw, alfa_rif_dw*I, beta_rif_dw*I, gamma_rif_dw*I)
            
            tau_tra_dw, alfa_tra_dw, beta_tra_dw, gamma_tra_dw = jones_class.jones_transmission()
            self.biquaternions.loc[i, 'h_tra_dw'] = Quaternion(tau_tra_dw, alfa_tra_dw*I, beta_tra_dw*I, gamma_tra_dw*I)
            
            theta_start = jones_class.theta1
            jones_class = jones(n1, n0, theta_start, wsuc, spessore)
            
            tau_rif_up, alfa_rif_up, beta_rif_up, gamma_rif_up = jones_class.jones_reflection()
            self.biquaternions.loc[i, 'h_rif_up'] = Quaternion(tau_rif_up, alfa_rif_up*I, beta_rif_up*I, gamma_rif_up*I)

            tau_tra_up, alfa_tra_up, beta_tra_up, gamma_tra_up = jones_class.jones_transmission()
            self.biquaternions.loc[i, 'h_tra_up'] = Quaternion(tau_tra_up, alfa_tra_up*I, beta_tra_up*I, gamma_tra_up*I)
            

    def propagazione(self):

        for k in range(self.nraggi):        
            
            #raggio che considero
            s_ = self.s_quaternions.loc[k, 's']
            ii_ = self.s_quaternions.loc[k, 'provenienza']
            jj_ = self.s_quaternions.loc[k, 'arrivo']
            h_part = self.s_quaternions.loc[k, 'h_partial'] #quaternione relativo alla propagazione precedente
            
            print()
            print('Raggio ', k+1)
            print('  ii = ', ii_)
            print('  jj = ', jj_)
            print('  intensità = ', abs(complex(s_.a)))
            
            if abs(complex(s_.a)) > self.precisione: #CONTROLLA
                
                #tratto a parte il caso di raggi vicini alla superficie
                if ii_ == 0 and jj_ == 1:
                    
                    #####################viene riflesso, e il primo cammino finisce qui!#################
                    #non considero effetti di propagazione nell'aria (si può aggiungere, con h_mat indice 0)
                    h = self.biquaternions.loc[0, 'h_rif_dw']
                    hdaga = conjugate(h)
                    shdaga = s_.mul(hdaga)
                    sfin	= h.mul(shdaga)
                    
                    self.s_h_fin = self.s_h_fin.append({'s_fin': sfin, 'h_fin': h}, ignore_index=True)
                    #####################################################################################
                    
                    #######################viene trasmesso e produce un nuovo cammino####################
                    #trasmissione
                    h_tra = self.biquaternions.loc[0, 'h_tra_dw']
                    
                    #propagazione nel mezzo
                    h_mat = self.biquaternions.loc[1, 'h_mat']
                    
                    h, hdaga = multiplication([h_mat, h_tra])                    
                    shdaga = s_.mul(hdaga)
                    sfin = h.mul(shdaga)
                    
                    #questo lo aggiungo al DF dei raggi
                    self.s_quaternions = self.s_quaternions.append({'s': sfin, "provenienza": 1, "arrivo": 2, "h_partial": h}, ignore_index=True)
                    #print('s, h, final: ',self.s_h_fin)
                    #####################################################################################

                elif ii_ == 1 and jj_ == 0:
                    
                    ######################viene trasmesso, e il cammino finisce qui!#####################
                    #non considero effetti di propagazione nell'aria (si può aggiungere, con h_mat indice 0)
                    
                    h, hdaga = multiplication([self.biquaternions.loc[0, 'h_tra_up'], h_part]) #aggiungo l'ultimo pezzo di trasmissione

                    shdaga = s_.mul(hdaga)
                    sfin	= h.mul(shdaga)
                    
                    self.s_h_fin = self.s_h_fin.append({'s_fin': sfin, 'h_fin': h}, ignore_index=True)
                    #print('s, h, final: ',self.s_h_fin)

                    #####################################################################################                    
                    
                    ###################viene riflesso in basso e produce un nuovo cammino################
                    #riflessione
                    h_rif = self.biquaternions.loc[0, 'h_rif_up']
                    
                    #propagazione nel mezzo
                    h_mat = self.biquaternions.loc[1, 'h_mat']
                    
                    h, hdaga = multiplication([ h_mat, h_rif, h_part])      
                    shdaga = s_.mul(hdaga)
                    sfin = h.mul(shdaga)
                    
                    #questo lo aggiungo al DF dei raggi
                    self.s_quaternions = self.s_quaternions.append({'s': sfin, "provenienza": 1, "arrivo": 2, "h_partial": h}, ignore_index=True)
                    #####################################################################################                    

                else:
                    
                    #raggi che propagano verso il basso
                    if ii_ < jj_:
                        
                        if jj_ < self.campione.strati+1:
                            #############################produce un raggio trasmesso#############################
                            #trasmissione
                            h_tra = self.biquaternions.loc[ii_, 'h_tra_dw']
                            
                            #propagazione nel mezzo jj_
                            h_mat = self.biquaternions.loc[jj_, 'h_mat']
                            
                            h, hdaga = multiplication([h_mat, h_tra, h_part])      
                            shdaga = s_.mul(hdaga)
                            sfin = h.mul(shdaga)
                        
                            #lo aggiungo al DF dei raggi
                            self.s_quaternions = self.s_quaternions.append({'s': sfin, "provenienza": ii_ + 1, "arrivo": jj_ + 1, "h_partial": h}, ignore_index=True)
                        
                            #####################################################################################                    

                        #############################produce un raggio riflesso##############################
                        #riflessione
                        h_rif = self.biquaternions.loc[ii_, 'h_rif_dw']
                    
                        #propagazione nel mezzo ii_
                        h_mat = self.biquaternions.loc[ii_, 'h_mat']
                        
                        h, hdaga = multiplication([h_mat, h_rif, h_part])      
                        shdaga = s_.mul(hdaga)
                        sfin = h.mul(shdaga)
                        
                        #lo aggiungo al DF dei raggi
                        self.s_quaternions = self.s_quaternions.append({'s': sfin, "provenienza": ii_, "arrivo": jj_ - 2, "h_partial": h}, ignore_index=True)

                        #####################################################################################
                        
                    #raggi che propagano verso l'alto
                    elif ii_ > jj_:
                        
                        #############################produce un raggio trasmesso#############################
                        #trasmissione
                        h_tra = self.biquaternions.loc[jj_, 'h_tra_up']
                    
                        #propagazione nel mezzo jj_
                        h_mat = self.biquaternions.loc[jj_, 'h_mat']
                        
                        h, hdaga = multiplication([h_mat, h_tra, h_part])      
                        shdaga = s_.mul(hdaga)
                        sfin = h.mul(shdaga)
                        
                        #lo aggiungo al DF dei raggi
                        self.s_quaternions = self.s_quaternions.append({'s': sfin, "provenienza": ii_ - 1, "arrivo": jj_ - 1, "h_partial": h}, ignore_index=True)
                        
                        #####################################################################################                    

                        #############################produce un raggio riflesso##############################
                        #riflessione
                        h_rif = self.biquaternions.loc[jj_, 'h_rif_up']
                    
                        #propagazione nel mezzo ii_
                        h_mat = self.biquaternions.loc[ii_, 'h_mat']
                        
                        h, hdaga = multiplication([h_mat, h_rif, h_part])      
                        shdaga = s_.mul(hdaga)
                        sfin = h.mul(shdaga)
                        
                        #lo aggiungo al DF dei raggi
                        self.s_quaternions = self.s_quaternions.append({'s': sfin, "provenienza": ii_, "arrivo": jj_ + 2, "h_partial": h}, ignore_index=True)

                        #####################################################################################
                        
                        
        #elimino i raggi che ho usato nel ciclo precedente da s_quaternions, se no le riuso!
        self.s_quaternions = self.s_quaternions.drop(range(self.nraggi))
        self.s_quaternions = self.s_quaternions.reset_index(drop = True)
            
        self.nraggi = self.s_quaternions.shape[0]



    def interference(self):
        h_ = []
        for i in self.s_h_fin.index:
            #if i != 0:
            h_.append(self.s_h_fin['h_fin'][i])
        
        print('***********************************************')
        print('Interferenza fra i cammini: ')
        for i in h_:
            print(complex(i.a), ' + ', complex(i.b), 'i + ',complex(i.c), 'j + ',complex(i.d), 'k')
        
        htot = Quaternion(0, 0, 0, 0)
        for i in h_:
            htot += i
        
        return htot

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
	
	l = len(arr)
	if l==1:
		h = arr[0]

	else:
		h = arr[l-2].mul(arr[l-1])
		for i in range(0, l-2):
			index = l-3-i
			h = arr[index].mul(h)

	return [ h, conjugate(h) ]


#Calcolo le grandezze ellissometriche Psi, delta
#a partire dai quaternioni arr di un singolo raggio 
def grandell(arr, s, ellipsometry=True, svfinal=False, quatfinal=False):

	h, hdaga = multiplication(arr)
	shdaga = s.mul(hdaga)
	sfin	= h.mul(shdaga)
	rfin = stokes_vector( complex(sfin.a), complex(sfin.b*(-1j)), complex(sfin.c*(-1j)), complex(sfin.d*(-1j)) )  #vettore di Stokes finale
	rin = stokes_vector( complex(s.a), complex(s.b*(-1j)), complex(s.c*(-1j)), complex(s.d*(-1j)) )  #vettore di Stokes iniziale
    
	results = []

	if(ellipsometry==True):
		'''
        psi = rfin.ellipsometric_Psi().real

        hs = h.mul(s)			#prodotto hs tra quaternioni
		shs = scalar_prod(s, hs)	#prodotto scalare s.hs
		delta = cmath.phase(shs).real
		'''
        
		rho = (np.tan(rin.alfa())*np.exp(1j*rin.delta()))/(np.tan(rfin.alfa())*np.exp(1j*rfin.delta()))
		psi = np.arctan(abs(rho))
		delta = cmath.phase(rho)
        
		if delta < 0:
			delta += 2*math.pi #CONTROLLA

		results.append(psi)
		results.append(delta)

	if(svfinal==True):
		results.append(rfin)

	if(quatfinal==True):
		results.append(sfin)

	if len(results)==1:
		results = results[0]
	
	return results
