import numpy as np
import array as arr
import math
import cmath
from campione import campione
from sorgente import sorgente
import matplotlib.pyplot as plt

'''
Script per simulare la propagazione e l'interferenza di in un raggio solo 
che si splitta alla prima interfaccia in due nuovi raggi: il primo viene riflesso
dalla prima interfaccia, il secondo viene trasmesso dalla prima interfaccia, 
propaga in uno strato di mezzo, viene riflesso dalla seconda interfaccia,
propaga nuovamente nel mezzo, e infine viene trasmesso dalla prima interfaccia.
I due raggia prodotti vengono quindi fatti interferire.
'''

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)
pi = math.pi

############################# INIZIALIZZAZIONE ##############################
#inizializzo il campione
nc0 = 1
nc1 = 1.5

camp = campione(nc0, [nc1])
   
#inizializzo la sorgente	
sorg = sorgente(lenght=632.7193201669578)

theta = [0, 1.46] #0.7242567, 1.2, 1.435, pi/2]

nval = 2
I = []
Q = []
U = []
V = []

for i in range(nval):
    print('############## theta = ', theta[i],' ###############')
    
    theta0 = theta[i]
    
    theta1 = cmath.asin( camp.nc[0]*cmath.sin(theta0)/camp.nc[1] )
    
    #coefficienti di fresnel della prima interfaccia
    r = camp.nc[0]/camp.nc[1]
    a = (camp.nc[1]/camp.nc[0])*(cmath.cos(theta1)/cmath.cos(theta0))
    b = (r**2)*a

    rho1_sigma = (1-a)/(1+a);  #riflessività ortogonale
    rho1_pi = -(1-b)/(1+b);	  #riflessività parallela
		
    tau1_sigma = 2/(1+a)		  #trasmissività ortogonale
    tau1_pi = 2*r/(1+b);		  #trasmissività parallela
    
    '''print('rho_s = ',rho1_sigma)
    print('tau_s = ',tau1_sigma)'''
    #coefficiente di fresnel di trasmissione up
    r = camp.nc[1]/camp.nc[0]
    a = camp.nc[0]/camp.nc[1]*cmath.cos(theta0)/cmath.cos(theta1)
    b = (r**2)*a
    
    rho1_sigma_up = (1-a)/(1+a);  #riflessività ortogonale
    rho1_pi_up = -(1-b)/(1+b);	  #riflessività parallela
    
    tau1_sigma_up = 2/(1+a)		  #trasmissività ortogonale
    tau1_pi_up = 2*r/(1+b);		  #trasmissività parallela
    print('tau_s_up = ',tau1_sigma_up)
    ########################## PROPAGAZIONE DEI RAGGI ###########################

    #inizializzazione del raggio iniziale
    y_sigma = 1#0.18014218115842623#
    y_pi = 0 #polarizzazione lineare 45°

    intensità_in = nc0*np.cos(theta0)*y_sigma**2
    print('intensità del raggio iniziale = ', nc0*np.cos(theta0)*y_sigma**2)
    print()    
    

    #il PRIMO raggio viene riflesso dalla prima interfaccia
    y_sigma_fin1 = y_sigma*rho1_sigma
    y_pi_fin1 = y_pi*rho1_pi

    intensità_rif1 = nc0*np.cos(theta0)*y_sigma_fin1**2
    print('intensità del raggio subito riflesso = ', intensità_rif1)
    #print('campo componente s = ', y_sigma_fin1)
    #print('campo componente p = ', y_pi_fin1)
    print()
    #il SECONDO, viene trasmesso...
    y_sigma_int2 = y_sigma*tau1_sigma
    intensità_tra1 = nc1*abs(y_sigma_int2)**2*np.cos(theta1)
    print('intensità dopo prima trasmissione = ', intensità_tra1)
    
    y_pi_int2 = y_pi*tau1_pi  
    
    fattore = camp.nc[1]/camp.nc[0]*cmath.cos(theta1)/math.cos(theta0)#cmath.sqrt(1-(camp.nc[0]/camp.nc[1])**2*(math.sin(theta0))**2)
    tra_S = fattore*(abs(y_sigma_int2))**2/(abs(y_sigma))**2 #fattore*y_sigma_int2**2/y_sigma**2#
    rif_S = (abs(y_sigma_fin1))**2/(abs(y_sigma))**2 #y_sigma_fin1**2/y_sigma**2#
    
    #print('coefficiente di trasmissione onda s = ', tra_S)
    #print('coefficiente di riflessione onda s = ', rif_S)
    
    #print('somma = ', tra_S + rif_S)
    #print('campo componente s = ', y_sigma_int2) 
    print('somme delle intensità = ', intensità_tra1+intensità_rif1)
    print()
    
    phase = (sorg.wsuc*(0.00002)/np.cos(theta1))*nc1
    y_sigma_int2 *= np.exp(1j*phase)    
    intensità = nc1*abs(y_sigma_int2)**2*np.cos(theta1)
    print('intensità dopo prima propagazione = ', intensità)

    y_sigma_ausiliario = y_sigma_int2
    y_sigma_int2 *= -rho1_sigma #riflessione su interfaccia da 1.5 a 1 (cambio segno)
    intensità_rif_sub = nc1*abs(y_sigma_int2)**2*np.cos(theta1)
    intensità_tra_sub = nc0*abs(y_sigma_ausiliario*tau1_sigma_up)**2*np.cos(theta0)    
    print('intensità dopo riflessione sul substrato = ', intensità_rif_sub)  
    rif_S = (abs(y_sigma_int2))**2/(abs(y_sigma_ausiliario))**2 #fattore*y_sigma_int2**2/y_sigma**2#
    #print('coefficiente di riflessione onda s = ', rif_S)
    #print('campo componente s = ', y_sigma_int2)
    print('somme delle intensità = ', intensità_rif_sub+intensità_rif1+intensità_tra_sub)
    print()
    
    y_sigma_int2 *= np.exp(1j*phase)
    intensità = nc1*abs(y_sigma_int2)**2*np.cos(theta1)
    print('intensità dopo seconda propagazione = ', intensità)
    
    y_sigma_ausiliario = y_sigma_int2     
    y_sigma_int2 *= tau1_sigma_up
    intensità = nc0*abs(y_sigma_int2)**2*np.cos(theta0)
    intensità_rif_ulteriore = nc1*np.cos(theta1)*abs(y_sigma_ausiliario*rho1_sigma_up)**2
    print('intensità dopo trasmissione finale = ', intensità) 
    print('somme delle intensità = ', intensità+intensità_rif1+intensità_tra_sub+intensità_rif_ulteriore)
    #print('x_sigma[x] = ', y_sigma_int2) 
   
    fattore = camp.nc[0]/camp.nc[1]*cmath.cos(theta0)/cmath.cos(theta1)
    tra_S = fattore*(abs(y_sigma_int2))**2/(abs(y_sigma_ausiliario))**2    
    #print('coefficiente di trasmissione onda s = ', tra_S) 
    #print('campo componente s = ', y_sigma_int2)

    print()
    print('Raggio iniziale:') ###########
    
    I = (abs(y_sigma))**2+(abs(y_pi))**2
    Q = (abs(y_sigma))**2-(abs(y_pi))**2
    U = 2*(y_sigma*np.conj(y_pi)).real
    V = -2*(y_sigma*np.conj(y_pi)).imag
                        
    #print('riflettività = ', i_/complex(self.s_quaternions.loc[k, 's'].a))                        
                        
    print() ############
    print('I = ', I) ############
    print('Q = ', Q) ############
    print('U = ', U) ############
    print('V = ', V) ############
    print()

    print('Raggio riflesso:') ###########
    
    I = (abs(y_sigma_fin1))**2+(abs(y_pi_fin1))**2
    Q = (abs(y_sigma_fin1))**2-(abs(y_pi_fin1))**2
    U = 2*(y_sigma_fin1*np.conj(y_pi_fin1)).real
    V = -2*(y_sigma_fin1*np.conj(y_pi_fin1)).imag
                        
    #print('riflettività = ', i_/complex(self.s_quaternions.loc[k, 's'].a))                        
                        
    print() ############
    print('I = ', I) ############
    print('Q = ', Q) ############
    print('U = ', U) ############
    print('V = ', V) ############
    print()
 
    print('Raggio trasmesso:') ###########
    
    I = (abs(y_sigma_int2))**2+(abs(y_pi_int2))**2
    Q = (abs(y_sigma_int2))**2-(abs(y_pi_int2))**2
    U = 2*(y_sigma_int2*np.conj(y_pi_int2)).real
    V = -2*(y_sigma_int2*np.conj(y_pi_int2)).imag
                        
    #print('riflettività = ', i_/complex(self.s_quaternions.loc[k, 's'].a))                        
                        
    print() ############
    print('I = ', I) ############
    print('Q = ', Q) ############
    print('U = ', U) ############
    print('V = ', V) ############
    print()       

    #print('theta1 = ', theta1)
    
    print('interferenza: ')
    r_fin_sigma = y_sigma_fin1 + y_sigma_int2
    print('componenti s da sommare: ', y_sigma_fin1, y_sigma_int2)
    intensità_interf = nc0*np.cos(theta0)*abs(r_fin_sigma)**2
    print('intensità: ', nc0*abs(y_sigma_fin1)**2*np.cos(theta0), nc0*abs(y_sigma_int2)**2*np.cos(theta0))
    
    print('campo componente s = ', r_fin_sigma)
    r_fin_pi = 0

    I = (abs(r_fin_sigma))**2+(abs(r_fin_pi))**2
    Q = (abs(r_fin_sigma))**2-(abs(r_fin_pi))**2
    U = 2*(r_fin_sigma*np.conj(r_fin_pi)).real
    V = -2*(r_fin_sigma*np.conj(r_fin_pi)).imag
    
    print() ############
    print('I = ', I) ############
    print('Q = ', Q) ############
    print('U = ', U) ############
    print('V = ', V) ############
    print()     
    
    intensità = nc0*abs(r_fin_sigma)**2*np.cos(theta0)
    print('intensità raggio uscente = ', intensità)
    print('somme delle intensità = ', intensità_interf+intensità_tra_sub)