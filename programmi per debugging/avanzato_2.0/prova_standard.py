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

theta = [0, 1.2, 1.435, pi/2]

nval = 4
I = []
Q = []
U = []
V = []

for i in range(nval):
    print('##############theta = ', theta[i],'###############à')
    
    theta0 = theta[i]
    
    theta1 = cmath.asin( camp.nc[0]*cmath.sin(theta0)/camp.nc[1] )
    
    #coefficienti di fresnel della prima interfaccia
    r = camp.nc[0]/camp.nc[1]
    a = (camp.nc[1]/camp.nc[0])*(cmath.cos(theta1)/cmath.cos(theta0))
    b = (r**2)*a

    rho1_sigma = (1-a)/(1+a);  #riflessività ortogonale
    rho1_pi = (1-b)/(1+b);	  #riflessività parallela
		
    tau1_sigma = 2/(1+a)		  #trasmissività ortogonale
    tau1_pi = 2*r/(1+b);		  #trasmissività parallela
    print('rho_s = ',rho1_sigma)
    print('tau_s = ',tau1_sigma)
    #coefficiente di fresnel di trasmissione up
    r = camp.nc[1]/camp.nc[0]
    a = camp.nc[0]/camp.nc[1]*cmath.cos(theta0)/cmath.cos(theta1)
    b = (r**2)*a
    
    tau1_sigma_up = 2/(1+a)		  #trasmissività ortogonale
    tau1_pi_up = 2*r/(1+b);		  #trasmissività parallela

    ########################## PROPAGAZIONE DEI RAGGI ###########################

    #inizializzazione del raggio iniziale
    y_sigma = 1
    y_pi = 0 #polarizzazione lineare 45°

    #il PRIMO raggio viene riflesso dalla prima interfaccia
    y_sigma_fin1 = y_sigma*rho1_sigma
    y_pi_fin1 = y_pi*rho1_pi


    #il SECONDO, viene trasmesso...
    y_sigma_int2 = y_sigma*tau1_sigma
    y_pi_int2 = y_pi*tau1_pi    

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

    tra_S = camp.nc[1]/camp.nc[0]*y_sigma_int2**2/y_sigma**2#(abs(y_sigma_int2))**2/(abs(y_sigma))**2
    rif_S = y_sigma_fin1**2/y_sigma**2#(abs(y_sigma_fin1))**2/(abs(y_sigma))**2
    
    print('coefficiente di trasmissione onda s = ', tra_S)
    print('coefficiente di riflessione onda s = ', rif_S)
    print()