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

############################# INIZIALIZZAZIONE ##############################
#lettura dei dati del campione
f = open("campione.txt", 'r')
cont = f.readlines()
for i in range(len(cont)):
	cont[i] = complex(cont[i])

#inizializzo il campione
nc0 = cont[0]
nc1 = cont[1]
spessore = cont[2]
nc2 = cont[3]

camp = campione(nc0, [nc1, spessore, nc2])
   
#lettura dei dati della sorgente
f = open("sorgente.txt", 'r')
cont = f.readlines()
for i in range(len(cont)):
	cont[i] = float(cont[i])

#inizializzo la sorgente	
sorg = sorgente(cont[0], cont[1], cont[2])

############################# INTERFACCE ####################################

passo = 0.01

Theta = []

nval = int(math.pi/(2*passo))

I = []
Q = []
U = []
V = []

for i in range(nval):
    
    theta0 = passo*i
    Theta.append(theta0)
    
    theta1 = cmath.asin( camp.nc[0]*cmath.sin(theta0)/camp.nc[1] )
    
    #coefficienti di fresnel della prima interfaccia
    r = camp.nc[0]/camp.nc[1]
    a = camp.nc[1]/camp.nc[0]*cmath.cos(theta1)/cmath.cos(theta0)
    b = (r**2)*a

    rho1_sigma = (1-a)/(1+a);  #riflessività ortogonale
    rho1_pi = (1-b)/(1+b);	  #riflessività parallela
		
    tau1_sigma = 2/(1+a)		  #trasmissività ortogonale
    tau1_pi = 2*r/(1+b);		  #trasmissività parallela


    theta2 = cmath.asin( camp.nc[1]*cmath.sin(theta1)/camp.nc[2] )

    #coefficienti di frasnel seconda interfaccia
    r = camp.nc[1]/camp.nc[2]
    a = camp.nc[2]/camp.nc[1]*cmath.cos(theta2)/cmath.cos(theta1)
    b = (r**2)*a

    rho2_sigma = (1-a)/(1+a);  #riflessività ortogonale
    rho2_pi = (1-b)/(1+b);	  #riflessività parallela
		
    tau2_sigma = 2/(1+a)		  #trasmissività ortogonale
    tau2_pi = 2*r/(1+b);		  #trasmissività parallela
    #ATTENZIONE: IN QUESTO MODO PARTE DELL'INTENSITA' VIENE PERSA NEL SUBSTRATO

    ########################## PROPAGAZIONE DEI RAGGI ###########################

    #inizializzazione del raggio iniziale
    y_sigma = 1
    y_pi = y_sigma*math.tan(sorg.psi_0)*cmath.exp(1j*sorg.delta_0)

    #il PRIMO raggio viene riflesso dalla prima interfaccia
    y_sigma_fin1 = y_sigma*rho1_sigma
    y_pi_fin1 = y_pi*rho1_pi


    #il SECONDO, viene trasmesso...
    y_sigma_int2 = y_sigma*tau1_sigma
    y_pi_int2 = y_pi*tau1_pi
    
    #...propaga nel mezzo...
    omega=2*math.pi*(sorg.energia)/h;
    wsuc=omega/c;

    phase = (-wsuc*(camp.spessori[0])/np.cos(theta1))*camp.nc[1]
    y_sigma_int2 = y_sigma_int2*cmath.exp(phase*1j)
    y_pi_int2 = y_pi_int2*cmath.exp(phase*1j)

    #...viene riflesso dalla seconda interfaccia...
    y_sigma_int2 = y_sigma_int2*rho2_sigma
    y_pi_int2 = y_pi_int2*rho2_pi

    #...propaga di nuovo nel mezzo...
    y_sigma_int2 = y_sigma_int2*cmath.exp(phase*1j)
    y_pi_int2 = y_pi_int2*cmath.exp(phase*1j)

    #...viene trasmesso dalla prima interfaccia
    y_sigma_fin2 = y_sigma_int2*tau1_sigma
    y_pi_fin2 = y_pi_int2*tau1_pi #CONTROLLA

    ################################ INTERFERENZA ###############################

    y_sigma = y_sigma_fin1+y_sigma_fin2
    y_pi = y_pi_fin1+y_pi_fin2
    
    delta1 = -cmath.phase( y_pi/y_sigma );

    I.append((abs(y_sigma))**2+(abs(y_pi))**2)
    Q.append((abs(y_sigma))**2-(abs(y_pi))**2)
    U.append(2*abs(y_sigma)*abs(y_pi)*np.cos(delta1))
    V.append(2*abs(y_sigma)*abs(y_pi)*np.sin(delta1))
    
#plot
#Grafico di I, Q, U, V
plt.title("Metodo standard")
plt.xlabel("angolo incidente [°]")
plt.ylabel("$I$ | $Q$ [°]")
plt.plot(Theta, I, label='$I$')
plt.plot(Theta, Q, label='$Q$')
plt.plot(Theta, U, label='$U$')
plt.plot(Theta, V, label='$V$')
plt.legend()
plt.grid(True)
plt.show()
