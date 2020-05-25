import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from sympy import I
from sympy.algebras.quaternion import Quaternion
from stokes import stokes_vector
from sorgente import sorgente
from campione import campione
from interazione import *

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)

nc = 1

#definisco la sorgente
sorg = sorgente(1.95954488, 1, 0.78539163) #per ora non uso i due argomenti a dx

omega=2*math.pi*(sorg.energia)/h;
wsuc=omega/c;
lunghezza_donda = 2*math.pi*c/(nc*omega)

#print('freq ang: ', wsuc)
 
r = stokes_vector()			#inizializzo oggetto vettore di stokes

nvalues = 1
theta_vect = np.linspace(0, pi/2, nvalues)

spessore1 = 0.00001
spessore2 = 0.00001 + lunghezza_donda/2.

theta = 0  
    
#matrice di Mueller 1
phase1 = wsuc*nc*(spessore1/np.cos(theta))
#print('phase1: ',phase1)
attenuazione = abs(np.exp(-1j*phase1)) ######
#print('attenuazione: ', attenuazione)
M1 = np.zeros( (4,4) , dtype = complex)
M1[0,0] = attenuazione
M1[1,1] = attenuazione
M1[2,2] = attenuazione
M1[3,3] = attenuazione

#matrice di Mueller 2
phase2 = wsuc*nc*(spessore2/np.cos(theta))
#print('phase2: ',phase2)
#print('differenza phase2-phase1: ',phase2-phase1)
attenuazione = abs(np.exp(-1j*phase2)) ######
M2 = np.zeros( (4,4) , dtype = complex)
M2[0,0] = attenuazione
M2[1,1] = attenuazione
M2[2,2] = attenuazione
M2[3,3] = attenuazione

#######################################
#Calcolo col formalismo dei quaternioni	
r.generic_polarization(1, 1)	#vettore di Stokes incidente1, polarizzazione lineare 
s = Quaternion( r.I(), r.Q()*I, r.U()*I, r.V()*I )	#quaternione corrispondente al vett Stoke

#controllo
#print("vettore di Stokes iniziale")
#print(r.parameters[0])
#print(r.parameters[1])
#print(r.parameters[2])
#print(r.parameters[3])
    
#quaternioni dei materiali
M1_ = normalize(M1)
M2_ = normalize(M2)

H1 = covariance_matrix(M1_)
H2 = covariance_matrix(M2_)

h1 = covariance_vector(H1)*np.exp(-1j*phase1)
h2 = covariance_vector(H2)*np.exp(-1j*phase2)
    
h1q = Quaternion(h1[0], h1[1]*I, h1[2]*I, h1[3]*I, real_field = False)
h2q = Quaternion(h2[0], h2[1]*I, h2[2]*I, h2[3]*I, real_field = False)

    
#rinormalizzo
h1q *= cmath.sqrt(M1[0,0]) 
h2q *= cmath.sqrt(M2[0,0]) 
    
#print('normalizzazioni: ')
#print(cmath.sqrt(M1[0,0]))
#print(cmath.sqrt(M2[0,0]))
    
#quaternioni finali corrispondenti ai singoli materiali
h1qf, h1qfdaga = multiplication([h1q])
h2qf, h2qfdaga = multiplication([h2q])

print('phase: ', np.exp(-1j*phase1))
print('h1qf: ', h1qf)

print('phase: ', np.exp(-1j*phase2))
print('h2qf: ', h2qf)  

#quaternione finale del raggio attraverso il singolo materiale
#h, hdaga = multiplication(arr) CONTROLLA NEL PROGRAMMA PRINCIPALE
sh1qfdaga = s.mul(h1qfdaga)
s1fin	= h1qf.mul(sh1qfdaga)

#print(h1qfdaga)
#print(h1qf)
    
sh2qfdaga = s.mul(h2qfdaga)
s2fin	= h2qf.mul(sh1qfdaga)

'''   
#MODIFICA IMPORTANTE
h1qfs = h1qf.mul(s1fin)			#prodotto hs tra quaternioni
    
#print(h1qfdaga)
#print(s)
#print(sh1qfdaga)
    
shs1 = scalar_prod(s, h1qfs)	#prodotto scalare s.hs CONTROLLA!
    
#print(shs1)
    
h2qfs = h2qf.mul(s2fin)			#prodotto hs tra quaternioni
shs2 = scalar_prod(s, h2qfs)	#prodotto scalare s.hs
    
#print(shs2)

htot = h1qf*cmath.exp(1j*shs1) + h2qf*cmath.exp(1j*shs2) #CI VUOLE -1j?
'''
      
htot, htotdaga = multiplication([h1qf + h2qf])
sfin, sfindaga = multiplication([htot, s, htotdaga])
#print(sfin)   
    
rfin = stokes_vector( complex(sfin.a), complex(sfin.b*(-1j)), complex(sfin.c*(-1j)), complex(sfin.d*(-1j)) )  #vettore di Stokes finale

#vettore si Stokes finale
S = rfin.I().real
Q = rfin.Q().real
U = rfin.U().real
V = rfin.V().real 

print()
print("vettore di Stokes finale: ")
print(S)
print(Q)
print(U)
print(V)
