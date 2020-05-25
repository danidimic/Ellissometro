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
from jones import jones

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)

nc0 = 1
nc1 = 2

#definisco la sorgente
sorg = sorgente(1.95954488, 1, 0.78539163) #per ora non uso i due argomenti a dx

omega=2*math.pi*(sorg.energia)/h;
wsuc=omega/c;
lunghezza_donda = 2*math.pi*c/(nc0.real*omega)

#print('freq ang: ', wsuc)
 
r = stokes_vector()			#inizializzo oggetto vettore di stokes

#nvalues = 1
#theta_vect = np.linspace(0, pi/2, nvalues)

spessore1 = 0.00001
spessore2 = 0.00001 + lunghezza_donda/2.

theta = 0   
    
#matrice di Jones 1
jones1 = jones(nc0, nc1, theta, wsuc, spessore1)
    
#matrice di Jones 2
jones2 = jones(nc0, nc1, theta, wsuc, spessore2)

#######################################
#Calcolo col formalismo dei quaternioni	
r.generic_polarization(1, 1)	#vettore di Stokes incidente1, polarizzazione lineare 
s = Quaternion( r.I(), r.Q()*I, r.U()*I, r.V()*I )	#quaternione corrispondente al vett Stoke
    
#quaternioni dei materiali
tau1, alpha1, beta1, gamma1 = jones1.jones_propagation()

print('h1: ',tau1, alpha1, beta1, gamma1)

tau2, alpha2, beta2, gamma2 = jones2.jones_propagation()
    
print('h2: ',tau2, alpha2, beta2, gamma2)

h1q = Quaternion(tau1, alpha1*I, beta1*I, gamma1*I, real_field = False)
h2q = Quaternion(tau2, alpha2*I, beta2*I, gamma2*I, real_field = False)

print()
print('h1q: ', h1q)
print('h1q: ', h2q)

htot = h1q + h2q
    
htot, htotdaga = multiplication([htot])
    
sfin, sfindaga = multiplication([htot, s, htotdaga])   
    
rfin = stokes_vector( complex(sfin.a), complex(sfin.b*(-1j)), complex(sfin.c*(-1j)), complex(sfin.d*(-1j)) )  #vettore di Stokes finale

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