import math
import cmath
import numpy as np

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)
pi=math.pi

class MM:
    
    #   
    #      *
    #        *
    #          *
    #  __________*________________________________
    #              \         :         /
    #   S           \        :        /
    #   P            \theta0 :       /
    #   E             \      :      /
    #   S              \     :     /
    #   S               \    :    /           nc0
    #   O                \   :   /
    #   R                 \  :  /
    #   E                  \ : /
    #  _____________________\:/___________________
    #                        :*
    #                        :  *
    #                        :    *
    #                        :      *         nc1
    #                        : theta1 *
    #                        :          *

    
    def __init__(self, nc0, nc1, theta0, wsuc, spessore):
        self.nc0=nc0
        self.nc1=nc1
        self.theta0=theta0
        self.wsuc=wsuc
        self.spessore=spessore
        

    def mueller_reflection(self):
        nc0=self.nc0
        nc1=self.nc1
        theta0=self.theta0
        
        theta1 = math.asin( nc0.real*math.sin(theta0)/nc1.real )  #angolo di rifrazione
        r = nc0/nc1
        a = nc1/nc0*math.cos(theta1)/math.cos(theta0)
        b = (r**2)*a
        
        rho_sigma = (1-a)/(1+a);  #riflessività ortogonale
        rho_pi = (1-b)/(1+b);	  #riflessività parallela
        
        Rp= abs(rho_sigma)**2
        Rs= abs(rho_pi)**2
        
        # creo la matrice di Mueller
        Rhs=(Rp+Rs)/2
        Rhd=(Rp-Rs)/2
        Rps = (rho_pi * rho_sigma.conjugate())
        M = np.zeros( (4,4) )
        M[0,0] = Rhs
        M[1,1] = Rhs
        M[0,1] = Rhd
        M[1,0] = Rhd
        M[2,2] = Rps.real
        M[3,3] = Rps.real
        M[2,3] = -Rps.imag
        M[3,2] = Rps.imag
        return M
   

    def mueller_transmission(self):
        nc0=self.nc0
        nc1=self.nc1
        theta0=self.theta0
        
        theta1 = math.asin( nc0.real*math.sin(theta0)/nc1.real )
        r = nc0/nc1
        a = nc1/nc0*math.cos(theta1)/math.cos(theta0)
        b = (r**2)*a
        
        tau_sigma = 2/(1+a);  #coefficiente di trasmissione onda s
        tau_pi = 2*r/(1+b);	  #coefficiente di trasmissione onda p
        
        Ts= abs(tau_sigma)**2 #riflettività onda s
        Tp= abs(tau_pi)**2    #riflettività onda s
        
        # creo la matrice di Mueller
        Ths=(Ts+Tp)/2
        Thd=(Tp-Ts)/2
        Tps = (tau_pi * tau_sigma.conjugate())
        M = np.zeros( (4,4) )
        M[0,0] = Ths
        M[1,1] = Ths
        M[0,1] = Thd
        M[1,0] = Thd
        M[2,2] = Tps.real
        M[3,3] = Tps.real
        M[2,3] = -Tps.imag
        M[3,2] = Tps.imag
        return M
    

    def mueller_layer(self):
        #ATTENZIONE: uso theta0 e n0 in modo da poter considerare anche gli 
        #effetti del primo mezzo (aria). Eventualmente si può cambiare.
        
        phase = self.wsuc*self.nc0*(self.spessore/np.cos(self.theta0))
        attenuazione = abs(np.exp(-1j*phase))
        M = np.zeros( (4,4) )
        M[0,0] = attenuazione
        M[1,1] = attenuazione
        M[2,2] = attenuazione
        M[3,3] = attenuazione
        
        return M
