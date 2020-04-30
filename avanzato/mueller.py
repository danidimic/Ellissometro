import math
import cmath
import numpy as np

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)
pi=math.pi

class MM:
    
    def __init__(self,nc0,nc1,theta0):
        self.nc0=nc0
        self.nc1=nc1
        self.theta0=theta0
        
    def matrix(self):
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
   
    
m = MM(1,1.5-0.000002j,pi/4)
print(m.matrix())

        
