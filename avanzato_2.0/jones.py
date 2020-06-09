import cmath
import numpy as np

class jones:
    
    def __init__(self, nc0, nc1, theta0, wsuc, spessore):
        
        self.nc0=nc0
        self.nc1=nc1
        self.theta0=theta0
        self.wsuc=wsuc
        self.spessore=spessore
        self.theta1 = cmath.asin( nc0*cmath.sin(theta0)/nc1 )
    
    def jones_propagation(self):
        
        phase = (self.wsuc)*(self.nc0)*(self.spessore/np.cos(self.theta0))
        
        #print('cos(theta0): ',np.cos(self.theta0))
        J_11 = np.exp(-1j*phase)
        J_12 = 0
        J_21 = 0
        J_22 = np.exp(-1j*phase)
        #print('phase: ', np.exp(-1j*phase))
        
        return calcola_parametri(J_11, J_12, J_21, J_22)
    
    def jones_transmission(self):
        nc0=self.nc0
        nc1=self.nc1
        theta0=self.theta0
        theta1 = self.theta1

        r = nc0/nc1
        a = nc1/nc0*cmath.cos(theta1)/cmath.cos(theta0)
        b = (r**2)*a
        
        tau_sigma = 2/(1+a);  #coefficiente di trasmissione onda s
        tau_pi = 2*r/(1+b);	  #coefficiente di trasmissione onda p

        J_11 = np.conjugate(tau_sigma)       
        J_12 = 0
        J_21 = 0
        J_22 = np.conjugate(tau_pi)        
        
        return calcola_parametri(J_11, J_12, J_21, J_22)
    
    def jones_reflection(self):
        nc0=self.nc0
        nc1=self.nc1
        theta0=self.theta0
        
        theta1 = self.theta1  #angolo di rifrazione
        r = nc0/nc1
        a = nc1/nc0*cmath.cos(theta1)/cmath.cos(theta0)
        b = (r**2)*a
        
        rho_sigma = (1-a)/(1+a);  #riflessività ortogonale
        rho_pi = (1-b)/(1+b);	  #riflessività parallela

        J_11 = np.conjugate(rho_sigma)       
        J_12 = 0
        J_21 = 0
        J_22 = np.conjugate(rho_pi)           
        
        return calcola_parametri(J_11, J_12, J_21, J_22)
    
    
def calcola_parametri(J_11, J_12, J_21, J_22):

    tau = (J_11 + J_22)/2.
    alpha = (J_11 - J_22)/2.
    beta = (J_12 + J_21)/2.
    gamma = (J_21 - J_12)/(2.*1j)
    
    return tau, alpha, beta, gamma
