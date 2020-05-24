import cmath
import numpy as np

class jones:
    
    def __init__(self, nc0, nc1, theta0, wsuc, spessore):
        
        self.nc0=nc0
        self.nc1=nc1
        self.theta0=theta0
        self.wsuc=wsuc
        self.spessore=spessore
        self.theta1 = cmath.asin( nc0.real*cmath.sin(theta0)/nc1.real )
    
    def jones_propagation(self):
        
        phase = (self.wsuc)*(self.nc0)*(self.spessore/np.cos(self.theta0))
        
        print('cos(theta0): ',np.cos(self.theta0))
        J_11 = np.exp(-1j*phase)
        J_12 = 0
        J_21 = 0
        J_22 = np.exp(-1j*phase)
        
        return calcola_parametri(J_11, J_12, J_21, J_22)
        
def calcola_parametri(J_11, J_12, J_21, J_22):

    tau = (J_11 + J_22)/2
    alpha = (J_11 - J_22)/2
    beta = (J_12 + J_21)/2
    gamma = (J_21 - J_12)/(2*1j)
    
    return tau, alpha, beta, gamma