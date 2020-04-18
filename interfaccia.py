'''	Funzione interfaccia 
	
	Input   		[theta0, nc0, nc1]
    theta0 	    	angolo di incidenza (rispetto alla normale)
    nc0 			indice di rifrazione del mezzo dal quale incide il raggio
    nc1				indice di rifrazione del primo mezzo sul quale incide il raggio

	Output			[theta1, tau_sigma, rho_sigma, tau_pi, rho_pi]
	theta1			angolo del fascio rifratto
	rho_sigma		riflettività luce polarizzata ortogonalmente
	rho_pi			riflettività luce polarizzata parallelamente
	tau_sigma		trasmissività luce polarizzata ortogonalmente
	tau_pi			trasmissività luce polarizzata parallelamente

    Riflettività e trasmissività corrispondono all'angolo di incidenza theta0  per l'interfaccia 
	tra due mezzi con indice di rifrazione nc0 e nc1. x_pi e x_sigma indicano i valori 
	per la luce polarizzata rispettivamente parallelamente e ortogonalmente al piano di incidenza '''

import math

def interfaccia(theta0, nc0, nc1):

	theta1 = math.asin( nc0*math.sin(theta0)/nc1 )  #angolo di rifrazione

	r = nc0/nc1
	a = nc1/nc0*math.cos(theta1)/math.cos(theta0)
	b = (r**2)*a

	rho_sigma = (1-a)/(1+a);  #riflessività ortogonale
	rho_pi = (1-b)/(1+b);	  #riflessività parallela
	
	tau_sigma = 2/(1+a)		  #trasmissività ortogonale
	tau_pi = 2*r/(1+b);		  #trasmissività parallela

	return [theta1, tau_sigma, rho_sigma, tau_pi, rho_pi]
