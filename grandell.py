# function [psi1,delta1,r_pi,r_sigma]=grandell(theta0)

# [psi1,delta1,r_pi,r_sigma]=grandell(theta0)
#
#    Calcola le grandezze ellissometriche psi1, delta1, r_pi, r_sigma per la luce riflessa 
#    dal sistema definito tramite ellmod.m.
#    {r_pi/r_sigma=tan(psi1)*exp(i*delta1), dove r_pi e r_sigma sono le riflettivita' per 
#    la luce polarizzata rispettivamente parallelamente e ortogonalmente al piano di incidenza}
#
#    theta0 =    angolo di incidenza (rispetto alla normale)
#
#    Variabili globali (cfr. ellmod.m, propagazione.m): strati (numero di strati=numero di interfacce-1), 
#    wsuc (pulsazione su c), theta (vettore -lungo (strati+2)- che contiene, per ogni mezzo, 
#    gli angoli tra il raggio e la normale alle interfacce), nc (indici di rifrazione dei 
#    mezzi), spessori (spessori degli strati; vale zero -irrilevante- per primo mezzo e 
#    substrato), tau_x/rho_x (trasmissivita' e riflettivita' di ciascuna interfaccia -vettore 
#    lungo (strati+1)), somma_x (variabili per la costruzione del coefficiente di Fresnel 
#    globale tramite somma dei singoli contributi dai diversi raggi), soglia (valore del limite 
#    di ampiezza al di sotto del quale abbandonare un raggio).

import numpy as np
import array as arr

class ellissometro:
		
	def __init__(self, precisione):
		self.precisione = precisione

	def grandell(self):

		#lettura dei valori del campione
		#read = np.loadtxt("campione.txt", delimiter='  ', unpack=True,  dtype=np.complex128)
		f = open("campione.txt", 'r')
		cont = f.readlines()
		print(cont)
		#nc = read[0]


		#lettura dei valori della sorgente
		#nc, spessori = np.loadtxt("sorgente.txt", delimiter='  ', unpack=True,  dtype=np.complex128)
		

E = ellissometro(1)
E.grandell()




'''
global strati wsuc theta nc spessori tau_sigma rho_sigma tau_pi rho_pi somma_pi somma_sigma Energia psi_0 delta_0 

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)
omega=2*pi*Energia/h;
wsuc=omega/c;

theta=theta0;
for jj=1:strati+1
    [theta(jj+1,1),tau_sigma(jj,1),rho_sigma(jj,1),tau_pi(jj,1),rho_pi(jj,1)]=interfaccia(theta(jj),nc(jj),nc(jj+1)); 
end

somma_pi=0;
somma_sigma=0;
ii=1;jj=2;x_sigma=1;x_pi=x_sigma*tan(psi_0)*exp(i*delta_0);
while ~isempty(ii)
    [ii,jj,x_pi,x_sigma]=propagazione(ii,jj,x_pi,x_sigma);
end

r_pi=somma_pi/(tan(psi_0)*exp(i*delta_0));
r_sigma=somma_sigma;
psi1=atan(abs(r_pi./r_sigma));
delta1=-angle(r_pi./r_sigma);
# psi1=atan(abs(r_sigma./r_pi));
# delta1=angle(r_sigma./r_pi);

delta1=2*pi*(delta1<0)+delta1;
'''
