import math
import cmath
import numpy as np

from campione import campione 
import sorgente

c=299792458; #(m/s);
h=4.13566743e-15; #(eV s)

class interazione:

	def __init__(self, precisione, campione, sorgente):
		
		self.Theta = []		#vettore contenente tutti gli angoli di incidenza
		self.rho_S = [] 	#vettore contenente tutte le riflessività ortogonali
		self.rho_P = [] 	#vettore contenente tutte le riflessività parallele
		self.tau_S = [] 	#vettore contenente tutte le trasmissività ortogonali
		self.tau_P = [] 	#vettore contenente tutte le trasmissività parallele
	
		self.campione = campione	#campione da analizzare
		self.sorgente = sorgente	#sorgente di radiazione
		self.precisione = precisione

		self.somma_pi = 0
		self.somma_sigma = 0
		self.ii = [0]
		self.jj = [1]
		self.y_pi = [0]
		self.y_sigma = [0]



	#funzione per angoli rifratti e coefficienti rho, tau per ogni interfaccia 
	def interfaccia(self, theta0, nc0, nc1):
		'''
		Funzione interfaccia 
		
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
		per la luce polarizzata rispettivamente parallelamente e ortogonalmente al piano di incidenza 

		'''
		theta1 = math.asin( nc0.real*math.sin(theta0)/nc1.real )  #angolo di rifrazione

		r = nc0/nc1
		a = nc1/nc0*math.cos(theta1)/math.cos(theta0)
		b = (r**2)*a

		rho_sigma = (1-a)/(1+a);  #riflessività ortogonale
		rho_pi = (1-b)/(1+b);	  #riflessività parallela
		
		tau_sigma = 2/(1+a)		  #trasmissività ortogonale
		tau_pi = 2*r/(1+b);		  #trasmissività parallela
	
		return [theta1, tau_sigma, rho_sigma, tau_pi, rho_pi]



	#funzione per calcolare tutti i theta, rho, tau necessari per 'propagazione'
	def inizializza(self, theta0):	#TODO i parametri nc e strati devono essere passati dal campione
		"""
		Funzione che determina tutti gli angoli di rifrazione theta, 
		le riflessvità e trsmissività rho_X e tau_x necessari 
		per la funzione propagazione
		"""
		self.Theta.append(theta0)
		for i in range(self.campione.strati-1):
			res = self.interfaccia(theta0, self.campione.nc[i], self.campione.nc[i+1])
			self.Theta.append(res[0])	#angolo di rifrazione
			self.tau_S.append(res[1]) 	#trasmissività ortogonale
			self.rho_S.append(res[2])	#riflessività ortogonale
			self.tau_P.append(res[3])	#trasmissività parallela
			self.rho_P.append(res[4])	#riflessività parallela



	#funzione per propagare i raggi luminosi
	def propagazione():	#TODO parametri devono essere passati da campione o sorgente
		"""
		Parameters (molti saranno nella classe "campione")
		----------
		ii : array con gli indici dei materiali di provenienza.
		jj : array con gli indici dei materiali oltre le interfacce.
		x_pi : ampiezze (complessa) dei raggi polarizzati p.
		x_sigma : ampiezze (complessa) dei raggi polarizzati s.
		strati : numero di strati che compongono il materiale.
		wsuc : frequenza angolare del ragggio/c.
		theta : angoli di incidenza della luce tra i vari strati, calcolati tramite 
		        "interfaccia".
		nc : indice di rifrazione (complesso (?)) *c (il substrato ne ha uno non nullo (?)).
		spessori : array con gli spessori dei vari strati.
		somma_pi : contributo parziale, propveniente dalle chiamate precedenti 
		           all'ampiezza p in uscita.
		somma_sigma : contributo parziale, propveniente dalle chiamate precedenti 
		              all'ampiezza p in uscita.

		Returns
		-------
		ii_new : strati di provenienza dei nuovi raggi.
		jj_new : strati di arrivo dei nuovi raggi.
		y_pi : nuove ampiezze dei raggi p.
		y_sigma : nuove ampiezze dei raggi s.

		"""
		nc = self.campione.nc
		strati = self.campione.strati
		spessori = self.campione.spessori

		x_pi = self.y_pi
		x_sigma = self.y_sigma
		
		ii = self.ii
		jj = self.jj

		omega=2*math.pi*sorgente.energia/h;
		wsuc=omega/c;

		indici_raggi_su = [] #lista degli indici
		indici_raggi_giù = []
		
		nraggi = len(ii) #numero di raggi
		
		#INIZIALIZZO a 0 itau_pi ecc...
		irho_pi = [0 for x in range(nraggi)]
		irho_sigma = [0 for x in range(nraggi)]
		itau_pi = [0 for x in range(nraggi)]
		itau_sigma = [0 for x in range(nraggi)]
		
		ii_superficial_dt = []
		jj_superficial_dt = []
		ii_superficial_ur = []
		jj_superficial_ur = []
		
		yur_sup_pi = []
		yur_sup_sigma = []
		ydt_sup_pi = []
		ydt_sup_sigma = []
		
		#self.somma_pi = 0
		#self.somma_sigma = 0
		
		phase = 0 #solo di supporto per i conti
		
		#individuo i raggi che salgono e quelli che scendono
		for x in range (0, nraggi): #ci va un +1? 
		    
		    #elimino raggi troppo flebili
		    if abs(x_pi[x])+abs(x_sigma[x])>self.precisione:
		            
		        #se i raggi sono nel mezzo, ho un raggio riflesso, che risolvo subito, e uno trasmesso sotto
		        if ii[x] == 0 and jj[x] == 1:
		            self.somma_pi += x_pi[x]*self.rho_P[0]
		            self.somma_sigma += x_sigma[x]*self.rho_S[0]  
		            #devo aggiungere il raggio trasmesso che viene prodotto (avrà ii=1, jj=2)
		            ii_superficial_dt.append(1)
		            jj_superficial_dt.append(2)
		            
		            phase = (-wsuc*(spessori[1])/np.cos(self.Theta[1]))*nc[1] #CONTROLLA COME INDICIZZI GLI STRATI E THETA!!!!!!!!
		            
		            # -polarizzati p
		            ydt_sup_pi.append(x_pi[x]*self.tau_P[0]*cmath.exp(phase*1j))
		            # -polarizzati s
		            ydt_sup_sigma.append(x_sigma[x]*self.tau_S[0]*cmath.exp(phase*1j))
		            
		        #se sono nel primo strato e stanno salendo ho un raggio trasmesso, ...idem
		        elif ii[x] == 1 and jj[x] == 0:
		            self.somma_pi += x_pi[x]*self.tau_P[x]
		            self.somma_sigma += x_sigma[x]*self.tau_S[0] #vedi formule coefficienti di Fresnel!
		            #devo aggiungere il raggio che manca (come prima)
		            ii_superficial_ur.append(1)
		            jj_superficial_ur.append(2)
		            
		            phase = (-wsuc*(spessori[1])/np.cos(self.Theta[ii[1]]))*nc[1] #CONTROLLA COME INDICIZZI GLI STRATI!!!!!!!!
		            
		            # -polarizzati p
		            yur_sup_pi.append(x_pi[x]*(-self.rho_P[0])*cmath.exp(phase*1j))
		            # -polarizzati s
		            yur_sup_sigma.append(x_sigma[x]*(-self.rho_S[0])*cmath.exp(phase*1j))
		            
		            #forse si può migliorare il codice, evitando di usare irho_pi ecc...
		        else:
		            if ii[x] > jj[x]:
		                indici_raggi_su.append(x)
		                irho_pi[x] = -self.rho_P[jj[x]] #quando ho riflessione, ho interazione con l'interfaccia che sta sopra...
		                irho_sigma[x] = -self.rho_S[jj[x]] #...in direzione "opposta", quindi... (vedi formule coefficienti di Fresnel)
		                itau_pi[x] = self.tau_P[jj[x]] #anche qui: vedi formule e te ne convinci!
		                itau_sigma[x] = self.tau_S[jj[x]]
		            elif ii[x] < jj[x] and ii[x] < strati+1: #non analizzo più i raggi che vanno nel substrato
		                indici_raggi_giù.append(x)
		                irho_pi[x] = self.rho_P[ii[x]]
		                irho_sigma[x] = self.rho_S[ii[x]]
		                itau_pi[x] = self.tau_P[ii[x]]
		                itau_sigma[x] = self.tau_S[ii[x]]
		#else elimino le componenti dei raggi troppo flebili da x_pi e x_sigma???
		    
		nraggi_up = len(indici_raggi_su)
		yur_pi = [0 for x in range(nraggi_up)]
		yut_pi = [0 for x in range(nraggi_up)]
		yur_sigma = [0 for x in range(nraggi_up)]
		yut_sigma = [0 for x in range(nraggi_up)]
		    
		ii_ur = [0 for x in range(nraggi_up)]
		jj_ur = [0 for x in range(nraggi_up)]
		ii_ut = [0 for x in range(nraggi_up)]
		jj_ut = [0 for x in range(nraggi_up)]
		    
		nraggi_down = len(indici_raggi_giù)
		ydr_pi = [0 for x in range(nraggi_down)]
		ydt_pi = [0 for x in range(nraggi_down)]
		ydr_sigma = [0 for x in range(nraggi_down)]
		ydt_sigma = [0 for x in range(nraggi_down)]
		    
		ii_dr = [0 for x in range(nraggi_down)]
		jj_dr = [0 for x in range(nraggi_down)]
		ii_dt = [0 for x in range(nraggi_down)]
		jj_dt = [0 for x in range(nraggi_down)]
		    
		#print(indici_raggi_su)
		#print(indici_raggi_giù)
		
		#aggiorno ampiezze dopo la propagazione
		#raggi che salgono
		for x in range (0, nraggi_up):
		    #vengono prodotti raggi riflessi

		    phase = (-wsuc*(spessori[ii[indici_raggi_su[x]]])/np.cos(self.Theta[ii[indici_raggi_su[x]]]))*nc[ii[indici_raggi_su[x]]]
		 
		    # -polarizzati p
		    yur_pi[x] = x_pi[indici_raggi_su[x]]*irho_pi[indici_raggi_su[x]]*cmath.exp(phase*1j)
		    # -polarizzati s
		    yur_sigma[x] = x_sigma[indici_raggi_su[x]]*irho_sigma[indici_raggi_su[x]]*cmath.exp(phase*1j)
		    
		    #calcolo la direzione dei nuovi raggi riflessi (verso il basso)
		    ii_ur[x] = ii[indici_raggi_su[x]]
		    jj_ur[x] = ii[indici_raggi_su[x]]+1
		        
		    #vengono prodotti raggi trasmessi

		    phase = (-wsuc*(spessori[ii[indici_raggi_su[x]]-1])/np.cos(self.Theta[ii[indici_raggi_su[x]]-1]))*nc[ii[indici_raggi_su[x]]-1]
		    
		    # -polarizzati p
		    yut_pi[x] = x_pi[indici_raggi_su[x]]*itau_pi[indici_raggi_su[x]]*cmath.exp(phase*1j)
		    # -polarizzati s
		    yut_sigma[x] = x_sigma[indici_raggi_su[x]]*itau_sigma[indici_raggi_su[x]]*cmath.exp(phase*1j)
		    
		    #calcolo la direzione dei nuovi raggi trasmessi (verso l'alto)
		    ii_ut[x] = ii[indici_raggi_su[x]]-1
		    jj_ut[x] = ii[indici_raggi_su[x]]-2
		 
		#raggi che scendono
		for x in range(0, nraggi_down):
		    #vengono prodotti raggi riflessi
		    
		    phase = (-wsuc*(spessori[ii[indici_raggi_giù[x]]])/np.cos(self.Theta[ii[indici_raggi_giù[x]]]))*nc[ii[indici_raggi_giù[x]]]
		 
		    # -polarizzati p
		    ydr_pi[x] = x_pi[indici_raggi_giù[x]]*irho_pi[indici_raggi_giù[x]]*cmath.exp(phase*1j)
		    # -polarizzati s
		    ydr_sigma[x] = x_sigma[indici_raggi_giù[x]]*irho_sigma[indici_raggi_giù[x]]*cmath.exp(phase*1j)
		    
		    #calcolo la direzione dei nuovi raggi riflessi (verso l'alto)
		    ii_dr[x] = ii[indici_raggi_giù[x]]
		    jj_dr[x] = ii[indici_raggi_giù[x]]-1
		        
		    #vengono prodotti raggi trasmessi

		    phase = (-wsuc*(spessori[ii[indici_raggi_giù[x]]+1])/np.cos(self.Theta[ii[indici_raggi_giù[x]]+1]))*nc[ii[indici_raggi_giù[x]]+1]
		    
		    # -polarizzati p
		    ydt_pi[x] = x_pi[indici_raggi_giù[x]]*itau_pi[indici_raggi_giù[x]]*cmath.exp(phase*1j)
		    # -polarizzati s
		    ydt_sigma[x] = x_sigma[indici_raggi_giù[x]]*itau_sigma[indici_raggi_giù[x]]*cmath.exp(phase*1j)
		    
		    #calcolo la direzione dei nuovi raggi trasmessi (verso il basso)
		    ii_dt[x] = ii[indici_raggi_giù[x]]+1
		    jj_dt[x] = ii[indici_raggi_giù[x]]+2
		    
		#nuovi array (liste) con le direzioni dei nuovi raggi prodotti    
		self.ii = ii_ur + ii_ut + ii_dr + ii_dt + ii_superficial_dt + ii_superficial_ur
		self.jj = jj_ur + jj_ut + jj_dr + jj_dt + jj_superficial_dt + jj_superficial_ur
		
		#ampiezze dei nuovi raggi prodotti
		self.y_pi = yur_pi + yut_pi + ydr_pi + ydt_pi + ydt_sup_pi + yur_sup_pi
		self.y_sigma = yur_sigma + yut_sigma + ydr_sigma + ydt_sigma + ydt_sup_sigma + yur_sup_sigma


