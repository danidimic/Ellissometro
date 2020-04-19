# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:49:08 2020

@author: Davide Deca
"""
#variabili inizializzate

import numpy as np
import cmath

#funzione per propagare i raggi luminosi
def propagazione(ii, jj, x_pi, x_sigma, strati, wsuc, theta, nc, spessori, tau_sigma, rho_sigma, tau_pi, rho_pi, somma_pi, somma_sigma, precisione):
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
    tau_sigma : coefficienti di trasmissione dei raggi plarizzati s.
    rho_sigma : coefficienti di riflessione dei raggi plarizzati s.
    tau_pi : coefficienti di trasmissione dei raggi plarizzati p.
    rho_pi : coefficienti di riflessione dei raggi plarizzati p.
    somma_pi : contributo parziale, propveniente dalle chiamate precedenti 
               all'ampiezza p in uscita.
    somma_sigma : contributo parziale, propveniente dalle chiamate precedenti 
                  all'ampiezza p in uscita.
    precisione : soglia al di sotto della quale non considero più i raggi.

    Returns
    -------
    ii_new : strati di provenienza dei nuovi raggi.
    jj_new : strati di arrivo dei nuovi raggi.
    y_pi : nuove ampiezze dei raggi p.
    y_sigma : nuove ampiezze dei raggi s.

    """
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
    
    somma_pi = 0 #qui metterò le ampiezze del raggio totale uscente
    somma_sigma = 0
    
    phase = 0 #solo di supporto per i conti
    
    #individuo i raggi che salgono e quelli che scendono
    for x in range (0, nraggi): #ci va un +1? 
        
        #elimino raggi troppo flebili
        if abs(x_pi[x])+abs(x_sigma[x])>precisione:
                
            #se i raggi sono nel mezzo, ho un raggio riflesso, che risolvo subito, e uno trasmesso sotto
            if ii[x] == 0 and jj[x] == 1:
                somma_pi = somma_pi + x_pi[x]*rho_pi[0]
                somma_sigma = somma_sigma + x_sigma[x]*rho_sigma[0]  
                #devo aggiungere il raggio trasmesso che viene prodotto (avrà ii=1, jj=2)
                ii_superficial_dt.append(1)
                jj_superficial_dt.append(2)
                
                phase = (-wsuc*(spessori[1])/np.cos(theta[1]))*nc[1] #CONTROLLA COME INDICIZZI GLI STRATI E THETA!!!!!!!!
                
                # -polarizzati p
                ydt_sup_pi.append(x_pi[x]*tau_pi[0]*cmath.exp(phase*1j))
                # -polarizzati s
                ydt_sup_sigma.append(x_sigma[x]*tau_sigma[0]*cmath.exp(phase*1j))
                
            #se sono nel primo strato e stanno salendo ho un raggio trasmesso, ...idem
            elif ii[x] == 1 and jj[x] == 0:
                somma_pi = somma_pi + x_pi[x]*tau_pi[x]
                somma_sigma = somma_sigma + x_sigma[x]*tau_sigma[0] #vedi formule coefficienti di Fresnel!
                #devo aggiungere il raggio che manca (come prima)
                ii_superficial_ur.append(1)
                jj_superficial_ur.append(2)
                
                phase = (-wsuc*(spessori[1])/np.cos(theta[ii[1]]))*nc[1] #CONTROLLA COME INDICIZZI GLI STRATI!!!!!!!!
                
                # -polarizzati p
                yur_sup_pi.append(x_pi[x]*(-rho_pi[0])*cmath.exp(phase*1j))
                # -polarizzati s
                yur_sup_sigma.append(x_sigma[x]*(-rho_sigma[0])*cmath.exp(phase*1j))
                
                #forse si può migliorare il codice, evitando di usare irho_pi ecc...
            else:
                if ii[x] > jj[x]:
                    indici_raggi_su.append(x)
                    irho_pi[x] = -rho_pi[jj[x]] #quando ho riflessione, ho interazione con l'interfaccia che sta sopra...
                    irho_sigma[x] = -rho_sigma[jj[x]] #...in direzione "opposta", quindi... (vedi formule coefficienti di Fresnel)
                    itau_pi[x] = tau_pi[jj[x]] #anche qui: vedi formule e te ne convinci!
                    itau_sigma[x] = tau_sigma[jj[x]]
                elif ii[x] < jj[x] and ii[x] < strati+1: #non analizzo più i raggi che vanno nel substrato
                    indici_raggi_giù.append(x)
                    irho_pi[x] = rho_pi[ii[x]]
                    irho_sigma[x] = rho_sigma[ii[x]]
                    itau_pi[x] = tau_pi[ii[x]]
                    itau_sigma[x] = tau_sigma[ii[x]]
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

        phase = (-wsuc*(spessori[ii[indici_raggi_su[x]]])/np.cos(theta[ii[indici_raggi_su[x]]]))*nc[ii[indici_raggi_su[x]]]
     
        # -polarizzati p
        yur_pi[x] = x_pi[indici_raggi_su[x]]*irho_pi[indici_raggi_su[x]]*cmath.exp(phase*1j)
        # -polarizzati s
        yur_sigma[x] = x_sigma[indici_raggi_su[x]]*irho_sigma[indici_raggi_su[x]]*cmath.exp(phase*1j)
        
        #calcolo la direzione dei nuovi raggi riflessi (verso il basso)
        ii_ur[x] = ii[indici_raggi_su[x]]
        jj_ur[x] = ii[indici_raggi_su[x]]+1
            
        #vengono prodotti raggi trasmessi

        phase = (-wsuc*(spessori[ii[indici_raggi_su[x]]-1])/np.cos(theta[ii[indici_raggi_su[x]]-1]))*nc[ii[indici_raggi_su[x]]-1]
        
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
        
        phase = (-wsuc*(spessori[ii[indici_raggi_giù[x]]])/np.cos(theta[ii[indici_raggi_giù[x]]]))*nc[ii[indici_raggi_giù[x]]]
     
        # -polarizzati p
        ydr_pi[x] = x_pi[indici_raggi_giù[x]]*irho_pi[indici_raggi_giù[x]]*cmath.exp(phase*1j)
        # -polarizzati s
        ydr_sigma[x] = x_sigma[indici_raggi_giù[x]]*irho_sigma[indici_raggi_giù[x]]*cmath.exp(phase*1j)
        
        #calcolo la direzione dei nuovi raggi riflessi (verso l'alto)
        ii_dr[x] = ii[indici_raggi_giù[x]]
        jj_dr[x] = ii[indici_raggi_giù[x]]-1
            
        #vengono prodotti raggi trasmessi

        phase = (-wsuc*(spessori[ii[indici_raggi_giù[x]]+1])/np.cos(theta[ii[indici_raggi_giù[x]]+1]))*nc[ii[indici_raggi_giù[x]]+1]
        
        # -polarizzati p
        ydt_pi[x] = x_pi[indici_raggi_giù[x]]*itau_pi[indici_raggi_giù[x]]*cmath.exp(phase*1j)
        # -polarizzati s
        ydt_sigma[x] = x_sigma[indici_raggi_giù[x]]*itau_sigma[indici_raggi_giù[x]]*cmath.exp(phase*1j)
        
        #calcolo la direzione dei nuovi raggi trasmessi (verso il basso)
        ii_dt[x] = ii[indici_raggi_giù[x]]+1
        jj_dt[x] = ii[indici_raggi_giù[x]]+2
        
    #nuovi array (liste) con le direzioni dei nuovi raggi prodotti    
    ii_new = ii_ur + ii_ut + ii_dr + ii_dt + ii_superficial_dt + ii_superficial_ur
    jj_new = jj_ur + jj_ut + jj_dr + jj_dt + jj_superficial_dt + jj_superficial_ur
    
    #ampiezze dei nuovi raggi prodotti
    y_pi = yur_pi + yut_pi + ydr_pi + ydt_pi + ydt_sup_pi + yur_sup_pi
    y_sigma = yur_sigma + yut_sigma + ydr_sigma + ydt_sigma + ydt_sup_sigma + yur_sup_sigma
    
    #restituisci anche i contributi diretti alle ampiezze in somma_
    
    return ii_new, jj_new, y_pi, y_sigma, somma_pi, somma_sigma
