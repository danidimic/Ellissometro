% Funzioni per simulare grandezze ellissometriche
%
%
%  Definizioni del modello e del setup
%
%     ellmod - definisce il sistema su cui si effettua la misura
%     parset - imposta alcuni parametri del setup
%
%
%  Funzioni di livello basso
%
%     interfaccia - calcola angolo del raggio rifratto e coefficienti di
%                   Fresnel per l'interfaccia tra due mezzi
%     propagazione - simula la propagazione di un set di raggi in un
%                    sistema definito da ellmod
%     grandell - calcola le grandezze ellissometriche per un singolo angolo
%                di incidenza
%     elld - funzione distanza tra grandezze in ingresso e grandezze
%            simulate
%
%  Funzioni di livello alto
%
%     ell - calcola le grandezze ellissometriche per diversi angoli di
%           incidenza
%     ellplot - come ell ma con output grafico
%     ellfit - alias per fminsearch
%     
%     Autore: Paolo Piseri 05/2005