function parset(E,psi,delta,prec)

% parset(Energia,psi_0,delta_0,precisione)
%
%  Definisce i parametri dello strumento usato per la misura
%  ellissometrica.
%
%    Energia =       energia dei fotoni incidenti
%
%    psi_0, delta_0: definiscono la polarizzazione iniziale (in termini di 
%                    variazione rispetto alla polarizzazione lineare 
%                    orientata a pi/4 rispetto al piano di incidenza).
%
%    precisione:     definisce il limite per un "cut-off" sul fattore di 
%                    ampiezza usato da propagazione.m che determina in base
%                    a questo quando viene terminata la propagazione di 
%                    ciascun singolo raggio. In qualche modo e' legato alla 
%                    sensibilita' del rivelatore.
%
%  Variabili globali: Energia,psi_0,delta_0,precisione

global Energia psi_0 delta_0 precisione

Energia=E;
psi_0=psi;
delta_0=delta;
precisione=prec;