function [psi_rad,delta_rad,R_pi,R_sigma]=ell(theta0)

% [psi,delta,R_pi,R_sigma]=ell(theta0)
%
%    Usando grandell.m calcola per diversi angoli di incidenza le 
%    grandezze ellissometriche psi, delta, e le riflettanze R_pi, R_sigma per la 
%    luce riflessa dal sistema definito in base agli argomenti forniti in ingresso.
%
%    E =         energia dei fotoni incidenti
%    theta0 =    angoli di incidenza (rispetto alla normale)
%    psi, delta: definiscono la polarizzazione iniziale (in termini di variazione rispetto alla 
%                  polarizzazione lineare orientata a pi/4 rispetto al piano di incidenza).
%    precisione: definisce il limite per un "cut-off" sul fattore di ampiezza in base al quale 
%                  viene terminata la propagazione di ciascun singolo raggio.

n=length(theta0); 
for jj=1:n
    [psi_rad(jj),delta_rad(jj),r_pi(jj),r_sigma(jj)]=grandell(theta0(jj));
end

