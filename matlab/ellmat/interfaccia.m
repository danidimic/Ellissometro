function [theta1,tau_sigma,rho_sigma,tau_pi,rho_pi]=interfaccia(theta0,nc0,nc1)

% [theta1,tau_sigma,rho_sigma,tau_pi,rho_pi]=interfaccia(theta0,nc0,nc1)
%
%    Calcola l'angolo del fascio rifratto theta1, le riflettivita' rho_x e trasmissivita' tau_x
%    corrispondenti all'angolo di incidenza theta0  per l'interfaccia tra due mezzi con indice 
%    di rifrazione nc0 e nc1. x_pi e x_sigma indicano i valori per la luce polarizzata 
%    rispettivamente parallelamente e ortogonalmente al piano di incidenza}
%
%    theta0 =    angolo di incidenza (rispetto alla normale)
%    nc0 =       indice di rifrazione del mezzo dal quale incide il raggio
%    nc1 =       indice di rifrazione del primo mezzo sul quale incide il raggio


theta1=asin(nc0*sin(theta0)/nc1);
a=nc1/nc0*cos(theta1)/cos(theta0);
r=nc0/nc1;
b=(r^2)*a;
tau_sigma=2/(1+a);
rho_sigma=(1-a)/(1+a);
tau_pi=2*r/(1+b);
rho_pi=(1-b)/(1+b);