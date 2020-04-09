function [delta_rad,psi_rad,R_pi,R_sigma]=ellplot(theta0)

% [delta,psi,R_pi,R_sigma]=ellplot(theta0)
%
%    Usando grandell.m calcola per diversi angoli di incidenza (e rappresenta graficamente) le 
%    grandezze ellissometriche psi, delta (in gradi), e le riflettanze R_pi, R_sigma per la 
%    luce riflessa dal sistema definito tramite ellmod.m.
%
%    theta0 =    angoli di incidenza (rispetto alla normale)

n=length(theta0); 
for jj=1:n
    [psi_rad(jj),delta_rad(jj),r_pi(jj),r_sigma(jj)]=grandell(theta0(jj));
end

delta_deg=180*(delta_rad)/pi;
psi_deg=180*(psi_rad)/pi;
R_pi=abs(r_pi).^2;
R_sigma=abs(r_sigma).^2;
figure(1);plot(theta0*180/pi,psi_deg,theta0*180/pi,delta_deg)
xlabel('Angle | °'),ylabel('\Psi,\Delta | °'),title('ellipsometric data'),legend({'\Psi','\Delta'})
figure(2);plot(theta0*180/pi,R_pi,theta0*180/pi,R_sigma)
xlabel('Angle | °'),ylabel('R\pi,R\sigma | #'),title('reflectivity data'),legend({'R\pi','R\sigma'})
