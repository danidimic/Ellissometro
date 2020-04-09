function d=elld(p,maschera,psi_in,delta_in,theta0,varargin)

%  d=elld(p,maschera,psi_in,delta_in,theta0)
%
%    Usando ell.m calcola per diversi angoli di incidenza le grandezze 
%    ellissometriche psi, delta,....
%
%    theta0 =    angoli di incidenza (rispetto alla normale)

global Energia psi_0 delta_0 precisione iterazione

iterazione=iterazione+1;

pp=maschera;
pp(isnan(maschera))=p;

pp=abs(pp);
if length(pp)==4
    ellmod(pp(1)+i*pp(2),pp(3)+i*pp(4))
elseif length(pp)==7
    ellmod(pp(1)+i*pp(2),pp(3)+i*pp(4),pp(5),pp(6)+i*pp(7))
elseif length(pp)==10
    ellmod(pp(1)+i*pp(2),pp(3)+i*pp(4),pp(5),pp(6)+i*pp(7),pp(8),pp(9)+i*pp(10))
elseif length(pp)==13
    ellmod(pp(1)+i*pp(2),pp(3)+i*pp(4),pp(5),pp(6)+i*pp(7),pp(8),pp(9)+i*pp(10),pp(11),pp(12)+i*pp(13))
elseif length(pp)==16
    ellmod(pp(1)+i*pp(2),pp(3)+i*pp(4),pp(5),pp(6)+i*pp(7),pp(8),pp(9)+i*pp(10),pp(11),pp(12)+i*pp(13),pp(14),pp(15)+i*pp(16))
end


psi_0=pp(length(pp)-1);
delta_0=pp(length(pp));

nn=length(theta0);
theta_line=[min(theta0):(max(theta0)-min(theta0))/(5*nn-1):max(theta0)];
[psi_sim,delta_sim]=ell(theta0);
parset(Energia,psi_0,delta_0,20*precisione)
[psi_line,delta_line]=ell(theta_line);
parset(Energia,psi_0,delta_0,precisione/20)
if nargin>5
    figure(varargin{1})
else
    figure(3)   
end
plot(theta0*180/pi,psi_in*180/pi,'x',theta_line*180/pi,psi_line*180/pi,theta0*180/pi,psi_sim*180/pi,'o',theta0*180/pi,delta_in*180/pi,'x',theta_line*180/pi,delta_line*180/pi,theta0*180/pi,delta_sim*180/pi,'o')
xlabel('Angle | °'),ylabel('\Psi,\Delta | °'),title(['fitting ellipsometric data (' num2str(iterazione) ')']),legend({'\Psi_{exp}','draft \Psi','\Psi_{sim}','\Delta_{exp}','draft \Delta','\Delta_{sim}'})
drawnow
d=norm(min(abs([psi_sim,delta_sim]-[psi_in,delta_in]),abs([psi_sim,delta_sim]-[psi_in,delta_in]-2*pi)));