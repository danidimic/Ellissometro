function [iii,jjj,y_pi,y_sigma]=propagazione(ii,jj,x_pi,x_sigma)

% [iii,jjj,y_pi,y_sigma]=propagazione(ii,jj,x_pi,x_sigma)
%
%    Simula la propagazione di raggi luminosi attraverso una serie di interfacce.
%    I vettori ii e jj indicano mezzo di provenienza e mezzo oltre l'interfaccia
%    dove i mezzi sono indicizzati da 1 a N con 1 che corrisponde al mezzo in cui 
%    e' immerso il sistema e nell'ordine i mezzi che il fascio incontra al 
%    passaggio della successione di interfacce. Ad ogni interfaccia il fascio viene
%    diviso in componente trasmessa e componente riflessa. I vettori x_pi e x_sigma 
%    contengono i termini di ampiezza (complessa) che i diversi raggi (nelle due
%    componenti di polarizzazione) si portano dietro ed a cui, ad ogni passaggio 
%    attraverso un'interfaccia, si aggiunge un fattore secondo le relazioni di Fresnel.
%    Ad ogni passaggio vengono anche introdotti i fattori che rendono conto della 
%    propagazione attraverso gli strati relativi (a seconda che si segua il fascio
%    trasmesso o quello riflesso). La propagazione dei raggi termina quando i) il fascio 
%    emerge dalla prima interfaccia. In questo caso il fattore di ampiezza relativo
%    viene sommato nel coefficiente di Fresnel globale. ii) Il fattore di attenuazione 
%    (per entrambe le componenti pi e sigma) scende in modulo sotto un valore limite. 
%    iii) il raggio emerge dall'ultima interfaccia verso il substrato.
%
%    Variabili globali: strati (numero di strati=numero di interfacce-1), wsuc (pulsazione su c),
%    theta (vettore -lungo (strati+2)- che contiene, per ogni mezzo, gli angoli tra il raggio e 
%    la normale alle interfacce), nc (indici di rifrazione dei mezzi), spessori (spessori degli 
%    strati; vale zero -irrilevante- per primo mezzo e substrato), tau_x/rho_x (trasmissivita' 
%    e riflettivita' di ciascuna interfaccia -vettore lungo (strati+1)), somma_x (variabili per 
%    la costruzione del coefficiente di Fresnel globale tramite somma dei singoli contributi dai 
%    diversi raggi), soglia (valore del limite di cui al punto ii).

global strati wsuc theta nc spessori tau_sigma rho_sigma tau_pi rho_pi somma_pi somma_sigma precisione 

spessori=abs(spessori);

ydr_pi=[];ydt_pi=[];yur_pi=[];yut_pi=[];ydr_sigma=[];ydt_sigma=[];yur_sigma=[];yut_sigma=[];
iiidr=[];iiidt=[];iiiur=[];iiiut=[];jjjdr=[];jjjdt=[];jjjur=[];jjjut=[];

% determina le riflettivita'/trasmissivita' da applicare in funzione del verso in cui e' attraversata l'interfaccia 
ijd=find(ii<jj);
iju=find(ii>jj);
irho_sigma(ijd)=rho_sigma(ii(ijd));
irho_sigma(iju)=-rho_sigma(jj(iju));
irho_pi(ijd)=rho_pi(ii(ijd));
irho_pi(iju)=-rho_pi(jj(iju));
itau_sigma(ijd)=tau_sigma(ii(ijd));
itau_sigma(iju)=tau_sigma(jj(iju));
itau_pi(ijd)=tau_pi(ii(ijd));
itau_pi(iju)=tau_pi(jj(iju));
irho_pi=irho_pi';
irho_sigma=irho_sigma';
itau_pi=itau_pi';
itau_sigma=itau_sigma';

% raggi che scendono
ida=intersect(find((abs(x_pi)+abs(x_sigma))>precisione),ijd);
idb=intersect(find(ii<strati+1),ida);
idn=intersect(find(ii>1),ida);
id1=intersect(find(ii==1),ijd);
% -riflessi
ydr_pi=x_pi(idn).*irho_pi(idn).*exp(-i*wsuc*(spessori(ii(idn))./cos((theta(ii(idn))))).*nc(ii(idn)));
%ydr_pi=x_pi(idn).*irho_pi(idn).*exp(-i*wsuc*(spessori(ii(idn))./cos(real(theta(ii(idn))))).*nc(ii(idn)));
%ydr_sigma=x_sigma(idn).*irho_sigma(idn).*exp(-i*wsuc*(spessori(ii(idn))./cos(real(theta(ii(idn))))).*nc(ii(idn)));
ydr_sigma=x_sigma(idn).*irho_sigma(idn).*exp(-i*wsuc*(spessori(ii(idn))./cos((theta(ii(idn))))).*nc(ii(idn)));
% --prima interfaccia
somma_pi=somma_pi+sum(x_pi(id1).*irho_pi(id1));
somma_sigma=somma_sigma+sum(x_sigma(id1).*irho_sigma(id1));
iiidr=ii(idn);
jjjdr=ii(idn)-1;
% -trasmessi
ydt_pi=x_pi(idb).*itau_pi(idb).*exp(-i*wsuc*(spessori(ii(idb)+1)./cos((theta(ii(idb)+1)))).*nc(ii(idb)+1));
%ydt_pi=x_pi(idb).*itau_pi(idb).*exp(-i*wsuc*(spessori(ii(idb)+1)./cos(real(theta(ii(idb)+1)))).*nc(ii(idb)+1));
%ydt_sigma=x_sigma(idb).*itau_sigma(idb).*exp(-i*wsuc*(spessori(ii(idb)+1)./cos(real(theta(ii(idb)+1)))).*nc(ii(idb)+1));
ydt_sigma=x_sigma(idb).*itau_sigma(idb).*exp(-i*wsuc*(spessori(ii(idb)+1)./cos((theta(ii(idb)+1)))).*nc(ii(idb)+1));
iiidt=ii(idb)+1;
jjjdt=ii(idb)+2;

% raggi che salgono
iua=intersect(find((abs(x_pi)+abs(x_sigma))>precisione),iju);
iun=intersect(find(ii>2),iua);
iu2=intersect(find(ii==2),iju);
% -riflessi
yur_pi=x_pi(iua).*irho_pi(iua).*exp(-i*wsuc*(spessori(ii(iua))./cos((theta(ii(iua))))).*nc(ii(iua)));
%yur_pi=x_pi(iua).*irho_pi(iua).*exp(-i*wsuc*(spessori(ii(iua))./cos(real(theta(ii(iua))))).*nc(ii(iua)));
%yur_sigma=x_sigma(iua).*irho_sigma(iua).*exp(-i*wsuc*(spessori(ii(iua))./cos(real(theta(ii(iua))))).*nc(ii(iua)));
yur_sigma=x_sigma(iua).*irho_sigma(iua).*exp(-i*wsuc*(spessori(ii(iua))./cos((theta(ii(iua))))).*nc(ii(iua)));
iiiur=ii(iua);
jjjur=ii(iua)+1;
% -trasmessi
yut_pi=x_pi(iun).*itau_pi(iun).*exp(-i*wsuc*(spessori(ii(iun)-1)./cos((theta(ii(iun)-1)))).*nc(ii(iun)-1));
%yut_pi=x_pi(iun).*itau_pi(iun).*exp(-i*wsuc*(spessori(ii(iun)-1)./cos(real(theta(ii(iun)-1)))).*nc(ii(iun)-1));
%yut_sigma=x_sigma(iun).*itau_sigma(iun).*exp(-i*wsuc*(spessori(ii(iun)-1)./cos(real(theta(ii(iun)-1)))).*nc(ii(iun)-1));
yut_sigma=x_sigma(iun).*itau_sigma(iun).*exp(-i*wsuc*(spessori(ii(iun)-1)./cos((theta(ii(iun)-1)))).*nc(ii(iun)-1));
% --prima interfaccia
somma_pi=somma_pi+sum(x_pi(iu2).*itau_pi(iu2));
somma_sigma=somma_sigma+sum(x_sigma(iu2).*itau_sigma(iu2));
iiiut=ii(iun)-1;
jjjut=ii(iun)-2;


%fine
y_pi=[ydr_pi;ydt_pi;yur_pi;yut_pi];
y_sigma=[ydr_sigma;ydt_sigma;yur_sigma;yut_sigma];
iii=[iiidr;iiidt;iiiur;iiiut];
jjj=[jjjdr;jjjdt;jjjur;jjjut];
