function ellmod(nc0,varargin)

% ellmod(nc0,nc1,[h1,..,[ncx,hx],..,ncN])
%
%    Definisce il modello di un sistema su cui effettuare la simulazione di misure ellissometriche
%    Il modello e' definito in base agli argomenti forniti in ingresso.
%    Indici di rifrazione complessi nc=n-i*k
%
%    nc0 =       indice di rifrazione del mezzo dal quale incide il raggio
%    nc1 =       indice di rifrazione del primo mezzo sul quale incide il raggio
%    h1 =        spessore dello strato corrispondente (da omettere se si tratta dell'ultimo mezzo,
%                  cioe' del substrato)
%    [ncx,hx]:   indice di rifrazione e spessore di ciascuno strato che compone il sistema
%    ncN =       indice di rifrazione dell'ultimo mezzo, cioe' il subsrato (e' gia' nc1 nel caso
%                  di un singolo mezzo).

global strati nc spessori

spessori=0;
nc=nc0;
strati=(nargin-2)/2;
for jj=1:strati
    nc(jj+1,1)=varargin{2*jj-1};
    spessori(jj+1,1)=varargin{2*jj};
end
nc(strati+2,1)=varargin{nargin-1};
ii=find(real(nc)<0);
nc(ii)=-nc(ii);
ii=find(imag(nc)>0);
nc(ii)=conj(nc(ii));