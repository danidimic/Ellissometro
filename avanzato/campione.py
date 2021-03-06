#    campione(nc0,nc1,[h1,..,[ncx,hx],..,ncN])
#
#    Definisce il modello di un sistema su cui effettuare la simulazione di misure ellissometriche
#    Il modello e' definito in base agli argomenti forniti in ingresso.
#    Indici di rifrazione complessi nc=n-i*k
#
#    nc0 =       indice di rifrazione del mezzo dal quale incide il raggio
#    nc1 =       indice di rifrazione del primo mezzo sul quale incide il raggio
#    h1 =        spessore dello strato corrispondente (da omettere se si tratta dell'ultimo mezzo,
#                  cioe' del substrato)
#    [ncx,hx]:   indice di rifrazione e spessore di ciascuno strato che compone il sistema
#    ncN =       indice di rifrazione dell'ultimo mezzo, cioe' il subsrato (e' gia' nc1 nel caso
#                  di un singolo mezzo).

class campione:

    def __init__(self, nc0, varargin):

        self.spessori = [0]
        self.nc = [nc0]
        self.strati=( (len(varargin)-1) / 2 )

        self.strati=int(self.strati)

        for count in range(0,self.strati):
            self.nc.append(varargin[2*count])
            self.spessori.append(float(varargin[2*count+1].real))
        self.spessori.append(0) #MOSSA MAGICA

        self.nc.append(varargin[-1])

        for count in range(0,len(self.nc)):
            if self.nc[count].real < 0:
                self.nc[count] = - self.nc[count]

        for count in range(0,len(self.nc)):
        	if self.nc[count].imag < 0:
        		self.nc[count] = self.nc[count].conjugate()



def loadCampione():

    file = open("campione.txt","r")

    line = str(file.readline())
    nc0 = float(file.readline())
    line = str(file.readline())
    varargin = []

    while True:
        line = str(file.readline())

        if line == "Indice del substrato \n":
            break

        else:
            x = complex(line)
            varargin.append(x)

    line = str(file.readline())
    ncN = complex(line)
    varargin.append(ncN)
    
    file.close()

    camp = campione(nc0, varargin)
    return camp


