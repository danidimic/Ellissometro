'''classe parset(E, psi, delta)

Energia =			energia dei fotoni incidenti

psi_0, delta_0: 	definiscono la polarizzazione iniziale (in termini di 
					variazione rispetto alla polarizzazione lineare 
					orientata a pi/4 rispetto al piano di incidenza).

precisione:			definisce il limite per un "cut-off" sul fattore di 
					ampiezza usato da propagazione.m che determina in base
					a questo quando viene terminata la propagazione di 
					ciascun singolo raggio. In qualche modo e' legato alla 
					sensibilita' del rivelatore.'''

class parset:

	def __init__(self, E, psi, delta, precisione):
		
		self.energia = E
		self.psi_0 = psi
		self.delta_0 = delta
		self.precisione = precisione
