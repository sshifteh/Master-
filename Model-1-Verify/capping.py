from dolfin import * 
from numpy import where 

def capping(K):
	"""
	The isocontour is at 1/2
	"""	
	
	Kvec = K.vector() # vector of length 5000

	solid = where(Kvec >= 0.5)[0]
	fluid = where(Kvec < 0.5)[0]

	Kvec[solid] = 1.0
	Kvec[fluid] = 0			 
