from dolfin import * 
from numpy import where 

def capping(K):
	"""
	Creating a 'homogenous' fluid domain, 
	in the sense that the fluid domain is level set 0,
	while solid domain is level set 1.
	The isocontour is at 1/2
	"""	
	
	Kvec = K.vector() # vector of length 5000

	solid = where(Kvec >= 0.5)[0]
	fluid = where(Kvec < 0.5)[0]

	Kvec[solid] = 1.0
	Kvec[fluid] = 0			 
