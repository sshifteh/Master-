from dolfin import * 
from numpy import where 

def capping(K_Func):
	"""
	Creating a 'homogenous' fluid domain, 
	in the sense that the fluid domain is level set 0,
	while solid domain is level set 1.
	"""	
	
	K = K_Func.vector() # vector of length 5000

	
	solid = where(K >= 0.5)[0]
	fluid = where(K < 0.5)[0]
	

	capping = Function(K_Func.function_space())
	capping.vector()[solid] = 1.0
	capping.vector()[fluid] = 0			 

	
	K_ = interpolate(capping, K_Func.function_space())	
	plot(K_)
	
	return K_ 


	
	
