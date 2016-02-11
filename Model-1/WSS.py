from dolfin import * 

def WSS(U, DG1):
	"""
	Calculates the WSS over the entire domain.
	Takes as input the velocity vector U and the 
	DG1 space ShearStressSpace, and returns the WSS
	over the entire domain.

	""" 
	

	WSS_product = 0.5*(U[0].dx(1) + U[1].dx(0))  # object type --> algebra.Product 
	WSS = project(WSS_product, DG1) # In DG1, object type --> Function 
	#plot(WSS, interactive = True, title = 'WSS') # plots the WSS over the entire 	        
	
	return WSS 	
