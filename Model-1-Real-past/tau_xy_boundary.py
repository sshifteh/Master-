
from dolfin import * 

def tau_xy_bdry(K_Func, DPspace, WSS):

	"""

	The implicit function, K, creates separates the R**2. The zerocontour gives and interface.
	The outer region is assigned value 1, while the inner is assigned value 0.
	This can be expressed as a function of one variable by the Heaviside function
	By definition the directional derivative of the Heaviside function is the dirac delta function.
	The volume or surface integral of a funcion over a domain is defined as the product of the function and the drac delta function.
	The dirac delta function picks out the boundary.
	
	Something similar I am trying to do here. 

	"""

	
	bdry= grad(K_Func)**2 
	tau_xy_bdry = WSS*bdry # only if ur at the bdry AND the wss is nonzero 
	

	#print max(WSS_bdry)	
	# keept the bdry and wss separate 
	# if at the bdry AND wss  larger grow the domain
	# if at bdry AND wss is less then shrink the domain 
	

	WSS_bdry_Function =  project(WSS_bdry, DPspace)

	WSS_bdry_max = max(WSS_bdry_Function.vector().array())
	print WSS_bdry_max




	#plot(WSS_bdry_Function, interactive = True, title = 'WSS_bdry_Function') 
	
	return WSS_bdry_Function 
