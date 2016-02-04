from dolfin import * 
	

def indicator_function_conditions(DPspace,WSS_bdry):
	"""
	Builds the indicator function with conditions.
	For the 
	"""

	indicator_Function = Function(DPspace) # in DG0, 0 to begin with 	
	indicator_Function_values =  indicator_Function.vector().array()

	# Treshold values 
	d1 = 1e-5; d2 = 20 ; eps = 1.0 
	
	indicator_Function_values[WSS_bdry_Function.vector() < d1] = -eps
	indicator_Function_values[np.logical_and(WSS_bdry_Function.vector() > d1, WSS_bdry_Function.vector() < d2)] = 0
	indicator_Function_values[WSS_bdry_Function.vector() > d2] = eps 
	

	indicator_Function.vector().set_local(indicator_Function_values) #same as copying by [:]	
	plot(indicator_Function, interactive = True, title = 'Indicator function after having been updated')	
	
	return indicator_Function
