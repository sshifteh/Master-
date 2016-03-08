from dolfin import * 

def alpha(u, K):
	C = Constant(1e10) # high value outside the vessel
	alpha= C*K*u
	return alpha

