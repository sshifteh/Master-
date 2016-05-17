from dolfin import * 

def alpha(u, K):
	C = Constant(1E5) # high value outside the vessel
	alpha= C*K*u
	return alpha

