from dolfin import * 

def alpha(u, K):
	#C= 0.1	
	#C = 1 
	#C = 10
	#C= 100
	#C= 1000
	#C = 1E4	
	C = 1E5
	alpha= C*K*u
	return alpha

