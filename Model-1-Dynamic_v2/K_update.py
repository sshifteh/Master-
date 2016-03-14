
from dolfin import * 

def update_K(ind_func, K):

	assert len(K.vector()) == len(ind_func.vector())
	#K_Func.vector().axpy(-1, ind_func.vector()) #axpy work in paralel
	K.vector()[:] -= ind_func.vector()[:]

	return K

