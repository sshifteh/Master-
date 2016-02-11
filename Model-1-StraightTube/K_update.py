
from dolfin import * 

def K_update(ind_func, K_Func):

	assert len(K_Func.vector()) == len(ind_func.vector())
	#K_Func.vector().axpy(-1, ind_func.vector())


	K_Func.vector()[:] -= ind_func.vector()[:]
		


	return K_Func

