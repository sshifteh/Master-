
from dolfin import * 

def K_update(ind_func, K_Func):
	K_Func.vector().axpy(+1, project(ind_func, K_Func.function_space()).vector())
	#plot(K_Func, interactive = True, title = 'K update')
	

	return K_Func

