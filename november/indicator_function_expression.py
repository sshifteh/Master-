
from dolfin import * 



def ind_func_expr(w, Pspace):
	"""
	The indicator function created with an expression
	"""
	
	#ind_f_exp = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] <= 0.7 || x[0] > 0.1 && x[0]<= 0.3 && x[1] >= 0.7 && x[1]< 0.8)? 0.0:1.0' , w=w)     


	ind_f_exp =  Expression('(x[0] > 0.1 && x[0]<= 0.3 && x[1] >= 0.7 && x[1]< 0.8)? -1.0:0.0' , w=w) # bumbp geometry      
	ind_func= interpolate(ind_f_exp, Pspace)  
	#plot(ind_func, interactive = True, title = 'ind func created with an expre')	

	return ind_func



