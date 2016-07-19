from dolfin import * 
from analytical_velocityProfile import g 
from stokes_solver import stokes_solver
from main_iterative import w, mesh, Vspace, Pspace, Wspace, K  



def velocity_difference(analytical_function, experimental_function): 
	"""
	calculate errornorm for the analytical function and the experimental 		function 
	"""
	
	u_error = errornorm(analytical_function, experimental_function, degree_rise = 5)
	print u_error

if __name__ == "__main__":

	U,P = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K=K)
	diff = velocity_difference(g,U)
