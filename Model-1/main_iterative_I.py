
# DOLFIN LIBRARY AND PLOTTING LIBRARY 
from dolfin import * 

from plotting import plot 
parameters['plotting_backend'] = "matplotlib"
from matplotlib import pyplot

# LEVEL SET FUNCTIONS, WSS, INDICATOR
from K_update import K_update
from capping import capping
#from stokes_solver_I import stokes_solver 
from stokes_solver_T import stokes_solver
from WSS import WSS
from WSS_bdry import WSS_bdry  
from indicator_funciton_conditions2 import ind_func_expr


# FINITE ELEMENT MESH, FUNCTION SPACES
n = 20
w=0.40 # for straigh tube 
#w= 0.01 # For Y shaped tube , crossed mesh  
mesh = UnitSquareMesh(n,n)
mesh = refine(mesh)
#plot(mesh)
Pspace = FunctionSpace(mesh, 'CG',1)
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace]) 
DG1 =FunctionSpace(mesh, 'DG', 1)  # ShearStressSpace
DG0 = FunctionSpace(mesh, 'DG', 0) # DPspace CHANGE this name to DG0. Pressure space

# I - SHAPED TUBE 
phi_I  = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)   

# T - SHAPED TUBE
#phi_T = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1- w + DOLFIN_EPS && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 )? 0.0:1.0', w=w)




def iterative(phi_I):
	U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=phi_I, n=n )
	#plot(K_Func, interactive = True, title = 'Velocity')
				
	wss = WSS(U=U, DG1=DG1)
	wss_bdry = WSS_bdry(K_Func=K_Func, DPspace=DG0, WSS=wss)
 	ind_f = ind_func_expr(WSS_bdry = wss_bdry)
	K_new = K_update(ind_func = ind_f, K_Func = K_Func)
	K_capped = capping(K_new)
	stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_capped, n=n )

	return K_capped  
	


for i in range(10):
	K_updated = iterative(phi_I = phi_I)
	plot(K_updated,interactive = True, title = 'K')
	
	plot(K_updated)	
	pyplot.suptitle("bcs2_skewed_iteration_no_%g" %(i+1))
	pyplot.savefig("exp1a_bcs2_skewed_iteration_no_%g" %(i+1))
	#pyplot.figure()	
	pyplot.show()

	phi_I = K_updated 
	


