
from dolfin import * 
#from indicator_function_expression import ind_func_expr
from K_update import K_update
from capping import capping
from stokes_solver import stokes_solver 
from WSS import WSS
from WSS_bdry import WSS_bdry  
from indicator_funciton_conditions2 import ind_func_expr


n = 50
w=0.3
mesh = UnitSquareMesh(n,n)
Pspace = FunctionSpace(mesh, 'CG',1)
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace]) 
DG1 =FunctionSpace(mesh, 'DG', 1) # ShearStressSpace
DG0 = FunctionSpace(mesh, 'DG', 0)	  # DPspace CHANGE this name to DG0 much more intuitive than discont. pressure space


K_1 = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 )? 0.0:1.0', w=w)






def iterative(K_1):

	U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_1, n=n )


	# Verification: 

	# alternative 1:
	# construct an analytical solution 
	# method of manufactored solutions 
	# first come up with a random slution
	# plut into PDE 
	# doesnt solve original homogeneous PDE exactly 
	# residual is then considered the source term 
	
		
	# alternative 2:
	# really fine mesh , solve a blind model on that mesh 
	# this is reference(considered exact) solution
	
		
		
	# actual convergence test(after alt 1 or 2):
	# mesh refinement test refining the mesh see if rate of convergence goes to 2.
	# mesh.refine() OR refine(mesh)
	# give a mesh, creates new mesh that is finer. 
	# for each mesh solve PDE
	# compute difference of u_e - u
	# in a norm 
	# for each mesh 

	# start with a very coarse one 
	# after halfin 4 times,gives a really fine mesh and solving will get slow. 

	# example:
	# for p2, half the mesh size 
	# error reduction of factor 4 

	




	# do it at the beginning


	wss = WSS(U=U, DG1=DG1)
	wss_bdry = WSS_bdry(K_Func=K_Func, DPspace=DG0, WSS=wss)
 	ind_f = ind_func_expr(WSS_bdry = wss_bdry)
	K_new = K_update(ind_func = ind_f, K_Func = K_Func)
	plot(K_new, interactive = True, title = 'K new ')
	K_capped = capping(K_new)
	plot(K_capped, interactive = False, title = 'K capped')
	stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_capped, n=n )

	return K_capped 
	

iterative(K_1=K_1)

for i in range(5):
	K_updated = iterative(K_1 = K_1)
	K_1= K_updated # this didnt work 
	


