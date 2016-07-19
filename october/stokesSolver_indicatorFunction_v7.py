
	
show = False


def alpha(u, k):
	C = Constant(1e5) # high value outside the vessel
	return C*k*u




def stokes_solver(w, func, n = 50 ):

	# Import statement TODO these really ha ve to be here to work 
	import time 
	import numpy as np
	import matplotlib.pyplot as plt 
	
	# Mesh and  Function spaces 
	mesh = UnitSquareMesh(n,n)
	Tspace = TensorFunctionSpace(mesh, "DG", 1)	
	ShearStressSpace = FunctionSpace(mesh, "DG", 1)	
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	DPspace = FunctionSpace(mesh, 'DG', 0)	
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
	
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
		
	k   = interpolate(func, Pspace) # control                     # Function object 
	#plot(k, interactive = show, title = 'tube with CG elements')

	#plot(project(grad(k)**2, DPspace), interactive=True)

	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, k), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)







	# boundary conditions____________________________________________________________________________
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	 	
	def p_boundaryRightWall(x, on_bnd):
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd # why is it read on the entire unitsquare boundary? 
	
	bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)), u_boundaryBottom)
	
	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])
	bc_u_B  = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	bc_p_R  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	bcs = [bc_u_B, bc_p_R]








	# Solve system of linear equations_______________________________________________________________
	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
		
	U, P = UP.split()
	#plot(U, interactive = show, title= "Velocity")
	#plot(P, interactive = show, title= "Pressure")







	# WSS, Indicator function and the spacially dependent(wss) K-function  ___________________________________________________

	shear_stress = project(0.5*(U[0].dx(1) + U[1].dx(0)), ShearStressSpace)    
	plot(shear_stress, interactive = True, title = 'WSS')	                 
	WSS_bdry = project(abs(shear_stress*grad(k)**2), DPspace).vector()

	indicator_function = Function(DPspace) 			        
	indicator_function_values =  indicator_function.vector().array() # TODO 


	d1 = 1e-5; d2 = 20 ; eps = 1.0        
	
	# Alternative: vectorized version : 
	start_time = time.time()	 	 	                  	
	indicator_function_values[WSS_bdry < d1] = -eps
	indicator_function_values[np.logical_and(WSS_bdry > d1, WSS_bdry < d2)] = 0
	indicator_function_values[WSS_bdry > d2] = eps 
	end_time = time.time()
	elapsed_time = end_time - start_time
	print 'time for the vectorized version is %.5f' %elapsed_time
	

	# Put values of temp indicator function array into indicator Function object 		
	indicator_function.vector().set_local(indicator_function_values)	
	plot(indicator_function, interactive = False, title = ' indicator function') 
	
	# Update step: axpy works in parallel 
	k.vector().axpy(-1, project(indicator_function, k.function_space()).vector())
	plot(k, interactive = True, title = 'updated k')

	# I want to make this recursive now 

	
 	return U, P

	
	# the shrinking and the expansion should be limited to the bnd of the domai
	# by computng the gradient of the k function 
	# the gradient of k ||grad(k)|| is zero in the solid and zero in the vessel. but nonzero in the interphace.
	# can multipy the indicator function with the gradient of k , then you restrict the update to the boundary 
	# experiment with this, and this might be the crusial part of the project i.e to find a good update scheme. 



def main():
	w   = 0.3
	# Expression for a T shape 
	f3  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7  )? 0.0:1.0', w=w) #K in math.formulation	
	U, P = stokes_solver(w=w, func=f3)
	
	

if __name__ == '__main__':
	from dolfin import * 	
	main()

	
