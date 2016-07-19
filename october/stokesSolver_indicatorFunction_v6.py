
	
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
	WSS = shear_stress.vector().array() 	 	                  
      	
	#grad_k_sq = project(grad(k)**2, DPspace)
	#grad_k_array = grad_k.vector() # its not a copy 	
	
	# when multiplying with the grad_k which lives in the DG0 meaning it is a constant function 
	# and consider the indicator function which is a DG1, meaning a linear function over each element
	# consider just one element, one linear function multiplyed by one constant function, what will happen?
	# the linear function will be scaled by the constant function, fx x times c =  cx.
	# so the resulting space will be DG1 


	# use this 
	WSS_bdry = project(abs(shear_stress*grad(k)**2), DPspace).vector()
	indicator_function = Function(DPspace) 			        
	indicator_function_values =  indicator_function.vector().array() # TODO 


	d1 = 1e-5; d2 = 20 ; eps = 1.0         

	"""
	# Aternative one: for loop version                     
	start_time = time.time()                
	for i in range(len(WSS_bdry)):
		if abs(WSS_bdry[i]) < d1:  # shinking 
			#print "This is True" 
			indicator_function_values[i] = -eps

		elif  d1 < abs(WSS_bdry[i]) < d2: # no change 
			indicator_function_values[i] = 0
		
		elif abs(WSS_bdry[i]) > d2: # expanding   
			indicator_function_values[i] = eps
	
	end_time = time.time()
	time_elapsed = end_time - start_time 
	print 'time for the for loop is %.3f' %time_elapsed
	""" 

	
	# Alternative two: vectorized version : 
	start_time = time.time()	
	#WSS = np.abs(shear_stress.vector().array()) 	 	                  	
	indicator_function_values[WSS_bdry < d1] = -eps
	indicator_function_values[np.logical_and(WSS_bdry > d1, WSS_bdry < d2)] = 0
	indicator_function_values[WSS_bdry > d2] = eps 
	end_time = time.time()
	elapsed_time = end_time - start_time
	print 'time for the vectorized version is %.5f' %elapsed_time
	
	# later it should be in ufl conditionals. Dvs alle conditionals skal vaere i ufl form something.
	# beware this is research and we dont know really how it will go. But it is a master project and should go fine
	# This stage is experimentation	
	# Try writing out the math. 		

	# Put values of temp indicator function array into indicator Function object 		
	indicator_function.vector().set_local(indicator_function_values)	
	plot(indicator_function, interactive = False, title = ' indicator function') 
	
	# Update step: axpy works in parallel 
	k.vector().axpy(-1, project(indicator_function, k.function_space()).vector())
	print "::::::::::::::::::::::::", type(k) # Function.
		
	# Vectorized capping and visualizing 
	
	
	# The peaks in the figure are a result of not applying capping.
	# the indicator function expanded the solid domain which result in the same geometric description, because it s already solid.

	
	"""
	#plt.figure() 
	k_values = k.vector().array()  # make a copy to manipulate	
	#plt.plot(k_values, label = 'before')	
	k_values[k_values > 0.1] = 1
	k_values[k_values <= 0.1] = 0
	#plt.plot(k_values, label = 'after')
	#plt.show(False)
	k.vector().set_local(k_values) # interpolate the values back into the function space 
	""" 	
	plot(k, interactive = True, title = 'updated k')



	
 	return U, P

	

	
	# the shrinking and the expansion should be limited to the bnd of the domain
	# by computng the gradient of the k function 
	# the gradient of k ||grad(k)|| is zero in the solid and zero in the vessel. but nonzero in the interphace.
	# can multipy the indicator function with the gradient of k , then you restrict the update to the boundary 
	# experiment with this, and this might be the crusial part of the project i.e to find a good update scheme. 
	
	# beware shrinking and expansion is applied in the cases. So before the cases I should 
	# calculate the grad and use it ?!?  




def main():
	w   = 0.3
	# Expression for a T shape 
	f3  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7  )? 0.0:1.0', w=w) #K in math.formulation	
	U, P = stokes_solver(w=w, func=f3)
	
	

if __name__ == '__main__':
	from dolfin import * 	
	main()

	
