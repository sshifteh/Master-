
	
show = True


def alpha(u, k):
	C = Constant(1e5) # high value outside the vessel
	return C*k*u


def stokes_solver(w, func, n = 50 ):
	mesh = UnitSquareMesh(n,n)
	
	Tspace = TensorFunctionSpace(mesh, "DG", 1)	
	ShearStressSpace = FunctionSpace(mesh, "DG", 1)	
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
	
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
		
	k   = interpolate(func, Pspace) # control                     # Function object 
	plot(k, interactive = show, title = 'tube with CG elements')
	
	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, k), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)

	# boundary conditions____________________________________________________________________________
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	 	
	def p_boundaryRightWall(x, on_bnd):
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd # why is it read on the entire unitsquare boundary?S 
	
	bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)), u_boundaryBottom)
	
	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])
	bc_u_B  = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	bc_p_R  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	bcs = [bc_u_B, bc_p_R]


	# Solve system of linear equations_______________________________________________________________
	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
		
	U, P = UP.split()
	plot(U, interactive = show, title= "Velocity")
	plot(P, interactive = show, title= "Pressure")








	# Handling of the WSS ___________________________________________________________________________

	shear_stress = project(0.5*(U[0].dx(1) + U[1].dx(0)), ShearStressSpace)    
	plot(shear_stress, interactive = show, title = 'WSS')	                 
	#WSS = np.abs(shear_stress.vector().array()) 	 	                  
	WSS = shear_stress.vector().array() 	 	                  
      
	indicator_function = Function(ShearStressSpace) 			        
	#plot(indicator_function, interactive = show, title = ' indicator function') 

	# A test : pointer to that object
	for i in range(4):
		indicator_function_values = indicator_function.vector().array()
		print id(indicator_function_values)                 


	#print indicator_func_array is indicator_function.vector()#.array()
	# gamma is the indicator function
	
	d1 = 2; d2 = 3 ; eps = 50


                                              
	for i in range(len(WSS)):
		if abs(WSS[i]) < d1:  # shrinking 
			#print "This is True" 
			indicator_function_values[i] = -eps

		elif  d1 < abs(WSS[i]) < d2: # no change 
			indicator_function_values[i] = 0
		
		elif abs(WSS[i]) > d2: # expanding 
			indicator_function_values[i] = eps

	
		
	indicator_function.vector().set_local(indicator_function_values)	
	#indicator_func = project(indicator_func, Pspace)
	plot(indicator_function, interactive = True, title = ' indicator function') 


	# Update step:
	k.vector().axpy(-1, project(indicator_function, k.function_space()).vector())
	print "::::::::::::::::::::::::", type(k)
	#import numpy as np
	#print "max of ind-func", np.max(project(indicator_function, k.function_space()).vector().array())
	
	# remedy : interpolate k first 	
		

	plot(k, ineractive = True, title = 'updated k')




 	return U, P






def main():
	import numpy as np
	from dolfin import * 
	w   = 0.3
	# Expression for a T shape 
	f3  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7  )? 0.0:1.0', w=w) #K in math.formulation	
	U, P = stokes_solver(w=w, func=f3)
	
	

if __name__ == '__main__':
	from dolfin import * 	
	main()

	
