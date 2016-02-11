
	

	
def alpha(u, K):
	C = Constant(1e5) # high value outside the vessel
	return C*K*u


def u_boundaryBottom(x, on_bnd):
	return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	 	
def p_boundaryRightWall(x, on_bnd):
	return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd # why is it read on the entire unitsquare boundary?S 


def stokesSolver(w,n, K):
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
	
	K = interpolate(K, Pspace) # Kontroll  Using PSpace, beware with a cont func we have isoline, not with discont func
	plot(K, interactive = True, title = 'tube with CG elements')

	f = Constant([0.0,0.0]) 
	a = inner(alpha(u, k), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)

	# boundary conditions____________________________________________________________________________
	
	bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)), u_boundaryBottom)
	
	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])
	bc_u_B  = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	bc_p_R  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	bcs = [bc_u_B, bc_p_R]


	# Solve system of linear equations_______________________________________________________________
	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
		
	U, P = UP.split()	
	plot(U, interactive = True, title= "Velocity")
	plot(P, interactive = True, title= "Pressure")

	return U, P,K  




# Handling of the WSS ___________________________________________________________________________

def WSS(U,P, K):
	shear_stress = project(0.5*(U[0].dx(1) + U[1].dx(0)), ShearStressSpace) 
	plot(shear_stress, interactive = True, title = 'WSS')	                
	WSSarray = shear_stress.vector().array() #len is 15 000                                


	indicator_function = Function(ShearStressSpace) 
	indicator_func_array = indicator_function.vector().array() 
	indicator_func_array[:] = WSSarray                        
	
	# smallest WSS value 4.07399371391e-08 largest WWS value 21.7227328423
	# beware gamma is the indicator function	
	d1 = 1; d2 = 10 ; eps = 1e-5		
	
		
	for i in range(len(WSSarray)):
		if abs(WSSarray[i]) < d1:
			indicator_func_array[i] = -eps
	
		elif  d1 < abs(WSSarray[i]) < d2:
			indicator_func_array[i] = 0
			
		elif abs(WSSarray[i]) > d2:
			indicator_func_array[i] = eps
			
  
	K_array = K.vector().array()
	#initalize first the new k_updated array
	K_updated = K_array[:]
		
	for i in range(len(K_array)):
			K_updated[i] = K_array[i] - indicator_func_array[i]
		
	K_updated_projected = interpolate(K_updated, Pspace) 
	plot(K_updated_projected, interactive = True, title = 'K_updated')
	
	
	
 	
	# Update step:
	# K.vector()[:] -= project(indicator_function, K.function_space()).vector()

	# depending on funcitonspace the .vector() DOF array have different lengths.
	# P2 functions also info on nodes and edges, p1 on the nodes only.
	# extract array of Kn-1 numpyarray of all values. same array of indicator function
	# in order to make it valid, the array has to be the same lenght, 
	# therefore project first into the same space. 

	# ideally the optimization will contruct the values of the d's for us

	# do eq (8)
	# look at old geo and compare to new and se if it makes sense. 
	

	#epsilon = sym(grad(U)) # equivalent to the epsilon function
	#plot(epsilon) # Error, plotting of higher order functions are not supported
	


	

from dolfin import * 

#main______________________________________________
n= 50 
w   = 0.3
K  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7  )? 0.0:1.0', w=w) # K in math.formulation	
mesh = UnitSquareMesh(n,n)
	
Tspace = TensorFunctionSpace(mesh, "DG", 1)	
ShearStressSpace = FunctionSpace(mesh, "DG", 1)	

Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace]) 

U,P,K = stokesSolver(w,n,K)
WSS(U,P,K)




 	







	
