
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
	
	k   = interpolate(func, Pspace) # control                     # Functions object 
	plot(k, interactive = True, title = 'tube with CG elements')

	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, k), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)

	# boundary conditions : 
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	 	
	def p_boundaryRightWall(x, on_bnd):
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd # why is it read on the entire unitsquare boundary?S 
	
	bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)), u_boundaryBottom)
	
	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])
	bc_u_B  = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	bc_p_R  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	bcs = [bc_u_B, bc_p_R]


	# Solve system of linear equations
	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
		
	U, P = UP.split()
	plot(U, interactive = True, title= "Velocity, vector valued function, approximated by FE")
	plot(P, interactive = True, title= "Pressure, scalar valued function, appriximated by FE")

	stress = project(0.5*(grad(U) + transpose(grad(U))), Tspace) # or is it called strain? 
	#plot(stress.sub(0), title="stress", interactive=True)
	plot(stress[0,0], title="stress - du/dx", interactive=True) # gradient in the x-direction
	plot(stress[0,1], title="stress - du/dy", interactive=True)	
	plot(stress[1,0], title="stress-dv/dx", interactive=True)
	plot(stress[1,1], title="stress-dv/dy", interactive=True)

	shear_stress = project(0.5*(U[0].dx[1] + U[1].dx[0]),Tspace) 

	shear_stress = project((0.5*(U[0].dx[1] + U[1].dx[0])), ShearStressSpace) 
	array = shear_stress.vector.array() # get numpy vector of all values of the shrear stress function
        # do some computations with array
	indicator_function = Function(ShearStresSpace) # indicator function is a function over the DG 1 shearstressSpace
	indicator_function.vector()[:] = new_array     # and we copy all the values from the shear stress array into the indicator function.


	# d1 = something, d2 = something
	# eps = something 
	# if abs(new_array) > bla bla I should do some testes:
	#	y = eps 
	# else etc. 


 	
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
	


	





 	return U, P






def main():
	from dolfin import * 
	w   = 0.3
	
	# Expression for a Straight Tube 	
	f1  = Expression('(x[0] > w && x[0] < 1 - w ) ? 0.0 : 1.0', w=w)   
	# Expression for a T shape
	f2  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 && x[0] > w )? 0.0:1.0', w=w)	
	# Expression for a L shape 
	f3  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7  )? 0.0:1.0', w=w)	
	
	
	U, P = stokes_solver(w=w, func=f3)
	#epsilon(U)
	



if __name__ == '__main__':
	from dolfin import * 	
	main()

	
