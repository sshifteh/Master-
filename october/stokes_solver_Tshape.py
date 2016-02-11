
def alpha(u, k):
	C = Constant(1e5) # high value outside the vessel
	return C*k*u



def stokes_solver(w, func, n = 50 ):
	mesh = UnitSquareMesh(n,n)
	
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 

	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
	
	k   = interpolate(func, Pspace)                                               # Functions object 
	plot(k, interactive = True, title = 'tube with CG elements')

	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, k), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)


	# boundary conditions: 
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	 
	def p_boundaryRightWall(x, on_bnd):
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd # why is it read on the entire unitsquare boundary?S 
	
	def p_boundaryLeftWall(x, on_bnd):
		return x[0] < DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd
	
	bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)), u_boundaryBottom)
	bc_p_R  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	bc_p_L  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryLeftWall)
	bcs = [bc_u_B, bc_p_R, bc_p_L]		
	

	# Solve system of linear equations
	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
		
	U, P = UP.split()
	plot(U, interactive = True, title= "Velocity, vector valued function, approximated by FE")
	plot(P, interactive = True, title= "Pressure, scalar valued function, appriximated by FE")

	return U,P


if __name__ == '__main__':
	from dolfin import * 
	w   = 0.3
	
	# Expression for a Straight Tube 	
	f1  = Expression('(x[0] > w && x[0] < 1 - w ) ? 0.0 : 1.0', w=w)   
	# Expression for a T shape
	f2  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 && x[0] > w )? 0.0:1.0', w=w)	
	# Expression for a L shape 
	f3  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7  )? 0.0:1.0', w=w)	
	
	stokes_solver(w=w, func=f2)
	


