def alpha(u, k):
	C = Constant(0)
	return C*k*u

def first_part(w, n = 50 ):
	mesh = UnitSquareMesh(n,n)
	
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 

	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
	
	f2  = Expression('(x[0] > w && x[0] < 1 - w ) ? 0.0 : 1.0', w=w)     # CompiledExpression Object
	#f3 = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 && x[0] > w )? 0.0:1.0', w=w)	
	k   = interpolate(f2, Pspace)                                               # Functions object 
	plot(k, interactive = True, title = 'tube with CG elements')

	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, k), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)
	
	U, P = UP.split(True)

	# Making boundary conditions 
	# So I want a square with zero velocity on R and L walls of the UnitSquare 
	# zero pressure on the top 
	# and u = (0,1) at bottom, so the so the flow is pushed forward
	
	
	def u_boundaryLeftRightWall(x, on_bnd):
		return  x[0] < DOLFIN_EPS or x[0] > 1- DOLFIN_EPS and on_bnd
		 		
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and on_bnd
	 
	def p_boundaryTop(x, on_bnd):
		return x[1] < 1 - DOLFIN_EPS and on_bnd
	
	
	bcs_u_LR = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), u_boundaryLeftRightWall)
	bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)),u_boundaryBottom)
	bc_p_T  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryTop)
	
	bcs = [bcs_u_LR, bc_u_B, bc_p_T]


	return a, L, UP, bcs


def second_part(a,L,UP, bcs):

	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
		
	U, P = UP.split()
	plot(U, interactive = True, title= "Velocity, vector valued function, approximated by FE")
	plot(P, interactive = True, title= "Pressure, scalar valued function, appriximated by FE")


def third_part(w, n = 20):               # The tube 
	mesh = UnitSquareMesh(n,n)       # Mesh over a unit square
	V = FunctionSpace(mesh, 'DG', 1) # Funcspace w DG as basisfunctions, i.e V = span{Ui*phi_i}, Phi_i DG basisfuncs for i = 0,..,n 
	#f  = Expression('(x[0] > w && x[0] < 1 - w ) ? 3.0 : 10.0', w=w)   # This is a CompiledExpression Object
	k = interpolate(f,V)                                               # This is a Function object k = sum(Ui*phi_i) and can be plotted 
	plot(k ,interactive = True, title = 'Tube sobdomain')


if __name__ == '__main__':
	from dolfin import * 
	w = 0.3
	a,L, UP, bcs = first_part(w = w)
	second_part(a, L,UP,bcs)
	
	#third_part(w = w)

