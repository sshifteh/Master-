from dolfin import * 

def compute(nx,ny, degree):

	mesh = UnitSquareMesh(nx,ny)
	plot(mesh)	
	V = FunctionSpace(mesh, 'Lagrange', degree = degree)
	bc = DirichletBC(V, 0.0, 'on_boundary')

	u_exact = Expression('sin(pi*x[0])*sin(pi*x[1])')
	f = 2*pi*pi*u_exact

	v = TestFunction(V)
	u = TrialFunction(V)
	a = dot(grad(u), grad(v))*dx
	L = f*v*dx
		
	u = Function(V)
	solve(a ==L, u, bc)
	
	solve(a == L, u, bc)
	solver_parameters = {'linear_solver': 'cg', 'preconditioner': 'ilu'}
	info(parameters, True)

	plot(u)

	Ve = FunctionSpace(mesh, 'Lagrange', 5)
	u_e_Ve = interpolate(u_exact, Ve)
	u_Ve = interpolate(u, Ve)
	e_Ve = Function(Ve)
	e_Ve.vector()[:] = u_e_Ve.vector().array()- u_Ve.vector().array()
	error = e_Ve**2*dx
	E =  sqrt(assemble(error))

	interactive()
	#E = errornorm(u, u_exact, Ve)
	return E 


if __name__ == "__main__":
	degree = 1 
	h = []
	E = []
	for nx in [4, 8, 16]: #, 32, 64, 128, 264]:
		h.append(1.0/nx)
		E.append(compute(nx,nx, degree))


	# convergence rates 
	from math import log as ln 
	for i in range(1, len(E)):
		r = ln(E[i]/E[i-1]) / ln(h[i]/h[i-1])
		print "h=%10.2E  r=%.2f" %(h[i], r)	
	 
	
"""
#for i in [8,16,32]:
	mesh = refine(mesh)
	plot(mesh, interactive = True)
""" 
