from alpha import alpha
from dolfin import *


def stokes_solver(w, mesh, Vspace, Pspace, Wspace, Kspace, K_array, n):
	"""
	Solves the Stokes equation
	"""

	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)

	K_Func  = interpolate(K_array, Kspace) # Kontroll function
	#plot(K_Func, interactive = False, title = 'Control Domain to be solved over')

	f = Constant([0.0,0.0])
 	a = inner(alpha(u, K_Func), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)

	velocityFunc = Expression(["0","-(x[0]-0.45)*(x[0]-0.55)"])


	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd

	def restBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and  not (x[0]> w) and not (x[0] < 1-w) and on_bnd

    	def leftWall(x, on_boundary):
        	return near(x[0], 0.0)

	def rightWall(x, on_boundary):
        	return near(x[0], 1.0)

	def outlet(x, on_bnd):
		return x[1] < 1 - DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd


	# SIMON:
	# Project to the function space, plot it and see
	# plot the expression if provide the mesh argument
	# this didnt work
	# velocity_profile = project(velocityFunc, Pspace)
	# plot(velocity_profile, interactive = True, title = 'velocity parabolic profile')

	inlet  = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	restBottom = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), restBottom)
	right = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), rightWall)
	left = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), leftWall)

	#top = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), outlet)


	# Make dirichlet bcs overalt unntatt outletene
	#rightOutlet  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	#leftOutlet   = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryLeftWall)

	bcs = [inlet, restBottom, right, left] # right, left, top] #, rightOutlet, leftOutlet]


	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
	U, P = UP.split()
	n = FacetNormal(mesh)
	print "Mass conservation: ", assemble(inner(U, n)*ds)


	#file = File('Stokes_flow.xdmf')
	#file << U

	return U,P, K_Func

