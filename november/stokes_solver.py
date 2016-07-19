from alpha import alpha
from dolfin import * 
	

def stokes_solver(w, mesh, Vspace, Pspace, Wspace, K_array, n):
	"""
	Solves the Stokes equation
	"""


	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
		
	K_Func  = interpolate(K_array, Pspace) # Kontroll function
	#plot(K_Func, interactive = False, title = 'Control Domain to be solved over')

	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, K_Func), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)


	def u_boundaryBottom(x, on_bnd): 
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	def p_boundaryRightWall(x, on_bnd): 
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd 

	#bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)), u_boundaryBottom)
	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])
	

	# Project to the function space, plot it and see
	# plot the expression if provide the mesh argument 
	
	# this didnt work 
	#velocity_profile = project(velocityFunc, Pspace)
	#plot(velocity_profile, interactive = True, title = 'velocity parabolic profile')	

	bc_u_B  = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	bc_p_R  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	bcs = [bc_u_B, bc_p_R]


	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
	U, P = UP.split()
	plot(U, interactive = True, title= "Velocity")
	#plot(P, interactive = True, title= "Pressure")


	#file = File('Stokes_flow.xdmf')
	#file << U 


	return U,P, K_Func

