from alpha import alpha
from dolfin import * 


def stokes_solver(w, mesh, Vspace, Pspace, Wspace, K_array, n):

	#DEFINE TEST AND TRIALFUNCTIONS__________________________________________________________________________ 
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)	
	K_Func  = interpolate(K_array, Pspace) # Kontroll function
	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, K_Func), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)

	#BCS____________________________________________________________________________________________________	

	# ---------------------------Boundary conditions for the T -shaped vessel ---------------------------------- 

	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])	
	def top(x, on_boundary):
        	return near(x[1], 1.0)	
	
	def bottom(x, on_boundary):
		return near(x[1], 0.0)

	def inletBottom(x, on_bnd): 
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	
	def restBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and  not x[0]> w and not x[0] < 1-w and on_bnd


	def p_boundaryRightWall(x, on_bnd): 
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd 
	
	def p_boundaryLeftWall(x, on_bnd):
		return x[0] < DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd
	
	
	# Experiment 1 : velocity ONLY on inflow + p = 0 ONLY on outflow area, Neumann else
	inlet_bottom = DirichletBC(Wspace.sub(0), velocityFunc, inletBottom)
	leftOutlet = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryLeftWall)
	rightOutlet = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	rest_bottom = DirichletBC(Wspace.sub(0), Constant((0,0)), restBottom)			
	
	bcs_2 = [inlet_bottom, rest_bottom, leftOutlet, rightOutlet]

	
	# SOLVING A LINEAR SYSTEM OF EQUATIONS____________________________________________________ 
	A, b = assemble_system(a, L, bcs_2)

	# SOLVE THE INDEFINITE MATRIX BY INVERTING IT AND USING LU FACTORIZATION__________________
	solve(A, UP.vector(), b, "lu")
	U, P = UP.split()

	# PLOT OF VELOCITY AND PRESSURE 
	plot(U, interactive = True, title= "Velocity")
	#plot(P, interactive = True, title= "Pressure")
	#plot(K_Func,interactive = True, title = 'K_func')	

	

	return U,P, K_Func






	
