
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

	# ---------------------------Boundary conditions for the I -shaped vessel ---------------------------------- 
	
	# Defining left wall 
    	def left(x, on_boundary):
        	return near(x[0], 0.0)
	
	#Defining right wall 
	def right(x, on_boundary):
        	return near(x[0], 1.0)
	
	# Defining ENTIRE bottom 
	def bottom(x, on_boundary):
		return near(x[1], 0.0)

	# Defining inlet area, PART of bottom  
	def u_boundaryBottom(x, on_bnd): 
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	

	def rest_bottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and  not x[0]> w and not x[0] < 1-w and on_bnd


	# Defining ENTIRE top 
    	def top(x, on_boundary):
        	return near(x[1], 1.0)	

	# Defining the outlet area, PART of top  
	def p_boundaryTop_I(x, on_bnd):
		return x[1] < 1 - DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	
	velocityFunc = Expression(["0","-(x[0]-0.45)*(x[0]-0.55)"])
		
	# Experiment 0 : parabolic velocity on entire bottom + p = 0 on entire top 
	entire_bottom = DirichletBC(Wspace.sub(0), velocityFunc, bottom)
	entire_top = DirichletBC(Wspace.sub(1), Constant(0.0), top)
	bcs_0 = [entire_bottom, entire_top] 
	

	# Experiment 1 : velocity ONLY on inflow + p = 0 ONLY on outflow area, Neumann else
	inlet_bottom = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	outlet_top   = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryTop_I) 	
	bcs_1 = [inlet_bottom, outlet_top]		
	


	# Experiment 2 : velocity on inflow and bottom + 0 on rest +  P on outlet 
	inlet_bottom = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	rest_bottom = DirichletBC(Wspace.sub(0), Constant((0,0)), rest_bottom)	
	#oulet_top = DirichletBC(Wspace.sub(0), Constant(0.0), p_boundaryTop_I) 		
	leftWall = DirichletBC(Wspace.sub(0), Constant((0,0)), left)		
	rightWall = DirichletBC(Wspace.sub(0), Constant((0,0)), right)	

	bcs_2 = [inlet_bottom, rest_bottom, leftWall, rightWall]

	

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

