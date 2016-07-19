
from alpha import alpha
from dolfin import * 

from plotting import plot
parameters['plotting_backend'] = "matplotlib"
from matplotlib import pyplot

# litteratur review and background 
# expemperiments , neumann condition everyw 


def stokes_solver(w, mesh, Vspace, Pspace, Wspace, K_array, n):
	"""
	Solves the Stokes equation with different bcs according to given level set function
	"""
	#FIXME argument parser on bsc_choices = [bcs_I, bcs_L, bcs_T, bcs_T, bcs_Y] 


	# DEFINE TEST AND TRIALFUNCTIONS 
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
		
	K_Func  = interpolate(K_array, Pspace) # Kontroll function

	#plot(K_Func, interactive = False, title = 'Control Domain to be solved over')
	#plot(K_Func)
	#pyplot.suptitle("K-Yshape")
	#pyplot.savefig("K-Yshape")
	#pyplot.show()	
	
	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, K_Func), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)




	# PARABOLIC VELOCITY PROFILE - POISEUILLE FLOW 	
	# CONSIDER HAVING THE BOUNDARY CONDITION DEFINED OVER THE ENTIRE WALL OF INLET AND OUTLET 
	# INCASE OF EXPANSION IN THOSE PLACES	
	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])	
	def u_boundaryBottom(x, on_bnd): 
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd

	# BOUNDARY CONDITIONS FOR T-SHAPED TUBE: 	
	def p_boundaryRightWall_T(x, on_bnd): 
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd 



	# test -------------
	
    	def left(x, on_boundary):
        	return near(x[0], 0.0)
	def right(x, on_boundary):
        	return near(x[0], 1.0)


	# I think this needs to be changes to only apply on the inlet .... 



	def bottom(x, on_boundary):
		return near(x[1], 0.0)


	# then we have this : 

	#def u_boundaryBottom(x, on_bnd): 
	#	return x[1] < DOLFIN_EPS #and x[0] > w and x[0] < 1-w and on_bnd



    	def top(x, on_boundary):
        	return near(x[1], 1.0)	


	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])
	bc_u_B  = DirichletBC(Wspace.sub(0), velocityFunc, bottom)


	bc_p_left  = DirichletBC(Wspace.sub(1), Constant(0.0), left)
	bc_p_right = DirichletBC(Wspace.sub(1), Constant(0.0), right)
	bcs_T = [bc_u_B, bc_p_left,bc_p_right]
	


	#----------------------------------------------------------------	
	# BOUNDARY CONDITIONS FOR I-SHAPED TUBE
	
    	def Bottom (x, on_boundary):
        	return near(x[1], 0.0)

    	def Top(x, on_boundary):
        	return near(x[1], 1.0)

  	
	#def u_boundaryBottom(x, on_bnd): 
	#	return x[1] < DOLFIN_EPS #and x[0] > w and x[0] < 1-w and on_bnd
	
	#def p_boundaryTop_I(x, on_bnd):
		#return x[1] < 1 - DOLFIN_EPS #and x[0] > w and x[0] < 1-w and on_bnd
	
	bc_u_B   = DirichletBC(Wspace.sub(0), velocityFunc, Bottom)
	bc_p_T   = DirichletBC(Wspace.sub(1), Constant(0.0),Top) 	
	bcs_I = [bc_u_B, bc_p_T]		

	#-----------------------------------------------------------------
 	
	# BOUNDARY CONDITION FOR AN L-SHAPED TUBE	
	#def u_boundaryBottom(x, on_bnd): 
	#	return x[1] < DOLFIN_EPS #and x[0] > w and x[0] < 1-w and on_bnd
		

	#def p_boundaryRightWall_L(x, on_bnd):
	#	return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd #red on the entire unitsquare boundary?S 
	

	def Right(x, on_boundary):
        	return near(x[0], 1.0)

	bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)),Bottom)
	bc_p_T  = DirichletBC(Wspace.sub(1), Constant(0.0), Right)
	bcs_L = [bc_u_B, bc_p_T]

	#-----------------------------------------------------------------	



	# BOUNDARY CONDITION FOR Y - SHAPED TUBE - OBS ONLY ON VELOCITY ! 
	# ALLOW FOR A NULSPACE THAT CONTAINS MORE THAN THE ZERO VECTOR	
	bc_u_B = DirichletBC(Wspace.sub(0), velocityFunc, Bottom)
	bcs_Y = [bc_u_B]
	
	

	# SOLVING A LINEAR SYSTEM OF EQUATIONS 
	A, b = assemble_system(a, L, bcs_Y)


	# NULSPACE FOR THE Y- SHAPED TUBE 
	nullspace = VectorSpaceBasis([interpolate(Constant((0.,0.,1.)),Wspace).vector()])
	as_backend_type(A).set_nullspace(nullspace) # Nul(A) = { c[0,0,1]^T}



	# SOLVE THE INDEFINITE MATRIX BY INVERTING IT AND USING LU FACTORIZATION
	solve(A, UP.vector(), b, "lu")
	U, P = UP.split()


	# PLOT OF VELOCITY AND PRESSURE 
	#plot(U, interactive = True, title= "Velocity")
	
	#plot(U)	
	#pyplot.suptitle("Y-shape:U: ")
	#pyplot.savefig("Y-shape:U: ")
	#pyplot.figure()	
	#pyplot.show()	




	#plot(P, interactive = True, title= "Pressure")

	

	# SAVE TO PARAVIEW 
	#file = File('Stokes_flow.xdmf')
	#file << U 


	return U,P, K_Func

