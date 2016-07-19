from alpha import alpha
from dolfin import *

def stokes_solver(w, mesh, Vspace, Pspace, Wspace, K):
	"""
	Solves the Stokes equation with inverse permeability term. 
	"""

	# Test and trial spaces over the Taylor-Hood elements  
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
	# Source term
	f = Constant([0.0,0.0])


	
	# Parameters 
	mu = Constant(0.00345)  # Dynamic ?viscosity, units: Pascal x second 
	dpdz = Constant(-2.1)   # Pressure gradient set for the analytical solution 
	R = w #Constant(0.05)   # This is a major mistake !!!! look at the verify code 
	
		
	
	# Variational form 
 	a = inner(alpha(u, K), v)*dx + nu*inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)
	
	
	
	# Parabolic inlet velocty Dirichlet condition 
        speed = 10
	#velocityFunc = Expression(["0","-speed*(x[0]-0.45)*(x[0]-0.55)"], speed=speed)
	analyticalVelocity = Expression(["0"," 0.25*mu*dpdz*((x[0]-0.50)*(x[0]-0.50) - R*R) "], R=R, dpdz=dpdz)


	# Defining the boundaries of the vessel 
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > w - DOLFIN_EPS and x[0] < 1-w + DOLFIN_EPS and on_bnd




	def restBottom(x, on_bnd):

		return x[1] < DOLFIN_EPS and  not (x[0]> w -DOLFIN_EPS) and not (x[0] < 1-w + DOLFIN_EPS) and 			on_bnd

    	def leftWall(x, on_boundary):
        	return near(x[0], 0.0)

	def rightWall(x, on_boundary):
        	return near(x[0], 1.0)

	# should not enforece strong pressure boundary conditions
	# that is taken care of weakly in the pressure/lazy/ do nothing boundary condition 
	# def outlet(x, on_bnd):
	#	return x[1] < 1 - DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd


	def top_rest(x, on_bnd):
		return x[1] < 1 - DOLFIN_EPS and not (x[0] > w) and not (x[0] < 1-w) and on_bnd

	inlet  = DirichletBC(Wspace.sub(0), analyticalVelocity, u_boundaryBottom)
	restBottom = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), restBottom)
	
	right = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), rightWall)
	left = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), leftWall)
	top = DirichletBC(Wspace.sub(0), Constant((0.0, 0.0)), top_rest )

	#top = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), outlet)
	# Make dirichlet bcs overalt unntatt outletene
	#rightOutlet  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	#leftOutlet   = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryLeftWall)

	bcs = [inlet, restBottom, right, left, top] # right, left, top] #, rightOutlet, leftOutlet]



	# Solving the system of linear equations 
	A, b = assemble_system(a, L, bcs)
        solver = LUSolver(A)
	solver.solve(UP.vector(), b)   # instead of savaing the solution in U immediatly to plot it successively
	U, P = UP.split()


	# Calculating mass conservation over the mesh to machine precision 
	n = FacetNormal(mesh)
	print "Mass conservation: ", assemble(inner(U, n)*ds)


	# Calculating the flux 

	return U, P


if __name__ == "__main__":
	
	from main_iterative import w,mesh,Vspace,Pspace, Wspace, K 
	U,P = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K=K)

	plot(U)
	interactive()
