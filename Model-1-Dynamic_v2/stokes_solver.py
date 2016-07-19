from alpha import alpha
from dolfin import *

def stokes_solver(x0,R, mesh, Vspace, Pspace, Wspace, K):
	"""
	Solves the Stokes equation with a additional inverse permeability term. 

	"""
	# Test and trial spaces over the Taylor-Hood elements  
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)

	# Source term
	f = Constant([0.0,0.0])
	
	# Parameters
	nu = Constant(0.00345) # Pascal x second 
	dpdz = Constant(-2.1)
	
	# Variational form 
 	a = inner(alpha(u, K), v)*dx + nu*inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)

	# Parabolic inlet velocity
        speed = 1
	#velocityFunc = Expression(["0","-speed*(x[0]-0.45)*(x[0]-0.55)"], speed=speed)
	analyticalVelocity = Expression(["0"," speed*0.25*nu*dpdz*((x[0]-x0)*(x[0]-x0) - R*R)"], x0=x0, R=R, nu=nu, dpdz=dpdz, speed = 					speed)


	# Defining the boundaries of the vessel 
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > (x0-R) - DOLFIN_EPS and x[0] < (x0+R) + DOLFIN_EPS and on_bnd

	def restBottom(x, on_bnd):

		return x[1] < DOLFIN_EPS and  not (x[0]> (x0-R) -DOLFIN_EPS) and not (x[0] < (x0+R) + DOLFIN_EPS) and 				on_bnd

    	def leftWall(x, on_boundary):
        	return near(x[0], 0.0)

	def rightWall(x, on_boundary):
        	return near(x[0], 1.0)

	def top_rest(x, on_bnd):
		return x[1] < 1 - DOLFIN_EPS and not (x[0] > (x0-R)) and not (x[0] < (x0+R)) and on_bnd

	inlet  = DirichletBC(Wspace.sub(0), analyticalVelocity, u_boundaryBottom)
	restBottom = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), restBottom)
	right = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), rightWall)
	left = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), leftWall)
	top = DirichletBC(Wspace.sub(0), Constant((0.0, 0.0)), top_rest )
	bcs = [inlet, restBottom, right, left, top] 

	# Solving the system of linear equations 
	A, b = assemble_system(a, L, bcs)
        solver = LUSolver(A)
	solver.solve(UP.vector(), b) # instead of savaing the solution in U immediatly to plot it successively
	U, P = UP.split()


	# Calculating mass conservation over the mesh to machine precision 
	n = FacetNormal(mesh)
	print "Mass conservation: ", assemble(inner(U, n)*ds)

	return U, P


if __name__ == "__main__":
	
	from main_iterative import w,mesh,Vspace,Pspace, Wspace, K 
	U,P = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K=K)

	plot(U)
	interactive()
