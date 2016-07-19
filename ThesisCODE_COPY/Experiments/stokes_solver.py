
from dolfin import * 



def alpha(u, K):
	C = Constant(1E5) # high value outside the vessel
	alpha= C*K*u
	return alpha




def stokes_solver(x0,R,K, mesh, W, speed, mu, dpdy):
	
	
	u, p = TrialFunctions(W)
	v, q = TestFunctions(W)
		

	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, K), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(W)
	


	def inlet_boundary(x, on_bnd):
		return (x[1] - (-1)) < DOLFIN_EPS  and x[0] > - R  and x[0] <  R  
  

	
	# CHANGE MADE HERE : removed the minus infront x[0] fungerte fint det 
	analyticalVelocity = Expression(["0","speed*0.25*dpdy*(-(x[0])*(x[0]) + R*R)"], R=R, dpdy=dpdy, speed = speed)

	inlet_velocity  = DirichletBC(W.sub(0), analyticalVelocity, inlet_boundary)
	bc = [inlet_velocity]
	
	
	A, b = assemble_system(a, L, bc)
        solver = LUSolver(A)
	solver.solve(UP.vector(), b) 
	U, P = UP.split()  
	
	

	return U,P 


