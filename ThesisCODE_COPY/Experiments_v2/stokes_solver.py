
from dolfin import * 



def alpha(u, K):
	C = Constant(1E5) # high value outside the vessel
	alpha= C*K*u
	return alpha







def stokes_solver(x0,R,K, mesh, W, speed, dpdy):
	
	
	u, p = TrialFunctions(W)
	v, q = TestFunctions(W)
		

	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, K), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(W)
	nu = Constant(1)
	dpdy = 2.1






	def inlet_boundary(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > (x0 - R - DOLFIN_EPS) and x[0] < (x0 + R + DOLFIN_EPS) #and on_bnd
  
	
	speed = 1.
	analyticalVelocity = Expression(["0","speed*0.25*dpdy*(-(x[0]-x0)*(x[0]-x0) + R*R)"], x0 = x0, R=R, dpdy=dpdy, speed = speed)
	inlet_velocity  = DirichletBC(W.sub(0), analyticalVelocity, inlet_boundary)
	bc = [inlet_velocity]
	
	
	#A, b = assemble_system(a, L, bc)
	#solve(A, UP.vector(), b, "lu")
	#U, P = UP.split()
	




	A, b = assemble_system(a, L, bc)
        solver = LUSolver(A)
	solver.solve(UP.vector(), b) 
	U, P = UP.split()  
	
	#plot(P, interactive=True, title = 'pressure')
        #plot(U, interactive=True, title = 'velocity')




	return U,P

