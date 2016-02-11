
import time 
import numpy as np
import matplotlib.pyplot as plt 
from dolfin import * 	



class Mybase(object):
	def __init__(self, n):

		# Parameters 		
		self.n = n; 
		self.w = 0.3; 
		self.k = k 

		# Creating mesh and Function spaces 
		self.mesh = UnitSquareMesh(n,n)
		self.TensorSpace = TensorFunctionSpace(self.mesh, "DG", 1)	
		self.DPspace = FunctionSpace(self.mesh, 'DG', 0)	

		self.Vspace = VectorFunctionSpace(self.mesh, 'Lagrange', 2)
		self.Pspace = FunctionSpace(self.mesh, 'Lagrange', 1)
		self.Wspace = MixedFunctionSpace([self.Vspace, self.Pspace]) 
		
		# Creating test and trialfuncitons	
		self.u, self.p = TrialFunctions(self.Wspace)
		self.v, self.q = TestFunctions(self.Wspace)
		
	
		# The variational form   
		self.f = Constant([0.0,0.0]) 
	 	self.a = inner(self.alpha(self.u, k), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
		self.L = inner(self.f, self.v)*dx
		self.UP = Function(self.Wspace)

		
	def u_boundaryBottom(self, x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	 	
	def p_boundaryRightWall(self, x, on_bnd):
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd # why is it read on the entire unitsquare boundary? 
	
	
	def bcs(self):
		bc_u_B  = DirichletBC(self.Wspace.sub(0), Constant((0.0,1.0)), u_boundaryBottom)
		velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])
		bc_u_B  = DirichletBC(self.Wspace.sub(0), velocityFunc, u_boundaryBottom)
		bc_p_R  = DirichletBC(self.Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
		bcs = [bc_u_B, bc_p_R]
		return bcs 


	def solve(self):
		self.Wspace= Wspace
		self.a = a; self.L = L; self.UP = UP 	
		bcs = bcs()

		# System of linear equations 
		A, b = assemble_system(a, L, bcs)
		solve(A, UP.vector(), b, "lu")
		U, P = UP.split()
		plot(U, interactive = show, title= "Velocity")
		plot(P, interactive = show, title= "Pressure")

		return U,P
 
	
	def alpha(self, u, k):
		C = Constant(1e5) # high value outside the vessel
		#print type(k), type(u), 		
		return C*k*u



"""
	def  k_update(self, U):
			
		self.U = U #TODO perhaps U should be in the constructor
		#self.ShearStressSpace = ShearStressSpace
		self.DPspace = DPspace 				
 
		# The first K
		func  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7  )? 0.0:1.0', w=w)	
		k   = interpolate(func, self.Pspace) # control # Function object 
		#plot(k, interactive = show, title = 'Tube with CG elements')
		#plot(project(grad(k)**2, DPspace), interactive=True) # gradient of k is only nonzero on the bnd
		
		# The next K's based on WSS at bdry 
		shear_stress = project(0.5*(U[0].dx(1) + U[1].dx(0)), TensorSpace)    
		plot(shear_stress, interactive = True, title = 'WSS')	                 
		WSS_bdry = project(abs(shear_stress*grad(k)**2), DPspace).vector()

		indicator_function = Function(DPspace) 			        
		indicator_function_values =  indicator_function.vector().array() # TODO 

		d1 = 1e-5; d2 = 20 ; eps = 1.0        
	
		# Vectorized version : 
		start_time = time.time()	 	 	                  	
		indicator_function_values[WSS_bdry < d1] = -eps
		indicator_function_values[np.logical_and(WSS_bdry > d1, WSS_bdry < d2)] = 0
		indicator_function_values[WSS_bdry > d2] = eps 
		end_time = time.time()
		elapsed_time = end_time - start_time
		print 'time for the vectorized version is %.5f' %elapsed_time
	

		# Put values of temp indicator function array into indicator Function object 		
		indicator_function.vector().set_local(indicator_function_values)	
		plot(indicator_function, interactive = False, title = ' indicator function') # on the bdry 
	
		# Update step: axpy works in parallel 
		k.vector().axpy(-1, project(indicator_function, k.function_space()).vector())
		plot(k, interactive = True, title = 'updated k') # not yet capped 
	
		#plt.figure() 
		k_values = k.vector().array()  # make a copy to manipulate	
		#plt.plot(k_values, label = 'before')	
		k_values[k_values > 1] = 1
		k_values[k_values < 0] = 0
		#plt.plot(k_values, label = 'after')
		#plt.show(False)
		k.vector().set_local(k_values) # interpolate the values back into the function space 
	
""" 		
		



n = 50
Mybase(n = n)	
U, P = stokes_solver()
 	
	 










