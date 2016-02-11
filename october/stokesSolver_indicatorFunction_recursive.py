
	

def alpha(u, K):
	C = Constant(1e5) # high value outside the vessel
	alpha= C*K*u
	return alpha


	

def stokes_solver(w,mesh,Vspace,Pspace,Wspace,TensorSpace,DPspace,ShearStressSpace,K_array, n = 2 ):
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
		
	K_Func  = interpolate(K_array, Pspace) # control                     # Function object 
	plot(K_Func, interactive = True, title = 'Control Domain to be solved over')

	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, K_Func), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)
	# boundary conditions____________________________________________________________________________
	# function inside function called closure perhaps, good I think? #TODO  	
	def u_boundaryBottom(x, on_bnd): 
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd

	def p_boundaryRightWall(x, on_bnd): 
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd # why is it read on the entire unitsquare boundary? 

	bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)), u_boundaryBottom)
	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])
	bc_u_B  = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	bc_p_R  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	bcs = [bc_u_B, bc_p_R]
	# Solve system of linear equations_______________________________________________________________
	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
	U, P = UP.split()
	plot(U, interactive = True, title= "Velocity")
	plot(P, interactive = True, title= "Pressure")

	return U,P, K_Func

 
def domain(mesh, TensorSpace, DPspace,ShearStressSpace, Pspace, U,K_Func, n = 50):
       
	# TODO: why does it have to be here, inside the function to work?
	import time

	# TODO: 
	# Calculating the shear stress over the TensorSpace DG1
	# The difference with ShearStressSpace and TensorSpace
	# Why cant I use TensorSpace for the shear_stress ? 
	
	WSS_product = 0.5*(U[0].dx(1) + U[1].dx(0))  # object type --> algebra.Product 
	WSS = project(WSS_product, ShearStressSpace) # In DG1, object type --> Function 
	plot(WSS, interactive = True, title = 'WSS') # plots the WSS over the entire 	        
	
	bdry= grad(K_Func)**2 # object type --> tensoralgebra.Inner
	WSS_bdry = WSS*bdry   # object type --> algebra.Product
	WSS_bdry_Function =  project(WSS_bdry, DPspace)
	plot(WSS_bdry_Function, interactive = True, title = 'WSS_bdry_Function') # plots the WSS over the bdry 
	
	# TODO: In fact WSS_bdy is our new indicator function I think? 
	# TODO: indicator_function_values =  indicator_function.vector().array() 
	indicator_Function = Function(DPspace) # in DG0, 0 to begin with 	
	indicator_Function_values =  indicator_Function.vector().array()

	# Defining treshold values for which the WSS is 'significant' to be considered
	d1 = 1e-5; d2 = 20 ; eps = 1.0 
	
	# TODO: should have abs of wss, but it doesnt work 	       
	# Vectorized version if test:  
	start_time = time.time()	 	                  		
	indicator_Function_values[WSS_bdry_Function.vector() < d1] = -eps
	indicator_Function_values[np.logical_and(WSS_bdry_Function.vector() > d1, WSS_bdry_Function.vector() < d2)] = 0
	indicator_Function_values[WSS_bdry_Function.vector() > d2] = eps 
	end_time = time.time()
	elapsed_time = end_time - start_time
	print 'time for the vectorized version is %.5f' %elapsed_time

	indicator_Function.vector().set_local(indicator_Function_values) #same as copying by [:]	
	plot(indicator_Function, interactive = True, title = 'Indicator function after having been updated')	
	
	# The update step: axpy works in parallel 
	K_Func.vector().axpy(-1, project(indicator_Function, K_Func.function_space()).vector())
	plot(K_Func, interactive = True, title = 'The first updated K on the boundary')
	


	# TODO: destroys everything 
	# vectorized capping 	
	#K_Func.vector()[K_Func.vector() > 1] = 1
	K_Func.vector()[K_Func.vector() < 0.5] = 0
	K_updated = project(K_Func, K_Func.function_space())
	plot(K_updated, interactive = True, title = 'K updated in the domain')		

	# TODO: an idea
	#w = 0.3
	#K_1 = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7)? 0.0:1.0', w=w) #K inmath.formulation	
	#K_old = project(K_1, Pspace) 


	# TODO:
	# beware: mixing of meshes. wont work in paralell . base evr. on one mesh. 
	# this recursive works!
	# however, one problem is that there should be some capping perhaps,
	# beacuse now it calculated velocty only in the extended parts(bdry)

	return K_updated



def domain2(Vspace, mesh, TensorSpace, DPspace,ShearStressSpace, Pspace, U,K_Func, n ):
       
	_n = FacetNormal(Vspace.mesh())    # FacetNormal ufl 
	
        T = dot((grad(U) + grad(U).T), _n) # tensoralgebra ufl, finding the wss
        Tn = dot(T, _n)                    # tensoralgebra ufl, wss on bdry 
        Tt = T - Tn*_n                     # algebra sum ufl, ???   
	
	print ''
	print Tn, ':::::::::::::::::::::::::::::::::::::::::'
	print ''	
	print ''
	print Tt, '*************************'
	print ''



	print 'type _n=', type(_n),'type T=',type(T), 'Tn=',type(Tn), 'type Tt=', type(Tt)  
	
        #tau_form = dot(self.v, Tt)*ds()
        #assemble(tau_form, tensor=self.tau.vector())

	# capping 
        #self.b[self._keys] = self.tau.vector()[self._values] # FIXME: This is not safe!!!
        #get_set_vector(self.b, self._keys, self.tau.vector(), self._values, self._temp_array


	

def main():

	# Parameters
	n = 40
	w   = 0.3
	show = True 
	
	# Mesh and Functionspaces
	mesh = UnitSquareMesh(n,n)
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
		
	ShearStressSpace=FunctionSpace(mesh, 'DG', 1)
	TensorSpace = TensorFunctionSpace(mesh, "DG", 1) # Space for the shearstess,could be named ShearStressSpace			
	DPspace = FunctionSpace(mesh, 'DG', 0)	         # Space for the gradient of the wss 
	
	# Expression for the initial T shape 
	K_1 = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7  )? 0.0:1.0', w=w) #K in math.formulation	
	
	# Function calls
	U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace,Pspace=Pspace, Wspace=Wspace, TensorSpace=TensorSpace, 			DPspace=DPspace, ShearStressSpace=ShearStressSpace , K_array=K_1, n=n)

	K_Func = domain(mesh = mesh, TensorSpace=TensorSpace, DPspace=DPspace, ShearStressSpace=ShearStressSpace, Pspace=Pspace, U=U, K_Func = K_Func, n=n)

	stokes_solver(w,mesh,Vspace,Pspace,Wspace,TensorSpace,DPspace,ShearStressSpace, K_array=K_Func, n=n)
 	

	

	#domain2(Vspace=Vspace, mesh=mesh, TensorSpace=TensorSpace, DPspace=DPspace,ShearStressSpace=ShearStressSpace, Pspace=Pspace, U=U,K_Func=K_Func, n = 40)
	
	

if __name__ == '__main__':
	# Import statements	
	import time 
	import numpy as np
	import matplotlib.pyplot as plt 
	from dolfin import * 	
	main()

	
