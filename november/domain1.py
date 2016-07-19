from dolfin import * 
import numpy as np
import time 

def domain(mesh, TensorSpace, DPspace,ShearStressSpace, Pspace, U,K_Func, n = 50):
       
	# TODO: 
	# Calculating the shear stress over the TensorSpace DG1
	# The difference with ShearStressSpace and TensorSpace
	# Why cant I use TensorSpace for the shear_stress ? 
	
	WSS_product = 0.5*(U[0].dx(1) + U[1].dx(0))  # object type --> algebra.Product 
	WSS = project(WSS_product, ShearStressSpace) # In DG1, object type --> Function 
	plot(WSS, interactive = False, title = 'WSS') # plots the WSS over the entire 	        
	
	bdry= grad(K_Func)**2 # object type --> tensoralgebra.Inner
	WSS_bdry = WSS*bdry   # object type --> algebra.Product
	WSS_bdry_Function =  project(WSS_bdry, DPspace)
	plot(WSS_bdry_Function, interactive = False, title = 'WSS_bdry_Function') # plots the WSS over the bdry 
	
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
	
	K_Func.vector()[K_Func.vector() > 1] = 1
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


