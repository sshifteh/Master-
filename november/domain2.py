from dolfin import * 
import numpy as np
import time 

       

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


