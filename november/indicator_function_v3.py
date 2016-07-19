from dolfin import * 

def indicator_function(WSS, K_Func, DG0, DG1):
	
	

	# STEP ONE: considering the bdry --------------------------------------------------
	# the dirac delta function is in DG0, 1 node evaluation for each element 	
	
	dirac_delta = grad(K_Func)**2
	dirac_delta_Func = project(dirac_delta, DG0)
	#plot(dirac_delta_Func, interactive = True, title = 'dirac delta function viz in DG0')
	# The dirac delta function is not infinity because the step of the Heaviside function is smoothed out by the 
	# discretization. Its around 1000 . 

	# Since I want to do check each entry in the for loop, the dirac delta function has to be interpolated into the 
	# richer space, DG1
	dirac_delta_Func_DG1 = interpolate(dirac_delta_Func, DG1)
	dirac_delta_vector = dirac_delta_Func_DG1.vector()
	


	# STEP TWO: consider the wss ------------------------------------------------------ 	
	#plot(WSS, interactive = True, title = 'WSS in DG1')
	WSS_DG1_vector = WSS.vector()


	# Three : making treshold values and making the indicator function ---------------------------
	threshold_low = -9 #1e-8
	threshold_high = 4  #1250
			
	# Making the indicator function in DG1 space 	
	indicator = WSS_DG1_vector.array() # copy all the values into the indicator 
	dirac_delta_array = dirac_delta_Func_DG1.vector().array()

	# Testing lengths of basises 
	print 'no of basisfunction for the dirac delta vector', len(dirac_delta_vector)
	print 'no of nodes for the indicator ', len(indicator)
	print 'no of nodes for the wss vector', len(WSS_DG1_vector)
	# 15000 for all three
	



	# interpolate always into the richer space 
	# tips for debugging this for loop : 
	# make a little example in 1d 
	# interpolate the expression sin (x) , dont need to solve PDE, 
	# this is actually  a postprocessing step
	# plot the arrays and compare the values 


	#from IPython import embed;
	#embed()	

	# I think I found the problem, dirac_delta_vector is a vector and cannot do numpy operations 
	# over it. 
	
	for entry in range(len(indicator)):		
		

		if dirac_delta_array[entry] != 0 and abs(indicator[entry]) > threshold_high:	

			indicator[entry] = -1
		
		elif dirac_delta_array[entry] != 0 and abs(indicator[entry])< threshold_low:

			indicator[entry] = 1

		elif dirac_delta_array[entry] != 0: 

			indicator[entry] = 0

		else:
			pass  
			 
	# maybe stuff wasnt copied ? yes it should have been copied with this syntax 

	# now careful the indicator which is updated is an array object.
	# is was a copy of the wss data
	# project it so it becomes an Function object and plot it 
	
	

	indicator_function = Function(DG1)
	indicator_function.vector()[:] = indicator[:]
	
	plot(indicator_function, interactive = True, title = 'indicator function')

	return 



