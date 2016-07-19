from dolfin import * 

#def indicator_function(WSS, K_Func, DG0, DG1):
	
	# step one: you have to consider the bdry(dirac delta function) which is in DG0 
	# step two: you have to consider the WSS which is in DG1
	# step three: make the tresholds and the indicator function which is in DG1 


	# STEP ONE: considering the bdry --------------------------------------------------
	# the dirac delta function is in DG0, 1 node evaluation for each element 	

# make the mesh and spaces
n = 10  
mesh = intervalmesh(n)
DG0 = FunctionSpace(mesh, 'DG', 0)
DG1 = FunctionSpace(mesh, 'DG', 1) 
K_Func = Expression('sin(x[0])')
	
dirac_delta = grad(K_Func)**2
dirac_delta_Func = project(dirac_delta, DG0)
plot(dirac_delta_Func, interactive = True, title = 'dirac delta function viz in DG0')
# The dirac delta function is not infinity because the step of the Heaviside function is smoothed out by the 
# discretization. Its around 1000 . 

dirac_delta_vector = dirac_delta_Func.vector()
		
for c in dirac_delta_vector:
	if c != 0:
		print c # c is 2500  

# STEP TWO: consider the wss ------------------------------------------------------ 	
# wss is a Function object in the shear stress space = DG1, 3 node evaluations for each element 	

WSS_product = 0.5*(U[0].dx(1) + U[1].dx(0))
plot(WSS, interactive = True, title = 'WSS in DG1')
WSS_DG1_vector = WSS.vector()
	

#WSS_DG0 = project(WSS, DPspace)
#plot(WSS_DG0, interactive = True, title = 'WSS in DG0')	
#WSS_DG0_vector = WSS_DG0.vector()	

	
# Three : making treshold values and making the indicator function ---------------------------
threshold_low = 1e-8
threshold_high = 1250
			
# can make the indicator function in DG1 space 	
indicator = WSS_DG1_vector.array() # copy all the values into the indicator 



# interpolate always into the richer space 
# tips for debugging this for loop : 
# make a little example in 1d 
# interpolate the expression sin (x) , dont need to solve PDE, 
# this is actually  a postprocessing step
# plot the arrays and compare the values 




	for entry in range(len(indicator)):
		
		if dirac_delta_vector[entry] != 0 and abs(indicator[entry]) > threshold_high:
			indicator[entry] = -1
		
		elif dirac_delta_vector != 0 and abs(indicator[entry])< threshold_low:
			 indicator[entry] = 1
		elif dirac_delta_vector != 0: 
			indicator[entry] = 0
		else:
			pass  
			 
	indicator_function = Function(DG1)
	indicator_function.vector()[:] = indicator[:]
	
	plot(indicator_function, interactive = True, title = 'indicator function')



