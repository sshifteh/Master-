from dolfin import * 


def capping(K_Func):
	"""
	Creating a 'homogenous' fluid domain, 
	in the sense that the fluid domain is level set 0,
	while solid domain is level set 1.
	Not sure of the terminology though. 
	"""	
	

	
	alist = K_Func.vector().array()  # copi of the K_Func vector data 
	print alist # vector of size 2601 	

	for a in range(len(alist)):
		#print alist[a]
		if alist[a] >= 0.5:
			alist[a] = 1.0

		elif alist[a]  < 0.5:
			alist[a] = 0

		# dont think I even need this line, either they are above or below a half
		else:
			alist[a] = alist[a]

	# test for negative values 
	#for a in alist:
	#	if a<0:
	#		print a
		
	
	"""
	vectorized version FIXME 	
	K_Func.vector()[K_Func.vector() > 1.0]  = 1.0
	K_Func.vector()[K_Func.vector() < 0] = 0
	K_updated = project(K_Func, K_Func.function_space())
	#plot(K_updated, interactive = True, title = 'K capped')		
	"""	
	
	# Put the copi back into K_Func
	K_Func.vector()[:]= alist[:]

	# again test for negative values 	
	for k in K_Func.vector():
		if k <0:
			print k 
		
	K_updated = project(K_Func, K_Func.function_space())
	
	return K_updated 


	
	
