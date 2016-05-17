from dolfin import * 


def ind_func_expr(WSS_bdry):
	
	# take the WSS_bdry function and make it into vector 
	# atm of getting it from the attributes in in DPspace and thats fine because so 
	# is our ind function supposed to be 

	
	alist = WSS_bdry.vector().array() # copying the values into alist  


	tresh_L = 1e-8; tresh_H = 1000 # 1250 is the highest threshold. 1000 works as well. Everything equal or above 1300 doesnt give growth, i.e is to high a treshold.
	
	for a in range(len(alist)):
		if (abs(alist[a]) > tresh_H):
	
			alist[a] = -1 
					
		#elif (abs(alist[a]) < tresh_L):
		
		#	alist[a] = 1 
		
		#elif (tresh_L< alist[a] <tresh_H): 
		#	alist[a] = 0
		else:
			alist[a] = 0

	WSS_bdry.vector()[:] = alist[:]

	#plot(WSS_bdry, interactive = True, title = ' ind f calculated')	
	
	return WSS_bdry 
	
