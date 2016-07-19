from dolfin import * 

# take in wss and bdry instead of wss-bdry
# make a really simple function that takes k func and returns bdry

#def ind_func_expr(WSS_bdry):

def ind_func(bdry, WSS, DG1, interesting_domain, plotting = False):	
	# take the WSS_bdry function and make it into vector 
	# atm of getting it from the attributes in in DPspace and thats fine because so 
	# is our ind function supposed to be 

	#thresh_L = 1e-8; thresh_H = 2.0 #1000 # 1250 is the highest threshold. 1000 works as well. Everything equal or above 1300 doesnt give growth, i.e is to high a treshold.
	"""
	alist = WSS_bdry.vector().array() # copying the values into alist  	
	
	for a in range(len(alist)):
		if (abs(alist[a]) > thresh_H):
	
			alist[a] = -1 
					
		#elif (abs(alist[a]) < thresh_L):
		
		#	alist[a] = 1 
		
		#elif (tresh_L< alist[a] <thresh_H): 
		#	alist[a] = 0
		else:
			alist[a] = 0

	WSS_bdry.vector()[:] = alist[:]

	#plot(WSS_bdry, interactive = True, title = ' ind f calculated')	
	
	return WSS_bdry 
	""" 
	
	# originally functions made to arrays to be manipulated 
	# makes a copy 
	wss = WSS.vector().array()
	bdry_ = bdry.vector().array()
	interesting_domain_ = interesting_domain.vector().array()


	# 50% off: thresh_L = 0.025; thresh_H = 0.075 
	# 60% off: 	
	thresh_L = 0.005; thresh_H = 0.095 	

	# sizes math up
	print type(wss), type(bdry_)	
	print len(wss), len(bdry_), type(interesting_domain)

	# the copy is processed 
	M = bdry.vector().norm("linf") # Max verdien til bdry vector 
	for i in range(len(wss)):
		if abs(bdry_[i]) > 0.1*M :   # hvis storre enn 0.1 av max 
			if abs(wss[i]) > thresh_H:
				bdry_[i] = -1 * interesting_domain_[i]
				print 'hei' 

			elif abs(wss[i]) < thresh_L:
				bdry_[i] = 1 * interesting_domain_[i]
				print 'hadet' 
		
			else: 
				bdry_[i] = 0 * interesting_domain_[i]
	import numpy; print("{0:1.2f} % av bdry_ er null".format(100*float(len(numpy.where(bdry_ == 0)[0]))/len(bdry_)))
 	#from IPython import embed; embed()
	# it is put back into the original vector bdry	
	bdry.vector()[:] = bdry_[:]
	# it is projected onto functionspace DG1
	bdry_func = project(bdry, DG1) # denne projiseringen: tar gradienten av denne: den gradienten er mesh avhengig: mister oversikt over hvor grad er ulik null fordi testen var for streng. 
	# it is plotted
	print assemble((bdry-bdry_func)**2 * dx(domain = bdry.function_space().mesh()))
	if plotting: plot(bdry_func, interactive=True, title = 'new ind func')
	# it is returned 	
	return bdry_func 

	 



	
