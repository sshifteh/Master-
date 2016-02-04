from dolfin import * 

# take in wss and bdry instead of wss-bdry
# make a really simple function that takes k func and returns bdry

#def ind_func_expr(WSS_bdry):

def ind_func(bdry, WSS, DG1, interesting_domain, plotting = False):	
	# take the WSS_bdry function and make it into vector 
	# atm of getting it from the attributes in in DPspace and thats fine because so 
	# is our ind function supposed to be 

	
	# originally functions made to arrays to be manipulated 
	# makes a copy 
	wss = WSS.vector().array()/WSS.vector().norm('linf')
	bdry_ = bdry.vector().array()
	interesting_domain_ = interesting_domain.vector().array()
	
	
	#wss_max = WSS.vector().norm('linf')
	#print wss_max 
	#wss_normalized = wss/wss_max
	#print 'heeeeeeeeeeeeeeeeeeeeeer',wss_normalized 

	# 50% off: thresh_L = 0.025; thresh_H = 0.075 
	# 60% off: 	
	thresh_L = 0.00025; thresh_H = 0.13 	

	# sizes math up
	print type(wss), type(bdry_)	
	print len(wss), len(bdry_), type(interesting_domain)

	# the copy is processed 
	M = bdry.vector().norm("linf") # Max verdien til bdry vector 
	for i in range(len(wss)):
		if abs(bdry_[i]) > 0.1*M :   # hvis storre enn 0.1 av max 
			if abs(wss[i]) > thresh_H:
				bdry_[i] = -1 * interesting_domain_[i]
				#print 'hei' 

			elif abs(wss[i]) < thresh_L:
				print ' lower treshold', wss[i]
				bdry_[i] = 1 * interesting_domain_[i]
				#print 'hadet' 
		
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

	 



	
