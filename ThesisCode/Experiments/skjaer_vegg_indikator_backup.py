from dolfin import * 
import numpy as np 
#from dofmapper import cell_index, neighbouring_cells



def tau_xy(u):
	tau_xy = 0.5*(u[0].dx(1) + u[1].dx(0))	
	return tau_xy 	













def Dirac_delta(K):

	Dirac_delta = sqrt(abs(grad(K)**2))
	return Dirac_delta












def indikator_test(space, Dirac_delta, tau_xy, interesting_domain):
	

	teoretisk_skjaer_verdi = 0.0525
	modell_skjaer_verdi_plotOverLine = 0.023#0.042		
	prosent =float(10) #(30) 
	

	maks_terskel_verdi = modell_skjaer_verdi_plotOverLine + modell_skjaer_verdi_plotOverLine *(prosent/100) 
	min_terskel_verdi = modell_skjaer_verdi_plotOverLine - modell_skjaer_verdi_plotOverLine *(prosent/100) 	
	assert min_terskel_verdi < maks_terskel_verdi


	skjaer    = tau_xy.vector().array()
	vegg      = Dirac_delta.vector().array()
	
	assert len(skjaer) == len(vegg)
	
	indikator = Function(space)	
	indikator_array = indikator.vector().array()
	domene    = interesting_domain.vector().array()	

	
	
	epsilon = 2
	for i in range(0,len(indikator_array)):
                print i, vegg[i], skjaer[i], min_terskel_verdi, maks_terskel_verdi, domene[i] #, len(skjaer), len(indikator_array), len(vegg), len(domene), 

		if abs(vegg[i]) >  0.1: #0.1: #1: # 0.1: <-- bra for n = 16 x 16  #35: #1E-4:
			#from IPython import embed; embed()
			if abs(skjaer[i]) > maks_terskel_verdi:
				indikator_array[i] = epsilon * domene[i]
                                print "test 1"
	
			elif abs(skjaer[i])  < min_terskel_verdi:
				indikator_array[i] = -epsilon * domene[i]
                                print "test 2"                                 
	
			else:
				indikator_array[i] = 0 * domene[i]
 	                        print "test 3"
	
		
	indikator.vector()[:] = indikator_array[:]
	return indikator 

	






def wall_shear_stress(tau_xy, Dirac_delta):

	
	assert len(tau_xy.vector()) == len(Dirac_delta.vector())	

	#for i in Dirac_delta.vector().array():
	#	print i
	# 4096 
	##interface_index = np.where(Dirac_delta.vector() > 1E-3) # DOLFIN_EPS
	##print type(interface_index) # tuple av grad med index 0 saa dvs en liste som ikke er foranderlig  
	##print interface_index
	print len(Dirac_delta.vector()) # er over hele beregnings domenet # nesten 25 000 
	# 4096 for interface og ellers 0 . bra. 
	##print len(interface_index[0]) # nesten 800 
	##vegg_skjaeret = tau_xy.vector()[interface_index[0]]
	#print 'vegg_skjaeret', vegg_skjaeret 
	##print 'maks vegg skjaer = ', max(abs(vegg_skjaeret)), 'min vegg skjaer = ', min(abs(vegg_skjaeret)) 
		
	##return vegg_skjaeret 	
	


