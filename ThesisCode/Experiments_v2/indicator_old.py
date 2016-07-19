from dolfin import * 
import numpy as np 
#from dofmapper import cell_index, neighbouring_cells






def indicator_old(space, Dirac_delta, tau_xy, interesting_domain):
	

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
                #print i, vegg[i], skjaer[i], min_terskel_verdi, maks_terskel_verdi, domene[i] #, len(skjaer), len(indikator_array), len(vegg), len(domene), 

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

	




