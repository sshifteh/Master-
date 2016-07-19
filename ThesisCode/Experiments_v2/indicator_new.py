from dolfin import * 
import numpy as np 

from cell_index_tests import cell_index, neighbouring_cells, cell_dofs_by_entity_dim


def tau_xy(u):
	return 0.5*(u[0].dx(1) + u[1].dx(0))	
	

def interface(K):
	# absolute value taken because K is projected nad might have - 1E-16 values 
	return sqrt(abs(grad(K)**2))
	



def indicator(space, interface, tau_xy, interesting_domain, epsilon):
	

	teoretisk_shear_value = 0.0525
	model_shear_value_plotOverLine = 0.042 #0.023#0.042		
	percent =float(20) #(30) 
	

	maks_terskel_verdi = model_shear_value_plotOverLine + model_shear_value_plotOverLine *(percent/100) 
	min_terskel_verdi = model_shear_value_plotOverLine - model_shear_value_plotOverLine *(percent/100) 	
	assert min_terskel_verdi < maks_terskel_verdi



	skjaer = tau_xy.vector().array()
	vegg   = interface.vector().array()
	
	assert len(skjaer) == len(vegg)
	
	indikator = Function(space)	
	indikator_array = indikator.vector().array()
	domene    = interesting_domain.vector().array()	

	print 'model WSS'
	print 'maks_terskel_verdi', maks_terskel_verdi	
	print 'min_terskel_verdi' , min_terskel_verdi
	#epsilon = 2
	for i in range(0,len(indikator_array)):
                #print i, vegg[i], skjaer[i], min_terskel_verdi, maks_terskel_verdi, domene[i] #, len(skjaer), len(indikator_array), len(vegg), len(domene), 

		if abs(vegg[i]) >  0.1: #0.1: #1: # 0.1: <-- bra for n = 16 x 16  #35: #1E-4:
			#from IPython import embed; embed()
			if abs(skjaer[i]) > maks_terskel_verdi:
                                # We want to grow og vi har dofen i til cellen som skal vokse
				# vi indeksen til den cellen 
                                cell_id = cell_index(i, space)  # fant celle indeksen 
				# vi finner naboen til alle grensende celler fra meshet. Siden cellenes sidekanter i meshet grenser mot hverandre 
				# men ikke i funksjonen. 
				# vi finner indeksen paa nabo celler 
                                indeks_naboceller = neighbouring_cells(cell_id, space.mesh())
                                


				for cell in indeks_naboceller :  # la oss gaa gjennom alle nabo cellene vi har funnet  
                                    dofs = list(cell_dofs_by_entity_dim(cell,  2, space)) # of finne dofene til hver av de 
                                    dofs += list(cell_dofs_by_entity_dim(cell, 2, space)) # of putte de inni en liste
                                    dofs += list(cell_dofs_by_entity_dim(cell, 2, space))

                                    #print dofs

                                    for i in dofs: 
                                        indikator_array[i] = epsilon * domene[i]

				indikator_array[i] = epsilon * domene[i]
                                
	
			elif abs(skjaer[i])  < min_terskel_verdi:
				#indikator_array[i] = -epsilon * domene[i]
                                
				cell_id = cell_index(i,space)
				index_naboceller = neighbouring_cells(cell_id, space.mesh())
				for cell in index_naboceller:
					dofs = list(cell_dofs_by_entity_dim(cell,  2, space)) # of finne dofene til hver av de 
                                    	dofs += list(cell_dofs_by_entity_dim(cell, 2, space)) # of putte de inni en liste
                                    	dofs += list(cell_dofs_by_entity_dim(cell, 2, space))

				for i in dofs:
					indikator_array[i] = -epsilon * domene[i] 


     				indikator_array[i] = -epsilon * domene[i]

                                 
	
			else:
				indikator_array[i] = 0 * domene[i]
 	                        #print "test 3"
	
		
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
	

