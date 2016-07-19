from dolfin import *
import numpy as np 

def wall_shear_stress_indicator(tau_xy, Dirac_delta):

	
	assert len(tau_xy.vector()) == len(Dirac_delta.vector())	

	interface_index = np.where(Dirac_delta.vector() > DOLFIN_EPS)
	print type(interface_index) # tuple av grad med index 0 saa dvs en liste som ikke er foranderlig  
	#print interface_index
	print len(Dirac_delta.vector()) # er over hele beregnings domenet
	# 4096 for interface og ellers 0 . bra. 



	print len(interface_index[0])
	vegg_skjaeret = tau_xy.vector()[interface_index[0]]
	print 'vegg_skjaeret', vegg_skjaeret 
	print max(abs(vegg_skjaeret)), min(abs(vegg_skjaeret)) 
		
		

	
	teoretisk_skjaer_verdi = 0.0525
	# terskel verdiene 20 prosent over og under   
	maks_terskel_verdi = teoretisk_skjaer_verdi + teoretisk_skjaer_verdi *(20./100) 
	min_terskel_verdi = teoretisk_skjaer_verdi - teoretisk_skjaer_verdi *(20./100 	
	assert min_terskel_verdi < maks_terskel_verdi



	# maa lage ekspansjon og kontraksjons objekter foerst tror jeg 
	ekspansjon = np.where[abs(vegg_skjaeret) > maks_terskel_verdi]
	kontraksjon = np.where[abs(vegg_skjaeret) < min_terskel_verdi]


	indikator = Function(tau_xy.function_space())
	


	epsilon = 1 
	indikator.vector() = np.where[ekspansjon = epsilon]  
 	indikator.vector() = np.where[kontraksjon = -epsilon]
	




	


wall_shear_stress_indicator_2(tau_xy, Dirac_delta):

	maks_terskel_verdi = 0.05775
	min_terskel_verdi = 0.04725 

	ekspasjon = np.where[abs(tau_xy_vector()) > max_terskel) and (Dirac_delta.vector() > DOLFIN_EPS)]
	kontraksjon = np.where[abs(tau_xy_vector()) < min_terskel) and (Dirac_delta.() > DOLFIN_EPS)] 
	
	
	indicator_function = Function(tau_xy.function_space()) 
	indicator_function_values =  indicator_function.vector().array()


	epsilon = 1 
	indicator_function_values[np.logical_and(tau_xy.vector() > maks_terskel_verdi, Dirac_delta.vector() > DOLFIN_EPS  )] = epsilon 
	indicator_function_values[np.logical_and(tau_xy.vector() < min_terskel_verdi, Dirac_delta.vector() > DOLFIN_EPS  )] = - epsilon 
	indicator_function_values[np.logical_and(tau_xy.vector() < maks_terskel_verdi , tau_xy.vector() > min_terskel_verdi, Dirac_delta.vector() > DOLFIN_EPS)] = 0
	

	# tror feilen ligger her 
	# for det foeste det vi ser av indikator funksjonen naa er ikke paa randen 
	
	 
	indicator_function.vector()[:] = indicator_function_values[:]
	

	return indicator_function

















wall_shear_stress_indicator_3(tau_xy, Dirac_delta):

	# ALternative 1 : SYNTAX ERROR 	next to where line 

	#ekspasjon = where(abs(tau_xy_vector) > max_terskel) and (Dirac_delta_vector > DOLFIN_EPS))
	#kontraksjon = where(abs(tau_xy_vecotr) < min_terskel) and (Dirac_delta_vector > DOLFIN_EPS)) 
	
        #indicator_function = Function(WSS.function_space())
        #indicator_function.vector()[growth] = 1
        #indicator_function.vector()[shrink] = -1



















wall_shear_stress_indicator_3(tau_xy, Dirac_delta):

	# Alternative 2: too many values to unpack ERROR  = saa for loop gaar altsaa ikke da. maa gjoere en numpy variant.
	#ekspansjon = []
	#kontraksjon = []
	#for i,j in (tau_xy_vector, Dirac_delta_vector):	
	# 	if (abs(tau_xy_vector[i] > maks_terskel_verdi) and Dirac_delta_vector[j] > DOLFIN_EPS ):
	#		ekspansjon.append(i)
	#	if (abs(tau_xy_vector[i] < min_terskel_verdi) and Dirac_delta_vector[j] > DOLFIN_EPS) : 
	#		kontraksjon.append(i)
	#	else: 
	#		pass 
	 


			
	#epsilon = 1 
	#indicator_function = Function(tau_xy.function_space())
	#indicator_function.vector()[ekspansjon] = epsilon
	#indicator_function.vector()[kontraksjon] = -epsilon 
	
	#plot(indicator_function, interactive = True, title = 'indicator function ')
	# Too many values to unpack ERROR Message 


	



	#return indicator_Function
	#return ekspansjon, kontraksjon
	
