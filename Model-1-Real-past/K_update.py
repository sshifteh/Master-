
from dolfin import * 



# den opprinnelige K update som foelger den matematiske formelen 
def K_update(indikator, K):
	""" 
	For aa oppdatere K_new = K_gammel + skjaer_indikator 	
	"""
	
	#K.vector().axpy(-1, project(indikator, K.function_space()).vector())
	K_ny = Function(K.function_space())
	indikator_proj = project(indikator, K.function_space())
	K_ny.vector()[:] = K.vector()[:] - indikator_proj.vector()[:] 


	return K_ny



# en manuell update som til slutt ikke vil ha behov for kapping 
def K_update_v2(indikator, K):
	K_ny = Function(K.function_space())
	indikator_proj = project(indikator, K.function_space())
	
	K_ny_array      = K_ny.vector().array()
	K_gammel_array  = K.vector().array()
	indikator_array = indikator_proj.vector().array() 

	assert(len(K_ny_array) == len(K_gammel_array) == len(indikator_array))

	""" 
	for i in len(K_gammel_array):
		
		# fluid domenet skal vokse og allerede fluid elementer forblid uforandret
		# Da er vi i et fluid element og skal forbi et fluid element  0 <== -1 <== 0 - 1 
		if (K_gammel_array[i] = 0 and indikator_array[i] = 1) : 
		 		K_ny_array[i] = 0



 		# da er vi i et fluid element forsatt og har krymping 1 <== 0 - (-1) 
		if (K_gammel_array[i] = 0 and indikator_array[i] = -1) : 
				K_ny_array[i] = 1

		if (K_gammel_array[i] = 0 and indikator_array[i] = 0 ):
				K_ny_array[i] = 0
			 
	""" 		
	
	print 'len(k_gammel_array) = ',len(K_gammel_array)
	print 'type(K_gammel_array) = ',type(K_gammel_array)


	for i in range(0,len(K_gammel_array)):
			
			# Foerst er vi et element i fluid domenet 
			if K_gammel_array[i] == 0 :
				# 0 - 1 ==> -1 ==> 0 fluid element kan ikke krympe til noe mindre, altsaa fluid element skal 
				# forbli fluid element 
				if (indikator_array[i] == 1):
					K_ny_array[i] = 0

				# 0 - (-1) ==> 1 fluid element skal bli solid element 
				if (indikator_array[i] == -1):
					K_ny_array[i] = 1 


				# 0 - 0 ==> var og forble et fluid element  
				if (indikator_array[i] == 0):
					K_ny_array[i] = 0 	
				else:
					pass #print 'i fluid element, her skal vi aldri havne'			

			
			# Vi befinner oss i et solid element 
			if (K_gammel_array[i] == 1):
				
				# 1 - 1 ==> 0 solid element gaard over til aa bli fluid element 
				if (indikator_array[i] == 1):
					K_ny_array[i] = 0 


				# 1 - (-1) krymping av fluid domenet  = 1 +1 = 2 MEN
				# et solid element kan ikke bli MER solid  saa ==> 1 
				if (indikator_array[i] == -1):
					K_ny_array[i]  = 1 



				# dersom det er ingen endring i elementet 
				# 1 - 0 => 1 
				if (indikator_array[i] == 0):
					K_ny_array[i] = 1 
				else:
					pass #print 'i solid element, else setningen'

			else:
				pass #print 'siste else statementen '
		
		
	K_ny.vector()[:] = K_ny_array[:]
	return K_ny 



def K_kappet(K):
	
	"""
	The isocontour is at 1/2
	"""
	
	# Tror det skal vaere absolutt verdi av verdienen 
	# men faar bad operand type for abs her 
	#K.vector()[abs(K.vector()) >= 0.5]  = 1.0   # alt over 0.5 er solig 
	#K.vector()[abs(K.vector()) < 0.5]   = 0     # mens alt under 0.5 er fluid 

	
	# prover en enkel loop
	for i in range(0, len(K.vector())):
		if abs(K.vector()[i]) >= 0.60:
			K.vector()[i] = 1.0
		else:
			K.vector()[i] = 0.0



	return K 

