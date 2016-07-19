
from dolfin import * 



# den opprinnelige K update som foelger den matematiske formelen 
def K_update(indikator, K):
	""" 
	For aa oppdatere K_new = K_gammel + skjaer_indikator 	
	"""
	
	#K.vector().axpy(-1, project(indikator, K.function_space()).vector())
	K_ny = Function(K.function_space())
        
	plot(indikator, interactive=True, title="indikator inside k update before proj")

	fine_mesh = refine(K.function_space().mesh())
        Vfine = FunctionSpace(fine_mesh, "CG", 1)
        plot(interpolate(indikator, Vfine), interactive=True, title="indikator on fine grid")

	indikator_proj = project(indikator, K.function_space())
        plot(indikator_proj, interactive=True, title="indikator inside k update, projisert paa romme av cg1 funksjoner")
	K_ny.vector()[:] = K.vector()[:] - indikator_proj.vector()[:] 


	return K_ny





# IKKE I BRUK 
# en manuell update som til slutt ikke vil ha behov for kapping 
def K_update_v2(indikator, K):

	
	#for i in indikator.vector().array():
	#	print 'dg1 indikator verdier', i 

	print 'lengde paa indikator vekter i dg1' ,len(indikator.vector()) # ca 6 000

	K_ny_funksjon = Function(K.function_space())
	#Pga projeksjonen er det mulig it den ikke blir noyaktig 1 eller 0 
	indikator_proj = project(indikator, K.function_space())
	

	
	K_ny      = K_ny_funksjon.vector().array()
	K_gammel  = K.vector().array()
	indikator = indikator_proj.vector().array() 


	
	print 'lengde paa indikator vekter i cg1' ,len(indikator) # versus ca 1000 


	assert(len(K_ny) == len(K_gammel) == len(indikator))
	


	over = []
	under = []
	for i in indikator:
			if abs(i) > 1E-4:
				 over.append(i)
			else:
				under.append(i)

	print 'len(over) =',len(over), 'len(under) =', len(under) 



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
	
	#print 'len(k_gammel) = ',len(K_gammel)
	#print 'type(K_gammel) = ',type(K_gammel)

	

	null = 1E-3

	for i in range(0,len(K_gammel)):
			
			# Foerst er vi et element i fluid domenet 
			if K_gammel[i] == null:
				# 0 - 1 ==> -1 ==> 0 fluid element kan ikke krympe til noe mindre, altsaa fluid element skal 
				# forbli fluid element 
				if (abs(indikator[i]) == 1):
					K_ny[i] = 0

				# 0 - (-1) ==> 1 fluid element skal bli solid element 
				if ( -1 -null < indikator[i] < -1 + null):
					K_ny[i] = 1 


				# 0 - 0 ==> var og forble et fluid element  
				if (abs(indikator[i]) <= null):
					K_ny[i] = 0. 	
				else:
					 print 'i fluid element, her skal vi aldri havne'			

			




			# Vi befinner oss i et solid element 
			if (K_gammel[i] == 1):
				
				# 1 - 1 ==> 0 solid element gaard over til aa bli fluid element 
				if (1-null<indikator[i] <1+null):
					K_ny[i] = 0. 


				# 1 - (-1) krymping av fluid domenet  = 1 +1 = 2 MEN
				# et solid element kan ikke bli MER solid  saa ==> 1 
				if (-1-null<indikator[i] < -1+null):
					K_ny[i]  = 1. 



				# dersom det er ingen endring i elementet 
				# 1 - 0 => 1 
				if (abs(indikator[i]) <= null):
					K_ny[i] = 1. 


				else:
					print 'i solid element, else setningen'

			else:
				print 'siste else statementen '
		
		
	K_ny_funksjon.vector()[:] = K_ny[:]
	return K_ny_funksjon 





def K_kappet(K):
	
	"""
	The isocontour is at 1/2
	"""
	
	# Tror det skal vaere absolutt verdi av verdienen 
	# men faar bad operand type for abs her 
	#K.vector()[abs(K.vector()) >= 0.5]  = 1.0   # alt over 0.5 er solig 
	#K.vector()[abs(K.vector()) < 0.5]   = 0     # mens alt under 0.5 er fluid 

	
	# prover en enkel loop
	# naar den er 0.5 oppfoerer den seg riktig, vokser saa krymper 
	new_level_lowSpeed = 0.55 # funker for speed 0.5 tilfellet 
	new_level_highSpeed =  0.60
	for i in range(0, len(K.vector())):
		if abs(K.vector()[i]) >= new_level_highSpeed:
			K.vector()[i] = 1.0
		else:
			K.vector()[i] = 0.0



	return K 

