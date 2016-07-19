from dolfin import * 



# den opprinnelige K update som foelger den matematiske formelen 
def K_update(indikator, K):
	""" 
	For aa oppdatere K_new = K_gammel + skjaer_indikator 	
	"""
	
	#K.vector().axpy(-1, project(indikator, K.function_space()).vector())
	K_ny = Function(K.function_space())
	#plot(indikator, interactive = True, title = 'indikator inside k update before proj')	


	fine_mesh = refine(K.function_space().mesh())
        Vfine = FunctionSpace(fine_mesh, "CG", 1)
        #plot(interpolate(indikator, Vfine), interactive=True, title="indikator on fine grid")



	
	indikator_proj = project(indikator, K.function_space())
	#plot(indikator_proj, interactive=True, title="indikator inside k update, projisert paa romme av cg1 funksjoner")	
	
	K_ny.vector()[:] = K.vector()[:] - indikator_proj.vector()[:] 


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
		if abs(K.vector()[i]) >= 0.65: #0.60:
			K.vector()[i] = 1.0
		else:
			K.vector()[i] = 0.0



	return K 

