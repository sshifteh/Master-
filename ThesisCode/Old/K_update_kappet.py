from dolfin import * 



# den opprinnelige K update som foelger den matematiske formelen 
def K_update(indikator, K):
	
	K_ny = Function(K.function_space())

	fine_mesh = refine(K.function_space().mesh())
        Vfine = FunctionSpace(fine_mesh, "CG", 1)
        #plot(interpolate(indikator, Vfine), interactive=True, title="indikator on fine grid")

	indikator_proj = project(indikator, K.function_space())
	#plot(indikator_proj, interactive=True, title="indikator inside k update, projisert paa romme av cg1 funksjoner")	
	
	K_ny.vector()[:] = K.vector()[:] - indikator_proj.vector()[:] 
	return K_ny





def K_cap(K):
	#The new isocontour is at 1/2
	for i in range(0, len(K.vector())):
		if K.vector()[i]>= 0.50: #0.20:#50: #65:#0: #0.65: #1: #1.5: #0.65: #0.60:
			K.vector()[i] = 1.0
		else:
			K.vector()[i] = 0.0
	return K 






