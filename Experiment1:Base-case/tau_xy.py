from dolfin import * 
import numpy as np 

def tau_xy(u):
	
	tau_xy = 0.5*(u[0].dx(1) + u[1].dx(0))	
	return tau_xy 	





def Dirac_delta(K):
	# Alternativ 1 	
	Dirac_delta = grad(K)**2
	plot(Dirac_delta, interactive = True, title = 'test test ')
	absolutt_Dirac_delta = (abs(Dirac_delta))
	#print test  
	#plot(absolutt_Dirac_delta,interactive = True, title = 'test')
	#File('abssolutt_Dirac_delta.xdmf' ) << absolutt_Dirac_delta
	#FIXME : 
	#interface_proj = project(interface, dg1)	
	#Det er 3 alternativer for rommet til interface, eller vent kanskje konverteres til et annet rom for plotting 
	# dg0 funker 
	# cg1 funker 
	# dg1 funker 
	return Dirac_delta












def wall_shear_stress(tau_xy, Dirac_delta):

	
	assert len(tau_xy.vector()) == len(Dirac_delta.vector())	

	#for i in Dirac_delta.vector().array():
	#	print i
	# 4096

 
	##interface_index = np.where(Dirac_delta.vector() > 1E-3) # DOLFIN_EPS
	##print type(interface_index) # tuple av grad med index 0 saa dvs en liste som ikke er foranderlig  
	##print interface_index
	##print len(Dirac_delta.vector()) # er over hele beregnings domenet # nesten 25 000 
	# 4096 for interface og ellers 0 . bra. 

	

	##print len(interface_index[0]) # nesten 800 
	##vegg_skjaeret = tau_xy.vector()[interface_index[0]]
	#print 'vegg_skjaeret', vegg_skjaeret 
	##print 'maks vegg skjaer = ', max(abs(vegg_skjaeret)), 'min vegg skjaer = ', min(abs(vegg_skjaeret)) 
		
	##return vegg_skjaeret 	
	
	

# denne fungerer, men mangler absolutt verdi av skjaret . vi maa ta absolutt verdi av skjaret foer vi vurderer det 
def indikator_test1(space, Dirac_delta, tau_xy):


	teoretisk_skjaer_verdi = 0.0525
	modell_skjaer_verdi_plotOverLine = 0.042	
	prosent = float(50) 
	maks_terskel_verdi = modell_skjaer_verdi_plotOverLine + modell_skjaer_verdi_plotOverLine *(prosent/100) 
	min_terskel_verdi = modell_skjaer_verdi_plotOverLine - modell_skjaer_verdi_plotOverLine *(prosent/100) 	
	assert min_terskel_verdi < maks_terskel_verdi

	
	indikator = Function(space)
	epsilon = 1  

	
	for i in range(0,len(indikator.vector())):
		if Dirac_delta.vector()[i] > 1E-4:

			if tau_xy.vector()[i] > maks_terskel_verdi:
				indikator.vector()[i] = epsilon  


			if tau_xy.vector()[i]  < min_terskel_verdi:
				indikator.vector()[i] = -epsilon 


			else:
				indikator.vector()[i] = 0 
		

	return indikator 








# 3 endringer: 
# 1. legger til interesting domain 
# 2. gjoer om til vektorer 
# 3. tar absolutt verdi av skjaret  
def indikator_test2(space, Dirac_delta, tau_xy, interesting_domain):
	

	teoretisk_skjaer_verdi = 0.0525
	modell_skjaer_verdi_plotOverLine = 0.042		
	prosent = float(20) 

	maks_terskel_verdi = modell_skjaer_verdi_plotOverLine + modell_skjaer_verdi_plotOverLine *(prosent/100) 
	min_terskel_verdi = modell_skjaer_verdi_plotOverLine - modell_skjaer_verdi_plotOverLine *(prosent/100) 	
	assert min_terskel_verdi < maks_terskel_verdi


	skjaer    = tau_xy.vector().array()
	vegg      = Dirac_delta.vector().array()
	
	indikator = Function(space)	
	indikator_array = indikator.vector().array()
	domene    = interesting_domain.vector().array()	


	epsilon = 1 
	for i in range(0,len(indikator_array)):
		if abs(vegg[i]) >  1E-4: # eller 100 

			if abs(skjaer[i]) >= maks_terskel_verdi:
				indikator_array[i] = epsilon * domene[i]


			elif abs(skjaer[i]) <= min_terskel_verdi:
				indikator_array[i] = -epsilon * domene[i]


			else:
				indikator_array[i] = Constant(0) * domene[i]
 	
		else:
			pass  

		
	indikator.vector()[:] = indikator_array[:]

	return indikator 

	















