from dolfin import *

def WSS(U):

	WSS = 0.5*(U[0].dx(1) + U[1].dx(0)) #*interesting_domain  # object type --> algebra.Product
	#print 'max wss:', WSS_product.vector().max(), 'min wss:', WSS_product.vector().min()
	# Miro : this is not wall shear stress, it is just stress 
	

	return WSS
