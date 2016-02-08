from dolfin import *

def WSS(U, WSSspace, mesh ): #, interesting_domain):

	WSS_product = 0.5*(U[0].dx(1) + U[1].dx(0)) #*interesting_domain  # object type --> algebra.Product
	#print 'max wss:', WSS_product.vector().max(), 'min wss:', WSS_product.vector().min()

	WSS = project(WSS_product, WSSspace) # In DG1, object type --> Function
	#plot(WSS, interactive = True, title = 'WSS') # plots the WSS over the entire

	return WSS
