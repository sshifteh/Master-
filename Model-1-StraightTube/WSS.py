from dolfin import *

def WSS(U, WSSspace, mesh ): #, interesting_domain):
	"""
	Calculates the WSS over the entire domain.
	Takes as input the velocity vector U and the
	DG1 space ShearStressSpace, and returns the WSS
	over the entire domain.

	"""

	#


	WSS_product = 0.5*(U[0].dx(1) + U[1].dx(0)) #*interesting_domain  # object type --> algebra.Product
	#print 'max wss:', WSS_product.vector().max(), 'min wss:', WSS_product.vector().min()




	# Glem dette :
	# projecting into vector function space instead of into function space
	#V_g = VectorFunctionSpace(mesh, 'Lagrange', 2, dim = 2)
	#w = TrialFunction(V_g)
	#v = TestFunction(V_g)

	#a = inner(w,v)*dx
	#L = inner(U[0].dx(1) + U[1].dx(0), v)*dx
	#grad_u = Function(V_g)
	#solve(a==L, grad_u)

	# Notice there are no essential bcs on this problem
	# dette vil ikke fungere fordi i L formen saa tar vi indreprodukt mellom en
	# skalar og en vektor , det gir en vektor aa integrere over.
	# det gir ikke energi.
	# basically : man kan ikke ta prikkprodukt mellom en vektor og en skalar.

	WSS = project(WSS_product, WSSspace) # In DG1, object type --> Function
	#plot(WSS, interactive = True, title = 'WSS') # plots the WSS over the entire

	return WSS
