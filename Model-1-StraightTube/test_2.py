
def test_2():
	n = 50
	w = 0.45
	mesh = UnitSquareMesh(n,n,)
	Kspace = FunctionSpace(mesh, 'DG', 0)
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'CG', 1)	
	Wspace = MixedFunctionSpace([Vspace, Pspace])

	
	K  = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)

	u, p, K = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, Kspace=Kspace, K_array=K, n=n)
	
	plot(u)

	#u = Expression(["0", "-10*(x[0]-0.45)*(x[0]-0.55)"])	
	#u = Function(Vspace)
	
	WSS_product = 0.5*(u[0].dx(1) + u[1].dx(0))	
	WSS_space = FunctionSpace(mesh, 'CR', 1)
	
	
	WSS = project(WSS_product, WSS_space)
	
	facet_f = FacetFunction('size_t', mesh, 0)
	plot(facet_f)
	
	
	for facet in facets(mesh):
		#print facet # 16 facets for a 2 by 2 mesh 
		connected_cells = facet.entities(2)
		# facet.entities(2) is the act of looking up the cells that share a facet 
		# (2) means topological dimension 2 = cell		
		# returns a list of indices of facets. If a cell has two facets, the list contains two indices.
		# if the cell is a boundary, if contains one indice 
		#print connected_cells 
		dofs = [WSS_space.dofmap().cell_dofs(cell)[0] for cell in connected_cells]
		print dofs 




if __name__ == '__main__':
	from dolfin import * 	
	from stokes_solver import stokes_solver
	

	test_2()
	

	
