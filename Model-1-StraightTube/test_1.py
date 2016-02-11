def test_1():
	
	n = 10
	w= 0.45
	mesh = UnitSquareMesh(n,n,)
	Kspace = FunctionSpace(mesh, 'DG', 0)
	
		
	K  = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)
	K = interpolate(K, Kspace)				
	plot(K,interactive = True, title= 'K_old' )
	K_values = K.vector().array()

	# Which edges are such that the connected cells have different values
	facet_f = FacetFunction('size_t', mesh, 0)
	# Tabulate the indices correlating cell to facet (here: edge)	
	mesh.init(1, 2)

	# looping over all facets in the mesh  	
	for facet in facets(mesh):
		#print facet # 16 facets for a 2 by 2 mesh 
		connected_cells = facet.entities(2)
		# facet.entities(2) is the act of looking up the cells that share a facet 
		# (2) means topological dimension 2 = cell		
		# returns a list of indices of facets. If a cell has two facets, the list contains two indices.
		# if the cell is a boundary, if contains one indice 
		#print connected_cells 
		dofs = [Kspace.dofmap().cell_dofs(cell)[0] for cell in connected_cells]
		#print 'dofs:', dofs 
		# dofs actually contain the same info as connected_cells it seems.  
			
		if len(dofs) > 1:
			one, two = K_values[dofs]
			#print one, two # 1-1 is 0 its an inner in solid
			# 1-0 is on interface or on external boundary
			# no couldnt be on boundary because the len(dofs) > 1 
			# thus it is an internal point ! clever !  
			# with a 2 by 2 mesh this should give none of the facets. 
			#print abs(two-one)
			if abs(two-one) > 0:
				# this should give noe of the facets. 
				#print 'facet_f[facets]:', facet_f[facet]
				facet_f[facet] = 1
				
	
	# Now that I found them, which cells are connected to them
	# we have located the facets on the boundary.
	# which cells are connected to these facets.
	# connect cells to facets, i.e. transpose the mesh init table from above. 
	
	cell_f = CellFunction('size_t', mesh, 0)
	#plot(cell_f)
	# All cells in the mesh are given the value 0 
	facet_f_values = facet_f.array()
	#print 'facet_f_values',facet_f_values 
	#facet_f_calues is an vector with the value 1 only if the facet is in boundary facet(solid-fluid boundary not external boundary), 		else 0	
	mesh.init(2, 1)
	# transpose the table from before
	# facet indices as rows, cell indices in columns 

	for cell in cells(mesh):
		#print cell
		# prints mesh entity 198 of topological dimension one 
		# i.e. is accesses each cell in the mesh with a number 
		# each cell have 3 facets(in 2D = edges)
		# 
	 	#print cell.entities(1)
		# cell.entities(1) is the indices of the facets for each cell
		# so facet_f_values[cell.entities(2)] means go through the facet_f_values vector, and look up the indices for each facet 
		# correspond the value given to that facet with its number ... 
		#print facet_f_values[cell.entities(1)]
		#  ... sum the values give the indices of the facets. 
		# if they are > 0 , then they are a boundary facet. because 
		# only boundary facets are given value 0,
		# each cell has 3 facets, if the cell is in the solid or fluid domain it will have sum zero.
		# however if it is in the boundary it will have values 0 0 1. 
		# and the sum will be > 0
		if np.sum(facet_f_values[cell.entities(1)]) > 0:
			# when that is the case. give the cell_f which assigned 0 to all cells, 
			# give those cells the value 1.  
	        	cell_f[cell] = 1
	
	
	I = Function(Kspace)	
	assert len(K.vector().array())==len(I.vector().array()) 		
	I_kopi = I.vector().array()
	I_kopi[8]= -1  #shrinking 
	I_kopi[11]= -1 #shrinking 
	I_kopi[28]= -1 #shrinking 
	I_kopi[28]= -1 #shrinking		
	I_kopi[49]= 1  #growing 
	I.vector()[:] = I_kopi[:]
	print I.vector().array()
	plot(I,interactive = True, title= 'I')


	print cell_f.array()
	plot(cell_f,interactive = True, title = 'cell f')
	
	

	K.vector()[:] -= I.vector()[:]
	plot(K, interactive = True, title = 'new K')
		

	# So we have all cells connected to boundary equate 1. 
	# what if one cell was suppose to be included.
	# so one more cell would have to be given value 1. 
	# this is what I have to do. 	
	
	 	

if __name__ == '__main__':
	import numpy as np
	from dolfin import * 	
	test_1()

