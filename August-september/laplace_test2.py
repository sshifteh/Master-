from dolfin import * 
import numpy as np

n = 5
mesh = UnitSquareMesh(n,n)
active_cells = CellFunction('size_t', mesh,0)					

widt = 0.25 									
active_domain = AutoSubDomain(lambda x, on_boundary: abs(x[0]) < widt)         
active_domain.mark(active_cells, 1)						
plot(active_cells, interactive = True)

edge_f = FacetFunction('size_t', mesh, 0)
mesh.init(2,1)                                                                 # assigning cells 3 values that is their numbered edges 
mesh.init(1,2)                                                                 # 2 for edge 1 for vertex, assigning one edge two vertices

ac_values = active_cells.array()
print ac_values 
for cell in SubsetIterator(active_cells, 1):
	for edge in edges(cell):
		print ac_values[edge.entities(2)]
		if np.sum(ac_values[edge.entities(2)]) == 1:
		 	edge_f[edge] = 1
print edge_f.array()
plot(edge_f, interactive = True)
#print timer.stop()


