from dolfin import * 

n = 5
mesh = UnitSquareMesh(n,n)
plot(mesh, interactive = True, title = 'My mesh over the domain')              # Mesh over the domain, i.e. a unit square
active_cells = CellFunction('size_t', mesh,0)					
print 'type of active cells, created using CellFunction',type(active_cells)    # CellFunction

plot(active_cells, interactive = True, title = 'active cells, cellvalue: 0')   # All cells in the mesh given value 0
										
widt = 0.25 									
active_domain = AutoSubDomain(lambda x, on_boundary: abs(x[0]) < widt)         # swig object, defines what is boundary
print active_domain		
						
active_domain.mark(active_cells, 1)						
#plot(active_domain, interactive = True, title = ' active_domain')		
plot(active_cells, interactive = True, title = ' active_cells, cellvalue: 1')  # orange part. Autosub(ac_dom) gave bcs on active_cell to give 

edge_f = FacetFunction('size_t', mesh, 0)
plot(edge_f, interactive = True, title = 'edge_f')                             # every facet is marked originally
mesh.init(2,1)                                                                 # assigning cells 3 values that is their numbered edges 
mesh.init(1,2)                                                                 # 2 for edge 1 for vertex, assigning one edge two vertices
                                                                               # and assigning every vertex one edge and saving in memory
# recall, we have applied active_domain to active_cells, so active_cells are now the cells inside the width
# we make it an array 
print 'active_cells: ',active_cells 
ac_values = active_cells.array() # MeshFunction object so we have to make it into an array object 
print 'ac_values: ',ac_values    # [1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0]
#for cell in SubsetIterator(a)

"""
    active_domain = AutoSubDomain(lambda x, on_boundary: abs(x[0]) < width) # func defining the active subdomain base on width from origi
    active_domain.mark(active_cells, 1)                                       # cells inside this widt are assigned value 1  
    # plot(active_cells, interactive=True)                                    # plot and see for yourself 
     
    edge_f = FacetFunction('size_t', full_mesh, 0)                            # FacetFunc marks all facets to 0, just intially to get started
    full_mesh.init(2, 1)						      # connectivity ????
    full_mesh.init(1, 2)                                                      # 2 = 

    ac_values = active_cells.array()                                          # print or plot this  
    for cell in SubsetIterator(active_cells, 1):                              # loop through all active cells with index 1  
        for edge in edges(cell):                                              # loop through all the edges 
            if np.sum(ac_values[edge.entities(2)]) == 1:                      # if an egde is connected to just one cell
                edge_f[edge] = 1                                              # assign to it value 1 
            
    print timer.stop() 
"""

