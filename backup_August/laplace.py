from dolfin import *

# dXX is a measure or volume integral, when it makes sense,
# for a triangle its area
# for a line its length
# for a cute its volute
def solve_poisson(mesh, dXX, f):
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    
    a = inner(grad(u), grad(v))*dXX + inner(u, v)*dXX # Using helmholtz equation so that it will be N.bcd in all around the rectangle
    L = inner(f, v)*dXX 
    
    A, b = assemble_system(a, L)                      # assebling the system in this manner avoid nan
    A.ident_zeros()                                   # because nan entries are substituted with zero 
    
    uh = Function(V)                                  # the solution of the Helmholtz equation by FE
    solve(A, uh.vector(), b)                          # the linear system solved by a lin.alg backend

    return uh                                         # the solution is returned 
    
    
# ---------------------------------------------------------------------------

# HW1: How do you handle boundary conditions? how to divide and give value to what is now marked
# HW2: What is the complexity of the algorithm. See WIKI for complexity! hos much does it cost, computational cost that is. 
# HW3: Can you make a more efficient algorithm? The for loops are not the problem here. There is other room for improvement in the script.

if __name__ == '__main__':
    import numpy as np
    
    N = 20
    full_mesh = RectangleMesh(Point(-1, -1), Point(1, 1), N, N)   # a rectangle mesh, the points are the two 'edges' 
    
    timer = Timer('marking')                                      # make a for loop with N in 8 16 32, refine the mesh and calc comp.cost
      
    active_cells = CellFunction('size_t', full_mesh, 0)           # give value 0 to all cells in the rectangle 
    
    width = 0.25                                                  # with for the tube
    active_domain = AutoSubDomain(lambda x, on_boundary: abs(x[0]) < width)   # func defining the active subdomain base on width from origin
    active_domain.mark(active_cells, 1)                                       # cells inside this widt are assigned value 1  
    # plot(active_cells, interactive=True)                                    # plot and see for yourself 
     
    edge_f = FacetFunction('size_t', full_mesh, 0)                            # FacetFunc marks all facets to 0, just intially to get started
    full_mesh.init(2, 1)						      # connectivity ????
    full_mesh.init(1, 2)                                                      # ?????

    ac_values = active_cells.array()                                          # print or plot this  
    for cell in SubsetIterator(active_cells, 1):                              # loop through all active cells with index 1  
        for edge in edges(cell):                                              # loop through all the edges 
            if np.sum(ac_values[edge.entities(2)]) == 1:                      # if an egde is connected to just one cell
                edge_f[edge] = 1                                              # assign to it value 1 
            
    print timer.stop()                                                        # clock how long that took 
    # print edge_f.array()               
    # plot(edge_f, interactive=True)
    
    
    # dXX = Measure('dx', domain=full_mesh, subdomain_data=active_cells, subdomain_id=1)    
    
    # print assemble(1*dXX)
    
    # f = Expression('1 + sin(4*x[0])*x[1]')
    # uh = solve_poisson(full_mesh, dXX, f)
    
    # plot(uh)
    # interactive()
    
    # print uh.vector().array()
    
