from dolfin import *

dolfin.parameters.reorder_dofs_serial = False
dolfin.parameters.allow_extrapolation = True

def vertex_normal(mesh, boundary_parts, i):
    ver = Vertex(mesh, i)
    n = Point(0.0, 0.0)
    div = 0.0
    for fac in entities(ver, mesh.geometry().dim()-1):
        f = Facet(mesh, fac.index())
        if f.exterior()==True:
            n+=f.normal()
            div+=1
    n/=div
    return n

def compute_normal_field(omega, mf):
    #Boundary Mesh
    gamma = BoundaryMesh(omega, "exterior")
    mapa = gamma.entity_map(0)
    
    V = FunctionSpace(gamma, "CG", 1)
    N = VectorFunctionSpace(gamma, "CG", 1)
    normal_field = Function(N)
    
    normal_x = Function(V)
    normal_y = Function(V)
    
    
    for ver in entities(gamma, 0):
        i = mapa[ver.index()]
        point = vertex_normal(omega, mf, i)
        normal_x.vector()[ver.index()] = -1.0*point[0]
        normal_y.vector()[ver.index()] = -1.0*point[1]
    
    (x, y) = TestFunctions(N)
    solve(inner(normal_field[0], x)*dx + inner(normal_field[1], y)*dx - inner(normal_x, x)*dx -inner(normal_y, y)*dx == 0, normal_field)
    
    #Dirichlet Boundary must be initialized by Function in Vector Space over Omega (only on Gamma fails)
    V_vec = VectorFunctionSpace(omega, "CG", 1)
    normal_fieldV = Function(V_vec)
    for ver in entities(gamma, 0):
        i = mapa[ver.index()]
        normal_fieldV.vector()[i] = normal_field.vector()[ver.index()]
        normal_fieldV.vector()[i+omega.num_vertices()] = normal_field.vector()[ver.index()+gamma.num_vertices()]
    
    deform = TrialFunction(V_vec)
    v = TestFunction(V_vec)
    a = 0.01*inner(nabla_grad(deform), nabla_grad(v))*dx + 0.01*inner(deform,v)*ds(4)
    L = inner(Constant((0.0,0.0)),v)*dx
    bc1 = DirichletBC(V_vec, normal_fieldV, mf, 4)
    bc2 = DirichletBC(V_vec, Constant((0,0)), mf, 2)
    bc3 = DirichletBC(V_vec, Constant((0,0)), mf, 3)
    bc4 = DirichletBC(V_vec, Constant((0,0)), mf, 1)
    bc = [bc1, bc2, bc3, bc4]
    deform = Function(V_vec)
    solve(a==L, deform, bcs=bc)
    return deform
