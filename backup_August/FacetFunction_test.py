from dolfin import * 

# Using the facetfunction 

n = 2 
mesh = UnitSquareMesh(n,n)
plot(mesh, interactive = True, title = 'My mesh ')
V = FunctionSpace(mesh, 'Lagrange', 1)
print 'The dimension of the FunctionSpace V:', V.dim()

def Shifteh_boundary(x):
	return x[0] < DOLFIN_EPS 


bc = DirichletBC(V, Constant(59.0), Shifteh_boundary)
boundary = FacetFunction('size_t',bc)

plot(boundary, interactive = True, title = 'My FaceFunc')
# This didnt work!!!!!!!!!!!!!!1
