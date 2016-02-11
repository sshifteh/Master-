from dolfin import *

 

n = 2
mesh = UnitSquareMesh(n,n)
mesh2 = UnitIntervalMesh(n)

V = FunctionSpace(mesh, 'Lagrange', 1)
W = VectorFunctionSpace(mesh, 'Lagrange', 2)

print 'dimension of W, n = 2, dim: ', W.dim()
print 'the dofmap of V, p1, n=2', V.element()


#V = VectorFunctionSpace(mesh, 'Lagrange', 1)
#W = VectorFunctionSpace(mesh, 'Lagrange', 2)
#P = FunctionSpace(mesh, 'Lagrange',1)





#print P.dim() # 3 nodes(x) times 3 nodes(y) = 9  
#print V.dim() # 3 nodes(x) times 3 nodes(y) times 2 evaluations = 18 
#print W.dim() # 50 

 






 
