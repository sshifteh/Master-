from dolfin import * 

n = 5
mesh = UnitIntervalMesh(n)
#plot(mesh, interactive= True, title = 'myMesh')

V = FunctionSpace(mesh, 'Lagrange', 1) # V = span{p0,p1,p2,...,p5 }
# which is first order lagrange polynomial on the form, ax +b 
# p0 = c0 + c1x
# p1 = c2 + c3x etc 
# These are our basis functions which span out room for finding the solution geometrically.

P = FunctionSpace(mesh, 'DG', 1)
D = FunctionSpace(mesh, 'Lagrange', 2)


u = TrialFunction(V) # are out basis functions 
v = TestFunction(V)  # Are also out basis functionas perhaps for now

print V.__contains__(u) # False 	
print V.__eq__(P)       # False
print V.__eq__(D)       # False 
print V.__contains__(D) # False 
print D.__contains__(V) # False, this is strange
# D is a space of square functions, 
# V a space of linear functions 
# No, they dont overlap actually, its not a complete set. 


print V.ufl_element()
print 'components of the functionaspace V', V.component() # []
print 'dimension of the functionspace V', V.dim() # 6 

f = 1
#print FunctionSpace_dofmap(V )

#print V.FunctionSpace_element(1)

#interpolate = FunctionSpace_interpolate(V, f)


