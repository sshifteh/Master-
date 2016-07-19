from dolfin import * 
n = 50 
w = 0.3

mesh = UnitSquareMesh(n,n)
V = FunctionSpace(mesh, 'CG', 1)
f = Expression('(x[0] > w && x[0] < 1- w )? 0.0:1.0', w=w)
f = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 && x[0] > w )? 0.0:1.0', w=w)
k = interpolate(f,V)
plot(k, interactive = True) 
