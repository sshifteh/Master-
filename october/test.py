from dolfin import *

mesh = UnitIntervalMesh(10)

V = FunctionSpace(mesh, "CG", 1)

f = Function(V)

f.vector()[:] = 0
f.vector()[1] = 1

plot(f, interactive=True)

