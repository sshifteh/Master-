from dolfin import *

from plotting import plot

parameters['plotting_backend'] = "matplotlib"

from matplotlib import pyplot

    
mesh = UnitSquareMesh(16,16)
V = FunctionSpace(mesh, "CG", 1)
u = interpolate(Expression("x[0]"), V)

W = VectorFunctionSpace(mesh, "CG", 1)
w = interpolate(Expression(("x[0]", "1-x[1]")), W)

"""
plot(u)
pyplot.savefig('u')
pyplot.figure()
plot(mesh)
pyplot.savefig('mesh')
pyplot.figure()
plot(w)
pyplot.savefig('w')

pyplot.show()
"""

for i in range(3):
	plot(u)
	pyplot.suptitle("velocity: %g t" %i)
	pyplot.savefig("velocity: %g t" %i)
	pyplot.figure()

	plot(w)
	pyplot.savefig("mesh: %g t" %i)
	pyplot.show()
