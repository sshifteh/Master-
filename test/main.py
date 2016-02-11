# -*- coding: utf-8 -*-
from dolfin import *
from MyStokes import *
from ShapeOpt import *


dolfin.parameters.reorder_dofs_serial = False
dolfin.parameters.allow_extrapolation = True

#omega = UnitSquare(6,4)
omega = Mesh('stationary_circle.xml')
plot(omega)
mf=MeshFunction('size_t', omega, 'stationary_circle_facet_region.xml')

[J, u, p] = StokesSolve(omega, mf)
plot(u, interactive=True)

N = compute_normal_field(omega, mf)
plot(N, interactive=True)

print "Done"
