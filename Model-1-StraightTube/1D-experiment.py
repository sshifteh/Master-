from dolfin import *
import numpy

parameters["reorder_dofs_serial"]=False
mesh = UnitIntervalMesh(8)
DG1 = FunctionSpace(mesh, 'DG', 1)
DG1. dofmap().tabulate_all_coordinates(mesh) # ? 

u = Function(DG1)
#u.vector()[0:8] = 1 # fungerte ikke 

a = numpy.zeros(16)
a[0:8] = 1
u.vector()[:] = a 
plot(u, interactive = True, title = 'u')
CG1 = FunctionSpace(mesh, 'CG',1)
v = project(u, CG1)
plot(v, interactive = True, title = 'v = u proj into CG1')
w = project(grad(v)**2, DG1)

plot(w, interactive = True, title = 'w proj into DG1')

print 'w proj into DG1', w.vector().array()

mesh_fine = UnitIntervalMesh(100)
V = FunctionSpace(mesh_fine, 'CG',1)
plot(interpolate(w,V), interactive = True, title='w is proj onto CG1')
print ''
print 'w proj into cG1', w.vector().array()
print ''

print 'normalized '
print w.vector().array()/w.vector().norm('linf')

y = Function(DG1)
a = numpy.zeros(16)
a[numpy.where(w.vector().array()/w.vector().norm('linf') > 0.1)] = 1 
y.vector()[:] = a 

plot(interpolate(y,V), interactive = True, title = 'y')
print''
print y.vector().array()

# now I do the same with the WSS 
#wss = 0.5*(v.dx(0)) #*interesting_domain  # object type --> algebra.Product 
#print wss.vector().array()
# not possible to calculate wss in 1D. 


