from dolfin import *


N = 40
w = 0.25

mesh = UnitSquareMesh(N, N, "crossed")
subdomains = CellFunction('size_t', mesh, 0)
# The CellFunction gives the each cell the value 0 
# When we make subdomains we give different values to the cells 

class Stem(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < 0.5 + DOLFIN_EPS and \
               0.5-w-DOLFIN_EPS < x[0] < 0.5+w+DOLFIN_EPS 

class LeftBranch(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 0.5 - DOLFIN_EPS and x[0] < 0.5 + DOLFIN_EPS\
               and not (x[1] < -x[0] + 1 - w)\
               and not (x[1] > -x[0] + 1 + w)

class RightBranch(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 0.5 - DOLFIN_EPS and x[0] > 0.5 - DOLFIN_EPS\
               and not (x[1] < x[0] - w)\
               and not (x[1] > x[0] + w)

Stem().mark(subdomains, 1)
LeftBranch().mark(subdomains, 1)
RightBranch().mark(subdomains, 1)


volume_stem = Constant(1)*dx(1, domain=mesh, subdomain_data=subdomains)
print assemble(volume_stem)

plot(subdomains)

y_mesh = SubMesh(mesh, subdomains, 1)

plot(y_mesh)

interactive()