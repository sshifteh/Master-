from dolfin import *


N = 40
w = 0.25

mesh = UnitSquareMesh(N, N, "crossed")
subdomains = CellFunction('size_t', mesh, 1)
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


class SubDomainIndicator(Expression):
    def __init__(self, subdomain_list):
		self.subdomain_list = subdomain_list

    def eval(self, values, x):
        values[0] = 1.
        for subdomain in self.subdomain_list:
            if subdomain.inside(x, 0):
                values[0] = 0.
                break

    
Stem().mark(subdomains, 0)
LeftBranch().mark(subdomains, 0)
RightBranch().mark(subdomains, 0)


t = SubDomainIndicator([Stem(), LeftBranch(), RightBranch()])
plot(t, mesh = mesh, interactive)

volume_stem = Constant(1)*dx(1, domain=mesh, subdomain_data=subdomains)
print assemble(volume_stem)

plot(subdomains)

#y_mesh = SubMesh(mesh, subdomains, 1)
#plot(y_mesh)

interactive()
