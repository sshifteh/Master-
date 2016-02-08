
from dolfin import *
from K_update import K_update
from capping import capping
from stokes_solver import stokes_solver
from WSS import WSS
from WSS_bdry import WSS_bdry, bdry
from indicator_funciton_conditions2 import ind_func


n = 50
w= 0.45
mesh = UnitSquareMesh(n,n) #, 'crossed')

DG1 =FunctionSpace(mesh, 'DG', 1)
DG0 = FunctionSpace(mesh, 'DG', 0)

Pspace = FunctionSpace(mesh, 'CG',1)
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace])

Kspace = DG0
WSSspace = DG1

K_1  = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)


# INTERESTING DOMAIN:
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 10./50 + DOLFIN_EPS and \
               10.0/50-DOLFIN_EPS < x[1] < 40.0/50+DOLFIN_EPS

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < 40./50 + DOLFIN_EPS and \
               10.0/50-DOLFIN_EPS < x[1] < 40.0/50+DOLFIN_EPS


# GENARAL CLASS FOR MAKING AN EXPRESSION OF ITEMS IN LIST
class SubDomainIndicator(Expression):
    def __init__(self, subdomain_list):
		self.subdomain_list = subdomain_list

    def eval(self, values, x):
        values[0] = 0.
        for subdomain in self.subdomain_list:
            if subdomain.inside(x, 0):
                values[0] = 1.
                break

subdomains = CellFunction('size_t', mesh, 1)
Bottom().mark(subdomains, 0)
#Top().mark(subdomains, 0)

interesting_domain = SubDomainIndicator([Bottom(), Top()])
plot(interesting_domain, mesh=mesh, interactive=True, title="Interesting domain")
interesting_domain_proj = project(interesting_domain, DG1)


def iterative(K_1):

	U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace,
                Pspace=Pspace, Wspace=Wspace, Kspace=Kspace, K_array=K_1, n=n )

	wss = WSS(U, WSSspace, mesh)
	plot(wss, interactive = True, title = 'wss')	
	bdry_ = bdry(mesh, K_Func, DG0)
	plot(bdry_, interactive = True, title = 'bdry_')

	ind_f = ind_func(bdry_, wss, interesting_domain_proj)
	        
	# Simon: indicator_function is currently in DG1, but your K_Func in DG0.
        # Hence you might need to interpolate indicator_function to DG0
	plot(ind_f, interactive=True, title='Indicator function')
	
	K_new = K_update(ind_func = ind_f, K_Func = K_Func)
	plot(K_new, interactive = True, title = 'K_updated ')
	K_capped = capping(K_new)
	plot(K_capped, interactive = False, title = 'K capped')
	stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace,
                Wspace=Wspace, Kspace=Kspace, K_array=K_capped, n=n )

	return K_capped


iterative(K_1=K_1)
for i in range(5):
	K_updated = iterative(K_1 = K_1)
	K_1= K_updated 



