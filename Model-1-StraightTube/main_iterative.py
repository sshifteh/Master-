
from dolfin import *
from K_update import update_K
from capping import capping
from stokes_solver import stokes_solver
from WSS import WSS
from WSS_bdry import bdry
from indicator_funciton_conditions2 import ind_func

show_plot = True
def plot(*args, **kwargs):
    if not show_plot:
        return
    dolfin.plot(*args, **kwargs)

n = 40
w = 0.25
mesh = UnitSquareMesh(n,n) #'crossed')

CG1 =FunctionSpace(mesh, 'CG', 1)
DG1 =FunctionSpace(mesh, 'DG', 1)
DG0 = FunctionSpace(mesh, 'DG', 0)

Pspace = CG1
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace])

Kspace = DG0
WSSspace = DG1

K_expr = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)
K = interpolate(K_expr, Kspace)

u = Function(Vspace)

def iterative(K, iters):

    for i in range(iters):
	# Solve for velocity and pressure
	U, P = stokes_solver(w, mesh, Vspace, Pspace, Wspace, K)
        u.vector()[:] = project(U, Vspace).vector()
        plot(u, title="Velocity")

	# Compute wall shear stress
	wss_ = WSS(U) 
        
	# We project the wall shear stress to have it in the 
	# same space as the geometry space
	wss = project(wss_, CG1)
	wss = interpolate(wss, Kspace)
	#plot(wss, title = 'wss')	

	# Cmpute boundary	
	bdry_ = bdry(K, DG0)
	#plot(bdry_, title = 'bdry_')

	# Compute indicator function
	ind_f = ind_func(bdry_, wss) #, interesting_domain_proj)
	#plot(ind_f, title='Indicator function')	

	# Update K
	update_K(ind_f, K)
	
	# Cap K
	capping(K)
	plot(K, title='K')

	yield K

plot(K, title="K")
for K in iterative(K, 50):
    pass
