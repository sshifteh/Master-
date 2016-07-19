
# Import statements 
# Importing the UFL forms, mesh and Function Spaces from the PDE equation module dolfin 
from dolfin import *

# Importing the 6 helper functions 
from K_update import update_K
from capping import capping
from stokes_solver import stokes_solver
from WSS import WSS
from boundary import boundary
from indicator_function import indicator_function



# For plotting in succession like a movie  

show_plot = True 
def plot(*args, **kwargs): # *args unknown how many args the function takes
    if not show_plot:      # **kwarg for named arguments not defined in advance 
        return		   # if not True, then False, nothing will happen 		
    dolfin.plot(*args, **kwargs) # Plot will be activated with fexible number of arguments in all plot statements 




# Parameters and mesh 

n = 40
w = 0.25
mesh = UnitSquareMesh(n,n) #'crossed')

# Function spaces for the approximations to the functions  
CG1 =FunctionSpace(mesh, 'CG', 1)
DG1 =FunctionSpace(mesh, 'DG', 1)
DG0 = FunctionSpace(mesh, 'DG', 0)

Pspace = CG1
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace])

Kspace = DG0
WSSspace = DG1

# The vessel geometry and interpolation into the Kspace 
K_expr = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)
K = interpolate(K_expr, Kspace)

# The approximation of the velocity 
u = Function(Vspace)



# Defining an iterative function that will show the dynamics of the domain K 
def iterative(K, iters):

    for i in range(iters):
	# Solve for velocity and pressure
	U, P = stokes_solver(w, mesh, Vspace, Pspace, Wspace, K)
        u.vector()[:] = project(U, Vspace).vector() 
        plot(u, title="Velocity")

	# Compute wall shear stress
	wss_ = WSS(U) 
        
	# We project the shear stress to have it in the 
	# same space as the geometry space
	wss = project(wss_, CG1)
	wss = interpolate(wss, Kspace)
	#plot(wss, title = 'wss')	

	# Compute boundary	
	bdry_ = boundary(K, DG0)
	#plot(bdry_, title = 'bdry_')

	# Compute indicator function
	ind_f = indicator_function(bdry_, wss) #, interesting_domain_proj)
	#plot(ind_f, title='Indicator function')	

	# Update K
	update_K(ind_f, K)
	
	# Cap K
	capping(K)
	plot(K, title='K')

	# yield is a new style return statement with doesnt save old frames 
	yield K

plot(K, title="K")
for K in iterative(K, 50):
    pass
