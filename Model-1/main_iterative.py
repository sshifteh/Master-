
# DOLFIN LIBRARY AND PLOTTING LIBRARY 
from dolfin import * 
#from plotting import plot
#parameters['plotting_backend'] = "matplotlib"
#from matplotlib import pyplot

# LEVEL SET FUNCTIONS, WSS, INDICATOR
from K_update import K_update
from capping import capping
from stokes_solver import stokes_solver 
from WSS import WSS
from WSS_bdry import WSS_bdry  
from indicator_funciton_conditions2 import ind_func_expr


# FINITE ELEMENT MESH, FUNCTION SPACES
n = 50
# refine mesh dolfin and the arg that was none set ut to epsilon around the Y shape 
# also it can be done in parallel p 217 distributed mem paralel computing
# mpirun -n 16 16 no of processes 
# use y on coarse mesh, one iteration, check serial and paralel produce the same result
# plots 2 meshes cuz -- > expect 2 windows
# alternative save result to file , open up in paraview will not be split will give the full mesh 


# start with a coarse mesh on the non interesting areas and 
# refine 3 x near the interesting Y area 

#w=0.45 # for straigh tube 
w= 0.01
mesh = UnitSquareMesh(n,n, 'crossed')
mesh = refine(mesh)
#plot(mesh)
Pspace = FunctionSpace(mesh, 'CG',1)
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace]) 
DG1 =FunctionSpace(mesh, 'DG', 1)  # ShearStressSpace
DG0 = FunctionSpace(mesh, 'DG', 0) # DPspace CHANGE this name to DG0. Pressure space


# THE FOUR DIFFERENT LEVEL SET FUNCTIONS 

# T - SHAPED TUBE
phi_T = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1- w + DOLFIN_EPS && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 )? 0.0:1.0', w=w)

# I - SHAPED TUBE 
phi_I  = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)   


# L - SHAPED TUBE
phi_L = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 && x[0] > w )? 0.0:1.0', w=w)	
	
# Y - SHAPED TUBE
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


# GENARAL CLASS FOR MAKING AN EXPRESSION OF ITEMS IN LIST 
class SubDomainIndicator(Expression):
    def __init__(self, subdomain_list):
		self.subdomain_list = subdomain_list

    def eval(self, values, x):
        values[0] = 1.
        for subdomain in self.subdomain_list:
            if subdomain.inside(x, 0):
                values[0] = 0.
                break

subdomains = CellFunction('size_t', mesh, 1)
Stem().mark(subdomains, 0)
LeftBranch().mark(subdomains, 0)
RightBranch().mark(subdomains, 0)

phi_Y = SubDomainIndicator([Stem(), LeftBranch(), RightBranch()])
#plot(phi_Y, mesh = mesh, interactive = True)


def iterative(phi_Y):
	U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=phi_Y, n=n )

	#plot(K_Func,interactive = True, title = 'Geometry')				
	
	wss = WSS(U=U, DG1=DG1)
	wss_bdry = WSS_bdry(K_Func=K_Func, DPspace=DG0, WSS=wss)
 	ind_f = ind_func_expr(WSS_bdry = wss_bdry)
	K_new = K_update(ind_func = ind_f, K_Func = K_Func)
	
	#K_fig = plot(K_new, interactive = True, title = 'K new ')
	#K_fig.save	
	
	K_capped = capping(K_new)
	#plot(K_capped, interactive = False, title = 'K capped')
	stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_capped, n=n )

	return K_capped  
	


for i in range(6):
	K_updated = iterative(phi_Y = phi_Y)
	

	
	plot(K_updated,interactive = True, title = 'K')
	
	#plot(K_updated)	
	#pyplot.suptitle("Y-shape:K updated: %g t" %i)
	#pyplot.savefig("Y-shape:K updated: %g t" %i)
	#pyplot.figure()	
		
	
	#pyplot.show()

	phi_I = K_updated # this didnt work 
	


