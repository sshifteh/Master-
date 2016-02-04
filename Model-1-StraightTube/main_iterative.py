
from dolfin import * 
#from indicator_function_expression import ind_func_expr
from K_update import K_update
from capping import capping
from stokes_solver import stokes_solver 
from WSS import WSS
from WSS_bdry import WSS_bdry, bdry  
from indicator_funciton_conditions2 import ind_func


n = 50
w= 0.45
mesh = UnitSquareMesh(n,n) #, 'crossed')
Pspace = FunctionSpace(mesh, 'CG',1)
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace]) 
DG1 =FunctionSpace(mesh, 'DG', 1) 	# ShearStressSpace
DG0 = FunctionSpace(mesh, 'DG', 0)	# DPspace CHANGE this name to DG0 much more intuitive than discont. pressure space


K_1  = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)   


# interesting domain 
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
Top().mark(subdomains, 0)
#Right().mark(subdomains, 0)
#Left().mark(subdomains, 0)


interesting_domain = SubDomainIndicator([Bottom(), Top()])
plot(interesting_domain, mesh = mesh, interactive = True)
interesting_domain_proj = project(interesting_domain, DG1)

#interesting_domain = Expression(" (x[0] > 2.0/n && x[0] < 48.0/n) || (x[1] >2.0/n && x[1] < 48.0/n ) ? 1.0: 0.0", n=n)
#interesting_domain_proj = project(interesting_domain, DG1) 
#plot(interesting_domain_proj, interactive = True, title= 'intersting_domain_proj')

def iterative(K_1):

	U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, DG1=Pspace, K_array=K_1, n=n )


	# Verification: 

	# alternative 1:
	# construct an analytical solution 
	# method of manufactored solutions 
	# first come up with a random slution
	# plut into PDE 
	# doesnt solve original homogeneous PDE exactly 
	# residual is then considered the source term 
	
		
	# alternative 2:
	# really fine mesh , solve a blind model on that mesh 
	# this is reference(considered exact) solution
	
		
		
	# actual convergence test(after alt 1 or 2):
	# mesh refinement test refining the mesh see if rate of convergence goes to 2.
	# mesh.refine() OR refine(mesh)
	# give a mesh, creates new mesh that is finer. 
	# for each mesh solve PDE
	# compute difference of u_e - u
	# in a norm 
	# for each mesh 

	# start with a very coarse one 
	# after halfin 4 times,gives a really fine mesh and solving will get slow. 

	# example:
	# for p2, half the mesh size 
	# error reduction of factor 4 

	# do it at the beginning

	# Magne: 
	# WSS skal naturlig veare i DG1.
	# og det fungerte aa bytte Pspace m DG1 
	
	wss = WSS(U=U, DG1= DG1) #Pspace) #,interesting_domain = interesting_domain_proj)
	bdry_ = bdry(K_Func = K_Func, DG1 = DG1)

	# Simon: 
	# DG1 osgaa inneholder CG1 


	ind_f = ind_func(bdry=bdry_, WSS = project(wss, DG1), DG1 = DG0, interesting_domain = interesting_domain_proj)
	K_new = K_update(ind_func = ind_f, K_Func = K_Func)
	
	#wss_bdry = WSS_bdry(K_Func, DG0, WSS)


	plot(K_new, interactive = True, title = 'K new ')
	K_capped = capping(K_new)
	plot(K_capped, interactive = False, title = 'K capped')
	stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_capped, n=n )

	return K_capped 
	

iterative(K_1=K_1)

for i in range(5):
	K_updated = iterative(K_1 = K_1)
	K_1= K_updated # this didnt work 
	


