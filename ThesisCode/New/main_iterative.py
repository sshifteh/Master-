
import numpy
from matplotlib import pyplot
from dolfin import * 
# LEVEL SET FUNCTIONS, WSS, INDICATOR
from stokes_solver import stokes_solver
from skjaer_vegg_indikator import tau_xy, interface, indicator   
from K_update_kappet import K_update, K_cap
#from skjaer_vegg_indikator_backup import indikator_test



# Parameters 
N = 32
speed = 20        

             
mesh = RectangleMesh(Point(-1.0,-1.0) ,Point( 1.0, 1.0), N, N, "right")
x0 = DOLFIN_EPS     
R = 0.20            
viz = True          
h = 1./N            
x1 = 0.2            
mu = 1              
dpdy = 2.1          
epsilon = 2         
alg = 'o'           


time_start = 0     # time start
time_stop = 7      # time stop


# Functionspaces 
V = VectorFunctionSpace(mesh, 'Lagrange', 2) 
P = FunctionSpace(mesh, 'Lagrange', 1)
W = MixedFunctionSpace([V, P]) 
DG1 = FunctionSpace(mesh, 'DG',1)


# The surface
phi_I = Expression('( x[0] > (x0 - R - DOLFIN_EPS) && x[0] < (x0 + R + DOLFIN_EPS)) ? 0.0 : 1.0', x0 = x0, R = R )
K  = interpolate(phi_I, P) # alternativt DG0 



#the interesting domain 
interesting_domain = Expression('(x[0] > - 1 + x1   &  x[1] > -1 + x1   &  x[0] < 1 -x1   & x[1] < 1 - x1) ? 1.0:0.0 ', x1 = x1) 
interesting_domain_interpolate = interpolate(interesting_domain, DG1)#P)          





# Iterating over the surface 
def time_stepping(K, time):

	print 'time level: ', time 
	

	# Calculating the velocity and the pressure  		
	u, p = stokes_solver(x0=x0, R=R, K=K, mesh=mesh, W=W, speed=speed, mu=mu, dpdy=dpdy)
	

	# Calculating the shear 
	tau_xy_ = tau_xy(u=u)
	tau_xy_proj = project(tau_xy_, DG1) #P)  
	#File('shear_N{}_t{}_s{}.xdmf'.format(N, time, speed)) << tau_xy_proj
	

	# Calculating the interface 
	interface_ = interface(K=K) 
	interface_proj = project(interface_, DG1) #P
	
	# Calculating the indicator
	#indicator_ = indicator(space = DG1, interface=interface_proj, tau_xy = tau_xy_proj, interesting_domain = interesting_domain_interpolate, epsilon = epsilon)

	indicator_ = indicator(space = DG1, interface =interface_proj, tau_xy=tau_xy_proj, interesting_domain= interesting_domain_interpolate, epsilon = epsilon)


	# Updating K : K_new = K_old - I  
	K_new = K_update(indikator = indicator_, K=K)
	

	K_cap_ = K_cap(K_new)
	K_cap_proj = project(K_cap_, K.function_space())
	plot(K_cap_proj, interactive = viz, title = 'K cap' )
	#File('Kcap_n{}_t{}_s{}_e{}_alg{}.xdmf'.format(N, time, speed, epsilon, alg)) << K_cap_proj 

	return K_cap_proj  
		 


for t in range(time_start, time_stop):
	K_cap_proj = time_stepping(K, t)
	K = K_cap_proj 

















 
