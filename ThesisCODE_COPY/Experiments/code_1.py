
import numpy
from matplotlib import pyplot
from dolfin import * 

# LEVEL SET FUNCTIONS, WSS, INDICATOR
from stokes_solver import stokes_solver
from skjaer_vegg_indikator import tau_xy, interface, indicator   
from K_update_kappet import K_update, K_cap
from skjaer_vegg_indikator_backup import indikator_test



# Parameters 
N = 64              # the mesh size              
mesh = RectangleMesh(Point(-1.0,-1.0) ,Point( 1.0, 1.0), N, N, "right")

x0 = DOLFIN_EPS     # the center point of the vessel 0.0 
R = 0.10            # the radius of the vessel 0.20 # 0.10 
viz = True          # vizualization parameter
speed = 1.0         # factor multiplied with the Poiseuille flow to alternatie the speed  
h = 1./N            # the mesh size  
x1 = 0.2            # the area to be cut at the inlet 
mu = 1              # the dynamic viscosity 
dpdy = 2.1          # the pressure gradient  
epsilon = 2         # the growth or shrinking factor 
alg = 'g'           # algorithm 'o': old vs 'n' : new


time_start = 0     # time start
time_stop = 1      # time stop




# Functionspaces 
V = VectorFunctionSpace(mesh, 'Lagrange', 2) 
P = FunctionSpace(mesh, 'Lagrange', 1)
W = MixedFunctionSpace([V, P]) 
DG1 = FunctionSpace(mesh, 'DG',1)






# The surface
phi_I = Expression('( x[0] > (x0 - R - DOLFIN_EPS) && x[0] < (x0 + R + DOLFIN_EPS)) ? 0.0 : 1.0', x0 = x0, R = R ) 
K  = interpolate(phi_I, P) # alternativt DG0 
plot(K, interactive = viz, title = 'K' )	
#File('K_c105_n{}_s{}_e{}_alg{}.xdmf'.format(N, speed, epsilon, alg)) << K




#the interesting domain 
interesting_domain = Expression('(x[0] > - 1 + x1   &  x[1] > -1 + x1   &  x[0] < 1 -x1   & x[1] < 1 - x1) ? 1.0:0.0 ', x1 = x1) 
interesting_domain_interpolate = interpolate(interesting_domain, DG1)#P)          
#plot(interesting_domain_interpolate,interactive = True, title = 'interesting domain ')
# CG1 viser seg aa vaere feil valg av funksjonsrom. Selv med vokst mapping saa fungerer det ikke. 







# Plotting the shear with high resolution
def dg1_line_plot(function, resolution=1000):
	x = numpy.linspace(-1, 1, resolution)      # Needs to be changed if the geometry changes
	y = [function((xval, 0.0)) for xval in x]  # Skjaret over domenet som regnet ut i modellen 
	pyplot.plot(x, y, linestyle="--")
	pyplot.show()









# Iterating over the surface 
def time_stepping(K, time):

	print 'time level: ', time 
	

	# Calculating the velocity and the pressure  		
	u, p = stokes_solver(x0=x0, R=R, K=K, mesh=mesh, W=W, speed=speed, mu=mu, dpdy=dpdy)
	#plot(u, interactive = True, title = 'Velocity field')
	#File('u_c105_n{}_t{}_s{}_e{}_alg{}.xdmf'.format( N, time, speed, epsilon, alg)) << u




	# Calculating the shear 
	tau_xy_ = tau_xy(u=u)
	tau_xy_proj = project(tau_xy_, DG1) #P)  
	#plot(tau_xy_proj, interactive = True, title = 'Shear')	
	dg1_line_plot(tau_xy_proj)	
	File('code1_shear_c105_n{}_t{}_s{}_e{}_alg{}.xdmf'.format(N, time, speed, epsilon, alg)) << tau_xy_proj
	



	# Calculating the interface 
	interface_ = interface(K=K) 
	interface_proj = project(interface_, DG1) #P
	#plot(interface_proj, interactive = viz, title = 'Interface')
	#File('interface_c105_n{}_t{}_s{}_e{}_alg{}.xdmf'.format(N, time, speed, epsilon, alg)) << interface_proj




	# Calculating the indicator
	#indicator_ = indicator(space = DG1, interface=interface_proj, tau_xy = tau_xy_proj, interesting_domain = interesting_domain_interpolate, epsilon = epsilon)

	indicator_ = indikator_test(space = DG1, Dirac_delta =interface_proj, tau_xy=tau_xy_proj, interesting_domain= interesting_domain_interpolate)
	plot(indicator_, interactive = viz, title = 'Indicator')
	#File('Indicator_c105_n{}_t{}_s{}__e{}_alg{}.xdmf'.format(N, time, speed,epsilon,alg)) << interface_proj 



	# Updating K : K_new = K_old - I  
	K_new = K_update(indikator = indicator_, K=K)
	#plot(K_new, interactive = viz, title = 'K_ny')
	#File('Knew_c105_n{}_t{}_s{}_e{}_alg{}.xdmf'.format(N, time, speed, epsilon, alg)) << K_new  




	K_cap_ = K_cap(K_new)
	K_cap_proj = project(K_cap_, K.function_space())
	plot(K_cap_proj, interactive = viz, title = 'K cap' )
	#File('Kcap_c105_n{}_t{}_s{}_e{}_alg{}.xdmf'.format(N, time, speed, epsilon, alg)) << K_cap_proj 

	return K_cap_proj  
		 









# The time loop 

for t in range(time_start, time_stop):
	K_cap_proj = time_stepping(K, t)
	K = K_cap_proj 






