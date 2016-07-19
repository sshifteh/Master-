

from dolfin import * 
from numpy import where 
# LEVEL SET FUNCTIONS, WSS, INDICATOR
from stokes_solver import stokes_solver
from tau_xy import tau_xy, Dirac_delta, wall_shear_stress, indikator_test1, indikator_test2   
from K_update import K_update, K_kappet, K_update_v2 





N = 64
x0 = 0.50 
R = 0.10 
viz = True
mesh = UnitSquareMesh(N,N)
#mesh = refine(mesh)
h = 1./N

V = VectorFunctionSpace(mesh, 'Lagrange', 2) 
P = FunctionSpace(mesh, 'Lagrange', 1)
W = MixedFunctionSpace([V, P]) 
dg0 = FunctionSpace(mesh, 'DG',0)
dg1 = FunctionSpace(mesh, 'DG',1)

phi_I = Expression('( x[0] > (x0 - R - DOLFIN_EPS) && x[0] < (x0 + R + DOLFIN_EPS)) ? 0.0 : 1.0', x0 = x0, R = R )
K  = interpolate(phi_I, P) # alternativt DG0 
plot(K, interactive = viz, title = 'K' )	
#File('K_iteration.xdmf') << K






#interesting domain 
x1 = 0.1 
interesting_domain = Expression('(x[0] > x1 & x[1] > x1 & x[0] < 1- x1 & x[1] < 1-x1) ? 1.0 : 0.0', x1 =x1)
interesting_domain_interpolate = interpolate(interesting_domain, dg1) 
plot(interesting_domain_interpolate,interactive = True, title = 'interesting domain ')


def iterative(K, i):

	print 'iteration', i



	u, p = stokes_solver(x0=x0, R=R, K=K, mesh=mesh, W=W)
	#plot(u,interactive = viz, title = 'veloctiy , u')
	#plot(p, interactive = viz, title = 'pressure, p')
	# FIXME rename 	
	#File('u_speed2_iteration={}.xdmf'.format(i)) << u
	#File('p_speed0.5_iteration={}.xdmf'.format(i)) << p


	
	# regner ut skjaeret 
	tau_xy_ = tau_xy(u=u)
	#plot(tau_xy_, interactive = viz, title = 'tau_xy ikke projec')	
	tau_xy_proj = project(tau_xy_, dg1) # kontonuerlig over hvert element, men diskontinuerlig fra element til element  
	plot(tau_xy_proj, interactive = viz, title = 'shear stress, tau_xy_proj')
	# FIXME rename	
	#File('skjaeret_proj_DG1_speed0.5_iteration={}.xdmf'.format(i)) << tau_xy_proj


	# regner ut interfacet 
	Dirac_delta_ = Dirac_delta(K=K) 	
	Dirac_delta_proj = project(Dirac_delta_, dg1) # naturlig ville vaert dg0 tror jeg 		
	#plot(Dirac_delta_proj, interactive = viz, title = 'interfacet ' )
	# FIXME rename
	#File("veggen_proj_DG1_speed2_iteration={}.xdmf".format(i)) << Dirac_delta_proj  
	


	##wall_shear_stress(tau_xy = tau_xy_proj, Dirac_delta= Dirac_delta_proj)


	#indikator_funksjon1 = indikator_test1(space = dg1, Dirac_delta=Dirac_delta_proj, tau_xy = tau_xy_proj)
	indikator_funksjon2 = indikator_test2(space = dg1, Dirac_delta=Dirac_delta_proj, tau_xy = tau_xy_proj, interesting_domain = interesting_domain_interpolate)
	plot(indikator_funksjon2, interactive = viz, title = 'indikator fra main i dg1')
	# FIXME rename	
	# File('indikatoren_speed2_iteratsjon={}.xdmf'.format(i)) << indikator_funksjon2

	#indikator_funksjon = indikator_test1(space = dg1, Dirac_delta=Dirac_delta_proj, tau_xy = tau_xy_proj)
	#plot(indikator_funksjon, interactive = viz, title = 'indikator')
	

	# for aa oppdatere K saa maa indikator funksjonen interpoleres til rommet av cg1 funcsjoner 
	# for naa er det i samme som som skjaeret of dirac delta funksjonen, nemlig dg1 		
	
	# vi projiserer det paa funksjons rommet til K i funksjonen men mulig det er bedre aa interpolere egentlig. mindre comp effort. 	
	#indicator_function_interpolate = interpolate(indicator_function, P)
	#plot(indicator_function_interpolate, interactive = True, title = 'indicator funksjnen interpolert til rommet av cont lin funksjoner ')	
	# fungerer helt fint 


	K_ny = K_update(indikator = indikator_funksjon2, K=K)
	plot(K_ny, interactive = viz, title = 'K_ny')
	# FIXME rename	
	# File('K_ny_speed2_iterasjon={}.xdmf'.format(i)) << K_ny	


	K_kappet_ = K_kappet(K_ny)
	K_kappet_proj = project(K_kappet_, K.function_space())
	plot(K_kappet_proj, interactive = viz, title = 'K kappet' )
	# FIXME rename	
	# File('K_kappet_proj_cg1_speed2_iterasjon={}.xdmf'.format(i)) << K_kappet_proj		
	
	return K_kappet_proj  
	#return K_new2	


for i in range(5):
	K_kappet_proj = iterative(K, i)
	K = K_kappet_proj 
	


