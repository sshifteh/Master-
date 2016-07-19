

from dolfin import * 
from numpy import where 
# LEVEL SET FUNCTIONS, WSS, INDICATOR
from stokes_solver_v2 import stokes_solver
from tau_xy import tau_xy, Dirac_delta, wall_shear_stress, indikator_test1, indikator_test2   
from K_update import K_update, K_kappet, K_update_v2 





N =32
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







#interesting domain 
x1 = 0.1 
interesting_domain = Expression('(x[0] > x1 & x[1] > x1 & x[0] < 1- x1 & x[1] < 1-x1) ? 1.0 : 0.0', x1 =x1)
interesting_domain_interpolate = interpolate(interesting_domain, dg1) 
plot(interesting_domain_interpolate,interactive = True, title = 'interesting domain ')


def iterative(K):

	print 'iteration', i



	u, p = stokes_solver(x0=x0, R=R, K=K, mesh=mesh, W=W)
	plot(u,interactive = viz, title = 'veloctiy , u')
	plot(p, interactive = viz, title = 'pressure, p')


	
	# regner ut skjaeret 
	tau_xy_ = tau_xy(u=u)
	#File("tau_xy.xdmf") <<  tau_xy_ 	
	#plot(tau_xy_, interactive = viz, title = 'tau_xy ikke projec')	
	tau_xy_proj = project(tau_xy_, dg1) # naturlig cg1 tror jeg  
	plot(tau_xy_proj, interactive = viz, title = 'shear stress, tau_xy_proj')
	File('tau_xy_proj.xdmf') << tau_xy_proj


	# regner ut interfacet 
	Dirac_delta_ = Dirac_delta(K=K) 	
	Dirac_delta_proj = project(Dirac_delta_, dg1) # naturlig ville vaert dg0 tror jeg 		
	plot(Dirac_delta_proj, interactive = viz, title = 'interfacet ' )
	#File("Dirac_delta_proj.xdmf") << Dirac_delta_proj  
	


	
	wall_shear_stress(tau_xy = tau_xy_proj, Dirac_delta= Dirac_delta_proj)
	



	indikator_funksjon1 = indikator_test1(space = dg1, Dirac_delta=Dirac_delta_proj, tau_xy = tau_xy_proj)
	indikator_funksjon2 = indikator_test2(space = dg1, Dirac_delta=Dirac_delta_proj, tau_xy = tau_xy_proj, interesting_domain = interesting_domain_interpolate)
	plot(indikator_funksjon2, interactive = viz, title = 'indikator')
	#FIXME problemet er indikatoren

	#indikator_funksjon = indikator_test1(space = dg1, Dirac_delta=Dirac_delta_proj, tau_xy = tau_xy_proj)
	#plot(indikator_funksjon, interactive = viz, title = 'indikator')
	

	# for aa oppdatere K saa maa indikator funksjonen interpoleres til rommet av cg1 funcsjoner 
	# for naa er det i samme som som skjaeret of dirac delta funksjonen, nemlig dg1 		
	
	# vi projiserer det paa funksjons rommet til K i funksjonen men mulig det er bedre aa interpolere egentlig. mindre comp effort. 	
	#indicator_function_interpolate = interpolate(indicator_function, P)
	#plot(indicator_function_interpolate, interactive = True, title = 'indicator funksjnen interpolert til rommet av cont lin funksjoner ')	
	# fungerer helt fint 

	#K_update_v2(indikator = indikator_funksjon, K=K)

	#K_new = K_update_v2(indikator = indikator_funksjon2, K=K)
	K_new = K_update(indikator = indikator_funksjon2, K=K)
	plot(K_new, interactive = viz, title = 'K_ny')
	
	K_kappet_ = K_kappet(K_new)
	K_kappet_proj = project(K_kappet_, K.function_space())
	plot(K_kappet_proj, interactive = viz, title = 'K kappet' )
		
	#stokes_solver(x0=x0, R=R, K= K_kappet_proj, mesh=mesh, W=W) # helt til slutt regner vi ut stokes over det nye domenet
	
	return K_kappet_proj  
	


for i in range(3):
	K_kappet_proj = iterative(K)
	K = K_kappet_proj 
	


