from alpha import alpha
from dolfin import *
from boundary import boundary
from numpy import where 	

# Faster C++ version    
u_analytical = '''  
class AnalyticalVelocity : public Expression
{
	// Parameters 
	public:
		double x0, R, nu, dpdz;

	// The eval function for the function expression 
        void eval(Array<double>& values, const Array<double>& x) const
	{
		if(((x0 - R - DOLFIN_EPS) <= x[0])   &&  (x[0] <= (x0 + R + DOLFIN_EPS)))
		{  
			double velocity = 0.25*nu*dpdz*((x[0]-x0)*(x[0]-x0) - R*R); 		
			values[0] = 0.0;
                        values[1] = velocity; 
		}
		else 
		{
			values[0] = 0.0; 
                        values[1] = 0.0;
		}
	}
	// Making the array 2 dimensional 
        size_t value_rank() const{
            return 1;
        }

        size_t value_dimension(size_t i) const{
            return 2;
        }

};''' 


# Experimental solver
def stokessolver(K, nu, dpdz, x0, R):  

	# Test and trial spaces over the Taylor-Hood elements  
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
	

	
	# Source term
	f = Constant([0.0,0.0])
	
		
	
	# Variational form on the geometry above 
 	a =  inner(alpha(u, K), v)*dx + nu*inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)
	
	
	
	# Parabolic inlet velocity Dirichlet condition 
	analyticalVelocity = Expression(["0","0.25*nu*dpdz*((x[0]-0.50)*(x[0]-0.50) - R*R) "], nu=nu, R=R, dpdz=dpdz)
	# Defining the boundaries of the vessel 
	def T_inlet_u(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > (x0 - R - DOLFIN_EPS) and x[0] < (x0 + R + DOLFIN_EPS) and on_bnd


	def T_notInlet(x, on_bnd):
		return x[1] < DOLFIN_EPS and  not (x[0]> (x0 - R - DOLFIN_EPS) ) and not (x[0] < (x0 + R + DOLFIN_EPS) ) and on_bnd


	def T_notOutlet(x, on_bnd):
		return x[1] < 1 - DOLFIN_EPS and not (x[0] > (x0 - R - DOLFIN_EPS)) and not (x[0] < (x0 + R + DOLFIN_EPS)) and on_bnd
	
	
    	def T_leftWall(x, on_boundary):
        	return near(x[0], 0.0)

	def T_rightWall(x, on_boundary):
        	return near(x[0], 1.0)



	inlet  = DirichletBC(Wspace.sub(0), analyticalVelocity, T_inlet_u)
	notInlet = DirichletBC(Wspace.sub(0), analyticalVelocity, T_notInlet)
	notOutlet = DirichletBC(Wspace.sub(0), analyticalVelocity, T_notOutlet)
	right = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), T_rightWall)
	left = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), T_leftWall)
	
	bcs = [inlet,notInlet, notOutlet, right, left] 



	# Solving the system of linear equations 
	A, b = assemble_system(a, L, bcs)
        solver = LUSolver(A)
	solver.solve(UP.vector(), b) 
	U, P = UP.split()
       
	
	File("U.xdmf") << U 	
	File('P.xdmf') << P
	File('u_analytical.xdmf') << interpolate(u_c, fine_V)
	

	# Plotting the experimental solution on a superfine mesh with 128 cells 
	plot(interpolate(U, fine_V), mesh, title = 'Experimental Velocity Field, N')
	plot(interpolate(P, fine_P), mesh, title = 'Experimental Pressure')	
	# Plotting the analytical solution on a fine mesh 
	plot(interpolate(u_c, fine_V), mesh, title = 'Analytical Velocity Field, N')
	interactive()

	return U,P
	


def shear(U):
	shear = 0.5*(U[0].dx(1) + U[1].dx(0)) 
	return shear




def tau_wall(bdry, Shear):
	
	shear_ = Shear.vector().array()
        bdry = interpolate(bdry, Shear.function_space())
	bdry_ = bdry.vector().array()

        assert len(shear_) == len(bdry_)

        bdry_idx = where(bdry_ > DOLFIN_EPS)[0]
        # Plotting wss_bdry is a pain, but iterpolation to DG0 seems to work
        # okish
        #DG0 = FunctionSpace(WSS.function_space().mesh(), "DG", 0)
        #plot(interpolate(wss_bdry, DG0), interactive=True)
        wss = shear_[bdry_idx]
	print wss 
	#return wss


def wss_analytical(dpdz, nu, R):
	return abs(-0.5*nu*(R-0.5))


if __name__ == '__main__':

	# Physical parameters 
	center = 0.50
	side_width = 0.1
	nu = 0.00345
	dpdz = -1 #-2.1
			
	u_c = Expression(u_analytical) 
	u_c.x0 = float(center);
	u_c.R  = float(side_width);
	u_c.nu = float(nu);
	u_c.dpdz = float(dpdz);




	# Create mesh and functionspaces for approximation of functions 
	N = 128
	mesh   = UnitSquareMesh(N,N)
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
	Kspace = FunctionSpace(mesh, 'DG', 0)
	CG1 = Pspace 
	DG0 = Kspace 	
	
	# Create a fine mesh for comparing analytical and numerical solution	
	fine_mesh = UnitSquareMesh(128,128)
	fine_V = VectorFunctionSpace(fine_mesh, 'Lagrange', 2)
	fine_P = FunctionSpace(fine_mesh, 'Lagrange',1)		


	# Introduse vessel geometry 
	K_expr = Expression('( x[0] > (x0 - R - DOLFIN_EPS) && x[0] < (x0 + R + DOLFIN_EPS)) ? 0.0 : 1.0', x0 = center, R = side_width )
	# Interpolate it into the Kspace
	K = interpolate(K_expr, Kspace)	
	
	# Program flow 
	U, P = stokessolver(K = K, nu = nu, dpdz = dpdz, x0 = center, R = side_width) 
	shear = shear(U)
	shear_proj = project(shear, CG1)
	shear_interpolated = interpolate(shear_proj, Kspace)
	plot(shear_interpolated, interactive = True) 		
	bdry_ = boundary(K, DG0) 
	
	wss_vector = tau_wall(bdry_, shear_interpolated)	
	wss_function = Function(shear_interpolated.function_space())
	wss_function.vector()[:] = wss_vector[:]		
	plot(wss_function, interactive=True)	



	# Comments 
	# wss_analytical maa gjoeres til like stor vektor som wss eksperimental og smamenliknes
	# og plottes 

