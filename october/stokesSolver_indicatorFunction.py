
def alpha(u, k):
	C = Constant(1e5) # high value outside the vessel
	return C*k*u


def stokes_solver(w, func, n = 50 ):
	mesh = UnitSquareMesh(n,n)
	
	Tspace = TensorFunctionSpace(mesh, "DG", 1)	
	ShearStressSpace = FunctionSpace(mesh, "DG", 1)	
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
	
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
		
	k   = interpolate(func, Pspace) # control                     # Function object 
	plot(k, interactive = True, title = 'tube with CG elements')
	
	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, k), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)

	# boundary conditions____________________________________________________________________________
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > w and x[0] < 1-w and on_bnd
	 	
	def p_boundaryRightWall(x, on_bnd):
		return x[0] > 1- DOLFIN_EPS and x[1] >= 0.5 and x[1] < 0.7 and on_bnd # why is it read on the entire unitsquare boundary?S 
	
	bc_u_B  = DirichletBC(Wspace.sub(0), Constant((0.0,1.0)), u_boundaryBottom)
	
	velocityFunc = Expression(["0","x[0]-x[0]*x[0]"])
	bc_u_B  = DirichletBC(Wspace.sub(0), velocityFunc, u_boundaryBottom)
	bc_p_R  = DirichletBC(Wspace.sub(1), Constant(0.0), p_boundaryRightWall)
	bcs = [bc_u_B, bc_p_R]


	# Solve system of linear equations_______________________________________________________________
	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
		
	U, P = UP.split()
	plot(U, interactive = True, title= "Velocity")
	plot(P, interactive = True, title= "Pressure")

	# Handling of the WSS ___________________________________________________________________________

	shear_stress = project(0.5*(U[0].dx(1) + U[1].dx(0)), ShearStressSpace) # calculate the WSS over the functionspace
	plot(shear_stress, interactive = True, title = 'WSS')	                # plot it 
	WSSarray = shear_stress.vector().array()                                # put all the calculated values in an array
        # do some computations with array                                       #just print out the values to see what they are
	print WSSarray
	for wss in WSSarray:
		print wss 
	print len(WSSarray) 							# 15 000, it means WSS is calculated in 15 000 DC nodes

	indicator_function = Function(ShearStressSpace) 
	plot(indicator_function, interactive = True, title = ' indicator function') # nothing yet  
	
	indicator_func_array = indicator_function.vector().array() 
	print type(indicator_func_array) # okei its an array with zeros. It should have same length as WSSarray 
	print indicator_func_array
	print 'size of indicator_array', len(indicator_func_array) # It is of size 15 000 
	indicator_func_array[:] = WSSarray                         # However, we make sure of it
	len(indicator_func_array)
	print indicator_func_array                        # Also the values of WSSarray have been copied into it. Is this desirable?

	

	# Now that WSSarray and indicator_func_array are two arrays of the same length we can continue the next step which is 
	# define d1 and d2 with random values, based on the numbers I can see, these will be choosen when the model is optimized later
	# lets consider the smalles value and the largest first before deciding on d1 and d2	
	# just to get a feel of things 
	
	smallest = min(abs(WSSarray))
	largest  = max(abs(WSSarray))
	print 'smallest WSS value',smallest, 'largest WSS value', largest 
	# smallest WSS value 4.07399371391e-08 largest WWS value 21.7227328423
	# The smallest WSS value the vessel wall is exposed to is 4.1x10-8, and the largest is 21.7. 
	# This is a big spectrum. Lets take the log and consider it. 
 
	# smallest_log = log(smallest)
	# largest_log  = log(largest) 
	# print 'log smallest and log largest', smallest_log, largest_log
	# it will not work cuz log means data log in fenics
	# with ln the smallest becomes -17 and largest 3. 
	# this might be easier to handle, no not really because we are taking abs values so this is useless.

	# Shrinking for WSS less than some value crazy random guess : Let d1 be  
	
	"""
	for wss in WSSarray:
		print wss
		if  20<abs(wss)<22 : 
			print 'found it!'
			
	
	for i in range(len(WSSarray)):          # its has index 205 in array 
		if  20 < abs(WSSarray[i])<22 : 
			print i, WSSarray[i]
			break
	"""
	
	# gamma is the indicator function	
	d1 = 1; d2 = 10 ; eps = 1e-5		
	
	print WSSarray[205]
	
	for i in range(len(WSSarray)):
		if abs(WSSarray[i]) < d1:
			indicator_func_array[i] = -eps

		elif  d1 < abs(WSSarray[i]) < d2:
			indicator_func_array[i] = 0
		
		elif abs(WSSarray[i]) > d2:
			indicator_func_array[i] = eps
		

	# Now we have to update k (Capital K in printouts)
	k   = interpolate(func, ShearStressSpace)  
	k_array = k.vector().array()
	
	#for k in k_array: here it is discrete, only taking values 0 or 1 
	#	print k  	
	#print len(k_array) #15 000
	


	#initalize first the new k_updated array
	k_updated = k_array[:]
	#print k_updated
	
	#print 'k array before', k_array
	for i in range(len(k_array)):
		k_updated[i] = k_array[i] - indicator_func_array[i]
	#print 'k array after', k_updated

	
	#for k in k_updated:
	#	print k
	# they have indeed changed	
	
	k_updated_projected = interpolate(k_updated, ShearStressSpace) 
	plot(k_updated_projected, interactive = True, title = 'k_updated')

	# make this a recursive algorithm
	f = Constant([0.0,0.0]) 
 	a = inner(alpha(u, k), v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)


	# Solve system of linear equations_______________________________________________________________
	A, b = assemble_system(a, L, bcs)
	solve(A, UP.vector(), b, "lu")
		
	U, P = UP.split()
	plot(U, interactive = True, title= "Velocity")
	plot(P, interactive = True, title= "Pressure")
		
 

	# be aware, with an cont func we have isolines, not with a discont func. 
	









	# d1 = something, d2 = something
	# eps = something 
	# if abs(new_array) > bla bla I should do some testes:
	#	y = eps 
	# else etc. 


 	
	# Update step:
	# K.vector()[:] -= project(indicator_function, K.function_space()).vector()

	# depending on funcitonspace the .vector() DOF array have different lengths.
	# P2 functions also info on nodes and edges, p1 on the nodes only.
	# extract array of Kn-1 numpyarray of all values. same array of indicator function
	# in order to make it valid, the array has to be the same lenght, 
	# therefore project first into the same space. 

	# ideally the optimization will contruct the values of the d's for us

	# do eq (8)
	# look at old geo and compare to new and se if it makes sense. 
	

	#epsilon = sym(grad(U)) # equivalent to the epsilon function
	#plot(epsilon) # Error, plotting of higher order functions are not supported
	


	





 	return U, P






def main():
	from dolfin import * 
	w   = 0.3
	
	# Expression for a Straight Tube 	
	f1  = Expression('(x[0] > w && x[0] < 1 - w ) ? 0.0 : 1.0', w=w)   
	# Expression for a T shape
	f2  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 && x[0] > w )? 0.0:1.0', w=w)	
	# Expression for a L shape 
	f3  = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7  )? 0.0:1.0', w=w) # K in math.formulation	
	
	
	U, P = stokes_solver(w=w, func=f3)
	#epsilon(U)
	



if __name__ == '__main__':
	from dolfin import * 	
	main()

	
