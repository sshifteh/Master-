
from dolfin import * 
from indicator_function_expression import ind_func_expr
from K_update import K_update
from capping import capping
from stokes_solver import stokes_solver 
from WSS import WSS
from WSS_bdry import WSS_bdry  
from indicator_funciton_conditions2 import ind_func_expr

n = 50
w=0.3
mesh = UnitSquareMesh(n,n)
Pspace = FunctionSpace(mesh, 'CG',1)
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace]) 
ShearStressSpace=FunctionSpace(mesh, 'DG', 1)
DPspace = FunctionSpace(mesh, 'DG', 0)	 



# Test 3 ------------------------------------
# I want to:::::::::::::::::::::::::::::::::: 
# step 1: need a K old 
# step 2: send it to Solver, 
# step 3: calculate the WSS over it and over its bdry 

# step 4:then we need to calculate an ind function with threshold values for addition domain to it 
# step 5:We need to add together the old and the new K functions 
# step 6:Then capping the new K function might be necessary 
# step 7:finally run stokes solver over the new K funciton 

# iterate 50 times 
# return the K_func 
# store it to xdmf file to open it paraview 




# step 1: Take K old 
K_1 = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 )? 0.0:1.0', w=w)

# step 2: Send it to Stokes solver 
U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_1, n=n )

# step 3: calculate the WSS over it and its bdry 
WSS = WSS(U=U, ShearStressSpace= ShearStressSpace)
WSS_bdry = WSS_bdry(K_Func=K_Func, DPspace=DPspace, WSS=WSS)

# Step 4: This is the diffucult step. we need to create an indicator function now 
ind_f = ind_func_expr(WSS_bdry = WSS_bdry)
# works!! finally 

# step 5: add together old and new K functions 
# ind_f is a function in DG0 space and old K function is inPspace
# have to solve this first 

K_new = K_update(ind_func = ind_f, K_Func = K_Func)
plot(K_new, interactive = True, title = 'K new ')
# doesnt look that great, lets try capping 

#step 6: capping 
K_capped = capping(K_new)
plot(K_capped, interactive = True, title = 'K capped')
# actually looks really good, 

# step 7: run fluid trough it 
U, P, K_Funcs = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_capped, n=n )






