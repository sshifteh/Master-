
from dolfin import * 
from indicator_function_expression import ind_func_expr
from K_update import K_update
from capping import capping
from stokes_solver import stokes_solver 
from WSS import WSS
from WSS_bdry import WSS_bdry  

n = 50
w=0.3
mesh = UnitSquareMesh(n,n)
Pspace = FunctionSpace(mesh, 'CG',1)
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace]) 
ShearStressSpace=FunctionSpace(mesh, 'DG', 1)
DPspace = FunctionSpace(mesh, 'DG', 0)	 


# strict inequality or not for 0.7 gives different results, but they are infact made the same by capping 
K_1 = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] <= 0.7 )? 0.0:1.0', w=w)


#------test 1--------------------------------

ind_func   = ind_func_expr(w=w, Pspace=Pspace) # interpolated to CG1
plot(ind_func, interactive = True, title = 'ind func created with an expre')

K_Func     = interpolate(K_1, Pspace) # control   
plot(K_Func, interactive = True, title = 'old K ') # old K  

K_Func_new = K_update(ind_func, K_Func) 
plot(K_Func_new, interactive = True, title = 'K updated') # updated K   

capping(K_Func_new)
K_capped   = capping(K_Func_new)
plot(K_capped, interactive = True, title = 'K capped') 



# Test 1 contionued with more calculations 
U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_capped, n=n )
WSS = WSS(U=U, ShearStressSpace= ShearStressSpace)


WSS_bdry = WSS_bdry(K_Func=K_Func, DPspace=DPspace, WSS=WSS)




#---- test 2------------------------------
# why isnt K_Func_new becoming this after being updated. 
# Begin by looking at the difference of the two plots?  The difference is that the one that doesnt work is projected 

#ind_f_exp = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] <= 0.7 || x[0] > 0.1 && x[0]<= 0.3 && x[1] >= 0.7 && x[1]< 0.8)? 0.0:1.0' , w=w)     
#stokes_solver(w=w,mesh=mesh,Vspace=Vspace,Pspace=Pspace,Wspace=Wspace,K_array=ind_f_exp, n=n )










