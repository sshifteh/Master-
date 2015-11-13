# -*- coding: utf-8 -*-

from dolfin import * 
from dolfin_adjoint import *


J = Functional(inner(u,u)*dx*dt[FINITH_TIME])
reduced_functional = ReducedFunctional(J,m)
m_optimized = minimize(reduced_functional)  # minimerer 
m_optimized = maximize(reduced_funcitonal)  # maksimerer 

print_optimization_methods()

m_optimized = minimize(reduced_functional, method ='CG')


# Callback-------------------------------------------------------------------------------------------------------
# A callback to be called after every optimization iteration
def eval_cb(j,m):
    print 'functional value', j, 'functional gradient', float(m)

def derivative_cb(j,dj,m):
    print 'functional values:',j, 'functional gradient:', dj, 'parameter value:', float(m)

reduced_functional = ReducedFunctional(J, ConstantControl('nu'), eval_cb = eval_cb, derivative_cb = derivative_cb)





# Advanced more options for optimization method--------------------------------------------------------------------
m_optimized = minimize(reduced_functional,method = 'SLSQP', tol = 1.0E-8, options = {'disp' = True} )

# method = optimization metoden, CG er et alternativ, for flere alternativer 
# print_optimization_methods()

# tol =  tolerance for termination

# options = {'maxiter': 3, 'disp': True, 'gtol': True?}
# maxiter: is the max no of iterations to be perfomed
# disp : True if I want to print convergence messages 
# gtol: tol for the iteration loop til aa stoppe hvis gradienten normen er mindre enn tol

reduced_functional = ReducedFunctional(J, [m1, m2,....])
m_optimized = minimize(reduced_functional, bounds = ([m1_lb, m1_ub],[m2_lb, m2_ub])

m_optimized = minimize(reduced_functional,bounds = (m_lb, m_ub), method = 'CG', tol = 1.0E-7, options = {'maxiter'= 10, 'disp' = True, 'gtol' = 1.0E-8})
m_optimized = minimize(reduced_functional, bounds = (0.0, 1.0))

# In case the gradient is incorrect by we are not made aware of it by an error message from the program
# Dette sørger for at gradient evelueringen er indeed korrekt
# ved å kjøre en Taylor remainder test for hver gradient evaluation during optimiseringen   
dolfin.parameters['optimization']['test_gradient'] = True
dolfin.parameters['optimization']['gradient_seed'] = 0.001






















