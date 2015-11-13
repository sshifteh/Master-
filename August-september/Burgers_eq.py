# -*- coding: utf-8 -*-

from dolfin import * 
from dolfin_adjoint import * # solve and assigne have been overloaded = erstattet tror jeg

# Making the mesh---------------------------------------------------
n = 30 
mesh =  UnitSquareMesh(n,n)
V = VectorFunctionSpace(mesh, 'CG', 2)

# Making the initial conditions ------------------------------------
ic = project(Expression(("sin(2*pi*x[0])", "cos(2*pi*x[1])")),  V) # initial condition discretized 
u = Function(ic)                                                   # Made to a Function Object
u_next = Function(V)                                               # an empty Function Object
v = TestFunction(V)                                                # TestFunction Object 

nu = Constant(0.0001)

timestep = Constant(0.01)

# Defining the variational problem---------------------------------
F = (inner((u_next - u)/timestep, v)
   + inner(grad(u_next)*u_next, v)
   + nu*inner(grad(u_next), grad(v)))*dx

bc = DirichletBC(V, (0.0, 0.0), "on_boundary")

t = 0.0
end = 0.1
while (t <= end):
    solve(F == 0, u_next, bc)
    u.assign(u_next)                                              # Dunno ?!
    t += float(timestep)                                          # iterating forward in time


# Using the dolfin-adjoint module for calculating the Gradient


J = Functional(inner(u,u)*dx*dt[FINISH_TIME])
dJdic = compute_gradient(J, Control(u), forget = False) # Forget False is cuz of deallocation which occurs with calling compute_gradient repeatedly
# Control does: Computes the Gradient(differentiates) wrt initial condition for that Function
# Assebmle each adjoint eq in turn 
# Uses the adjoint solution to compute the requested gradient 

dJdnu = compute_gradient(J, Control(nu))
# Computes gradient of J wrt nu


"""
shiftehs@shiftehs-Lenovo-IdeaPad-S300:~/Code$ python Burgers_eq.py 
Solving linear system of size 7442 x 7442 (PETSc Krylov solver).
No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.

Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 1.901e+00 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 1.162e-01 (tol = 1.000e-10) r (rel) = 6.115e-02 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 3.064e-03 (tol = 1.000e-10) r (rel) = 1.612e-03 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 4.929e-06 (tol = 1.000e-10) r (rel) = 2.593e-06 (tol = 1.000e-09)
  Newton iteration 4: r (abs) = 1.049e-11 (tol = 1.000e-10) r (rel) = 5.520e-12 (tol = 1.000e-09)
  Newton solver finished in 4 iterations and 4 linear solver iterations.
  
No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 1.033e-01 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 2.474e-03 (tol = 1.000e-10) r (rel) = 2.396e-02 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 3.787e-06 (tol = 1.000e-10) r (rel) = 3.667e-05 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 6.254e-12 (tol = 1.000e-10) r (rel) = 6.056e-11 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 9.208e-02 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 1.919e-03 (tol = 1.000e-10) r (rel) = 2.084e-02 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 1.412e-06 (tol = 1.000e-10) r (rel) = 1.533e-05 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 6.950e-13 (tol = 1.000e-10) r (rel) = 7.547e-12 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 8.375e-02 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 1.466e-03 (tol = 1.000e-10) r (rel) = 1.751e-02 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 4.557e-07 (tol = 1.000e-10) r (rel) = 5.441e-06 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 1.091e-13 (tol = 1.000e-10) r (rel) = 1.303e-12 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 7.844e-02 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 1.109e-03 (tol = 1.000e-10) r (rel) = 1.414e-02 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 3.213e-07 (tol = 1.000e-10) r (rel) = 4.096e-06 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 6.689e-14 (tol = 1.000e-10) r (rel) = 8.527e-13 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 7.509e-02 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 8.529e-04 (tol = 1.000e-10) r (rel) = 1.136e-02 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 2.660e-07 (tol = 1.000e-10) r (rel) = 3.542e-06 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 7.359e-14 (tol = 1.000e-10) r (rel) = 9.801e-13 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 7.272e-02 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 6.931e-04 (tol = 1.000e-10) r (rel) = 9.531e-03 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 2.546e-07 (tol = 1.000e-10) r (rel) = 3.501e-06 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 1.296e-13 (tol = 1.000e-10) r (rel) = 1.782e-12 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 7.093e-02 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 6.073e-04 (tol = 1.000e-10) r (rel) = 8.562e-03 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 2.689e-07 (tol = 1.000e-10) r (rel) = 3.791e-06 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 2.190e-13 (tol = 1.000e-10) r (rel) = 3.088e-12 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 6.958e-02 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 5.718e-04 (tol = 1.000e-10) r (rel) = 8.218e-03 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 3.139e-07 (tol = 1.000e-10) r (rel) = 4.512e-06 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 3.176e-13 (tol = 1.000e-10) r (rel) = 4.565e-12 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 6.862e-02 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 5.818e-04 (tol = 1.000e-10) r (rel) = 8.478e-03 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 4.676e-07 (tol = 1.000e-10) r (rel) = 6.814e-06 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 6.751e-13 (tol = 1.000e-10) r (rel) = 9.838e-12 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

No Jacobian form specified for nonlinear variational problem.
Differentiating residual form F to obtain Jacobian J = F'.
Solving nonlinear variational problem.
  Newton iteration 0: r (abs) = 6.806e-02 (tol = 1.000e-10) r (rel) = 1.000e+00 (tol = 1.000e-09)
  Newton iteration 1: r (abs) = 6.633e-04 (tol = 1.000e-10) r (rel) = 9.746e-03 (tol = 1.000e-09)
  Newton iteration 2: r (abs) = 7.493e-07 (tol = 1.000e-10) r (rel) = 1.101e-05 (tol = 1.000e-09)
  Newton iteration 3: r (abs) = 3.171e-12 (tol = 1.000e-10) r (rel) = 4.659e-11 (tol = 1.000e-09)
  Newton solver finished in 3 iterations and 3 linear solver iterations.

"""


