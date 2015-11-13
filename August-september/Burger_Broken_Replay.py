# -*- coding: utf-8 -*-

# Import statements------------------------------------------------------------------
from dolfin import *
from dolfin_adjoint import * # after dolfin, because it overrides some functionality

# The mesh over the geometry----------------------------------------------------------
n = 30
mesh = UnitSquareMesh(n, n)
V = VectorFunctionSpace(mesh, "CG", 2)

ic = project(Expression(("sin(2*pi*x[0])", "cos(2*pi*x[1])")),  V)

# The Variational problem-------------------------------------------------------------
def main(nu):
    u = Function(ic)
    u_next = Function(V)
    v = TestFunction(V)

    timestep = Constant(0.01)

    F = (inner((u_next - u)/timestep, v)
       + inner(grad(u_next)*u_next, v)
       + nu*inner(grad(u_next), grad(v)))*dx

    bc = DirichletBC(V, (0.0, 0.0), "on_boundary")

    parameters['adjoint']['test_derivative'] = True

    t = 0.0
    end = 0.1
    while (t <= end):
        solve(F == 0, u_next, bc)
        u.vector()[:] = u_next.vector()          # THIS IS THE ERROR , should be assign, dolfin cannot read this ?????????????
        t += float(timestep)

    return u

if __name__ == "__main__":
    nu = Constant(0.0001)
    u = main(nu)


# Replay function is supposed to detect this............................................
    success = replay_dolfin(tol=0.0, stop=True)
    # print f_21 Error message
"""

The output states that the norm of the difference(maybe the misfit) is greater then the tolerance, thus
there is an error in the code.THisis a local error I think. 

Comparing f_21:0:2:Forward against previously recorded value: norm of the difference is 6.092945e-02 (> tolerance of 0.000000e+00)


"""