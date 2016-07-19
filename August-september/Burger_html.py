# -*- coding: utf-8 -*-

# Import statements--------------------------------------------------------------------------------------
from dolfin import *
from dolfin_adjoint import *


# Making the code run in paralell with 8 theads in memory------------------------------------------------
parameters["num_threads"] = 8

# Checkpointing, dvs saving value of the forward step 11 time stages-------------------------------------
adj_checkpointing(strategy='multistage', steps=11,
                  snaps_on_disk=2, snaps_in_ram=2, verbose=True)
                  
# filename.py | grep Revolve to pipe the output from file to the program grep to search for Revolve                  
                  
                  
# Making the mesh over the geometry-----------------------------------------------------------------------
n = 30
mesh = UnitSquareMesh(n, n)
V = VectorFunctionSpace(mesh, "CG", 2)


# Solving the Variational Problem--------------------------------------------------------------------------
ic = project(Expression(("sin(2*pi*x[0])", "cos(2*pi*x[1])")),  V)

def main(nu):
    u = Function(ic)
    u_next = Function(V)
    v = TestFunction(V)

    timestep = Constant(0.01)

    F = (inner((u_next - u)/timestep, v)
       + inner(grad(u_next)*u_next, v)
       + nu*inner(grad(u_next), grad(v)))*dx

    bc = DirichletBC(V, (0.0, 0.0), "on_boundary")

    # Solving the PDE for each timestep, starting from time 0.0 ending at time 0.1
    t = 0.0
    end = 0.1
    while (t <= end):
        solve(F == 0, u_next, bc)
        u.assign(u_next)
        t += float(timestep)
        adj_inc_timestep()               # The incrementatin of timestep according to the checkpointing scheme ? 

    return u

if __name__ == "__main__":
    nu = Constant(0.0001)
    u = main(nu)


# Sending the code to html visualization-----------------------------------------------------------------------------
    adj_html("forward.html", "forward")
    adj_html("adjoint.html", "adjoint")

# Calculating the Functional 
    J = Functional(inner(u, u)*dx*dt[FINISH_TIME])
# Calculating the gradient of the functional, the sensitivity
    dJdnu = compute_gradient(J, Control(nu))
# Taking the integral
    Jnu = assemble(inner(u, u)*dx)

# This is dont know what does.
    parameters["adjoint"]["stop_annotating"] = True
    def Jhat(nu):
        u = main(nu)
        return assemble(inner(u, u)*dx)
# Performing a Taylor remainder convergence test 
    conv_rate = taylor_test(Jhat, Control(nu), Jnu, dJdnu)