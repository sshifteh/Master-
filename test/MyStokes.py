# -*- coding: utf-8 -*-
from dolfin import *

def StokesSolve(omega, mf):
    plot(omega)

    V = VectorFunctionSpace(omega, "CG", 2)
    Q = FunctionSpace(omega, "CG", 1)
    W = V*Q

    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)

    f1 = Constant((0.0, 0.0))

    a = -inner(nabla_grad(u),nabla_grad(v))*dx + p*div(v)*dx - div(u)*q*dx
    L = inner(f1,v)*dx

    bc1 = DirichletBC(W.sub(0), (0.0, 0.0), mf, 3)
    bc2 = DirichletBC(W.sub(0), (0.0, 0.0), mf, 4)
    bc3 = DirichletBC(W.sub(0), (1.0, 0.0), mf, 1)

    bc = [bc1, bc2, bc3]

    w = Function(W)
    solve (a == L, w, bc)
    (u,p) = w.split(True)


    J = 0.5*inner(grad(u),grad(u))*dx
    J = assemble(J)
    print "Objective function J = ", J
    return J, u, p
