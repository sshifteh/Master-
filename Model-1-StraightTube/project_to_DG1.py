from dolfin import *

mesh = UnitSquareMesh(10, 10)

cell_f = CellFunction('size_t', mesh, 0)
CompiledSubDomain('x[0] < 0.25+DOLFIN_EPS').mark(cell_f, 1) 
CompiledSubDomain('x[0] > 0.75-DOLFIN_EPS').mark(cell_f, 1) 

class MyExpression(Expression):
    def __init__(self, cell_f):
        self.cell_f = cell_f
    def eval_cell(self, value, x, cell):
        value[0] = self.cell_f[cell.index]

fexp = MyExpression(cell_f)

V = FunctionSpace(mesh, 'DG', 1)
f = interpolate(fexp, V)
# The claim is that for such f, (f.dx(i), v) = 0
v = TestFunction(V)
print 'Projecting from DG1 to DG1'
for i in range(2):
    b_form = inner(f.dx(i), v)*dx
    b = assemble(b_form)
    print 'dx(%d)' % i, b.norm('linf')

W = FunctionSpace(mesh, 'CG', 1)
f = interpolate(fexp, W)
# The claim is that for such f, (f.dx(i), v) != 0
print 'Projecting from CG1 to DG1'
for i in range(2):
    b_form = inner(f.dx(i), v)*dx
    b = assemble(b_form)
    print 'dx(%d)' % i, b.norm('linf')


def my_project(f, W):
    '''Project f to W function space'''
    # TODO
    return Function(W)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import numpy as np

    mesh = UnitSquareMesh(15, 15)
    V = FunctionSpace(mesh, 'CG', 2)
    f_vals = np.random.rand(V.dim())
    f = Function(V)
    f.vector()[:] = f_vals
    
    W = FunctionSpace(mesh, 'CR', 1)
    f0 = project(f, W)
    f1 = my_project(f, W)

    # Are they the same
    f0.vector().axpy(-1, f1.vector())
    print 'YES!' if near(f0.vector().norm('linf'), 1, 1E-13) else 'NO!'

# https://github.com/FEniCS/dolfin/blob/master/site-packages/dolfin/fem/norms.py




