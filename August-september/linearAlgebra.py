# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 16:33:24 2015

@author: shiftehs
"""

from dolfin import *

x = Vector()  # Creating a vector 
A = Matrix()  # Creating a matrix
#x2 = Vector(10) # a vector of sixe 10 

#print A 
#solve (A,x,b, 'lu', 'ilu')

list_krylov_solver_methods()

list_krylov_solver_preconditioners()