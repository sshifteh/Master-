# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 22:30:02 2015

@author: shiftehs
"""

from dolfin import * 

n = 5 
mesh = UnitSquareMesh(n,n)
plot(mesh, interactive = True)


class VesselTurn(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], )