# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 11:00:22 2015

@author: shiftehs
"""

from dolfin import * 


class LeftWall(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0)

class RightWall(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 1.0)


class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0) and not Stem().inside(x, on_boundary)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0) and not RightBranch().inside(x, on_boundary) and not LeftBranch().inside(x, on_boundary) 

































