# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 10:59:34 2015

@author: shiftehs
"""

from dolfin import * 


class RightBranch(SubDomain):
    """    
    Subclass of SubDomain.
    Object type CellFunction    
    """
    def inside(self, x, on_boundary):
        return  on_boundary and x[1] > 0.5 - DOLFIN_EPS and x[0] > 0.5 - DOLFIN_EPS\
               and not (x[1] < x[0] - w)\
               and not (x[1] > x[0] + w)


class LeftBranch(SubDomain):
    """   
    Subclass of SubDomain.
    Object type CellFunction    
    """
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 0.5 - DOLFIN_EPS and x[0] < 0.5 + DOLFIN_EPS\
               and not (x[1] < -x[0] + 1 - w)\
               and not (x[1] > -x[0] + 1 + w)
               
class Stem(SubDomain):
    """    
    Subclass of SubDomain.
    Object type CellFunction    
    """
    def inside(self, x, on_boundary):
        return on_boundary and x[1] < 0.5 + DOLFIN_EPS and \
                   0.5-w-DOLFIN_EPS < x[0] < 0.5+w+DOLFIN_EPS 
