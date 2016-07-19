# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:05:42 2015

@author: shiftehs
"""

from dolfin import * 
mesh = Mesh()
editor = MeshEditor()
editor.open(mesh,2,2)
editor.init_vertices(4)
editor.init_cells(2)
editor.add_vertex(0, 0.0, 0.0)
editor.add_vertex(1, 1.0, 0.0)
editor.add_vertex(2, 1.0, 0.0)
editor.add_vertex(3, 1.0, 1.0)
editor.add_cell(0, 0,1,2)
editor.add_cell(1, 0,2,3)
editor.close()
plot(mesh, interactive = True)
