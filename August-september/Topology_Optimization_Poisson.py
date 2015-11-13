# -*- coding: utf-8 -*-

# Import statements -------------------
from dolfin import * 
from dolfin_adjoint import * 

# Importing the python interface to IPOPT
try:
    import pyipopt
except ImportError:
    info_red("""This example depends on IPOPT and pyipopt. \
  When compiling IPOPT, make sure to link against HSL, as it \
  is a necessity for practical problems.""")
    raise
    
# turn off abundant output
parameters["std_out_all_processes"] = False