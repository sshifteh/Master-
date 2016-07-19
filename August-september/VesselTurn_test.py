# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 00:23:08 2015

@author: shiftehs
"""

import numpy as np 
import matplotlib.pyplot as plt 
import scitools.std as sci

x = np.linspace(1,10,100)
y = np.cos(x)
sci.plot(x,y)
sci.xaxis = [0,3]
sci.show()  

