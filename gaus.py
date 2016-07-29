# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 19:12:51 2016

@author: Mahdiye
"""

import numpy as np
h=0.05
N=int(2/h)
G = np.zeros((N,1))
for i in range(N):
      t_0=0.05/4
      t=1
      x= -1.0 + i*h
      g = np.sqrt(t_0/t)*np.exp(-x**2/(0.05))
      G[i] = g
        
        
print G