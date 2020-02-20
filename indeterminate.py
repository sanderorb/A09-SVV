# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 11:49:05 2020

@author: joost
"""
import numpy as np
from math import cos,sin
import matplotlib.pyplot as plt
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.tools import FigureFactory as FF

import numpy as np
import pandas as pd
import scipy
theta=28
d_1=0.1
d_3=0.2
# x=[F_z1,F_z2,F_z3,F_y1,F_y2,F_y3,H,C1,C2,C3,C4,C5]
a = np.matrix(""" 0       1      0   ;
                 -0.0071 -0.111  0.12;
                  0       0.07  -0.3""")
print (a)

b=np.array([d_1*cos(theta),0,
                0,
                d_3*cos(theta),
                d_1*sin(theta),
                0,
                0,
               d_3*sin(theta)] )
print (b)
#a=0
#x= np.linalg.solve(a,b)
x1=2
x2=5
My=[]
def MCLY(F1,F2,x):
    
    My= F1*toggle(x,x1)**1
    return My
    
def toggle(x,loc):
    if x>=loc:
        return (x-loc)
    else:
        return 0
x=np.linspace(0,6.0,200)

for i in range (len(x)):
    My.append(MCLY(2.,4.,x[i]))
    
plt.plot(x,My)

def integrate(f,a,b):
    x = np.linspace(a,b,len(f)) # N+1 points make N subintervals
    y = f
    
    dx = x[1]-x[0]
    Int=[]
    for i in range (0,len(f)-1):
        T = (dx/2) * np.sum(y[i+1] + y[i])
        if i==0:
            Int.append(T)
        else:
             Int.append(T+Int[i-1])
    return Int,x
D,x=integrate(My,0,6.0)

print (x)
plt.plot(x[0:199],D)
plt.show()
plt.show()
