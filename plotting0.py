# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:05:41 2020

@author: joost
"""

# libraries
import numpy as np
import matplotlib.pyplot as plt
 
# Fake dataset
height = [3, 12, 5, 18, 45]
bars = ('A', 'B', 'C', 'D', 'E')
y_pos = np.arange(len(bars))
 
# Create bars and choose color
plt.bar(y_pos, height, color = (0.5,0.1,0.5,0.6))
 
# Add title and axis names
plt.title('My title')
plt.xlabel('categories')
plt.ylabel('values')
 
# Limits for the Y axis
plt.ylim(0,60)
 
# Create names
plt.xticks(y_pos, bars)
 
# Show graphic
plt.show()
