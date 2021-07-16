import pandas as pd
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from ipywidgets import interactive, IntSlider, fixed, interact
import sys

data=np.loadtxt('mylin.dat')

for i in range(data.shape[0]):
    x=data[i,2]
    y=data[i,3]
    if data[data[:,0]==data[i,1]].shape[0]==1:
       xn=data[data[:,0]==data[i,1]][0,2]
       yn=data[data[:,0]==data[i,1]][0,3]
       plt.plot([x,xn],[y,yn],'-')
plt.show()
