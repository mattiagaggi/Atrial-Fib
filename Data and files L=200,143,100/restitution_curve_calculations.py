import numpy as np
from numpy import vectorize
import pickle
#import dill
import time
import types
import random
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import gridspec
from itertools import compress
import math
import copy


x=np.arange(0,800,3)
y=np.array([])
for elements in x:
    if elements<=210:
        y=np.append(y,90)
    if elements>210 and  elements<=400:
        y=np.append(y, ((130./190)*elements+(220-400*130/190)))
    if elements>400:
        y=np.append(y,220)

plt.plot(x,y)
plt.ylim(0,max(y)+0.1*max(y))
plt.xlabel("S1-S2 ms")
plt.ylabel("APD ms")
plt.title("APD vs S2-S1")


plt.figure()
z=y*150./220




kim=np.array([])
for elements in x:
    kim=np.append(kim,150)
plt.plot(x,z)
plt.plot(x,kim)
plt.ylim(0,max(z)+0.1*max(z))
plt.xlabel("S1-S2 ms")
plt.ylabel("refractory period ms")
plt.title("refractory period vs S2-S1")

plt.figure()
stepsx=x/3
stepsx=np.round(stepsx)
stepsx=stepsx.astype(int)
stepsz=z/3
stepsz=np.round(stepsz)
stepsz=stepsz.astype(int)

plt.plot(stepsx,stepsz)
plt.ylim(0,max(stepsz)+0.1*max(stepsz))
plt.xlabel("S1-S2 in steps")
plt.ylabel("refractory period in steps")
plt.title("refractory period vs S2-S1")

            
plt.show()
