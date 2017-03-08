import numpy as np
from scipy.sparse import csr_matrix
import operator
from itertools import chain
import time
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from matplotlib import collections  as mc
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import types
import pickle
import ico_prop_projection as sp
import Network_Projection_Appended_modified as c

s = sp.Sphere( recursion_level = 6 )
conn = c.Define_Connections(s)
colours, vconn, hconn, pent_ind = conn.define_connections()
myfiles=np.zeros((1,20480))
for i in range(10):
    namepickle='heat_maps'+str(i)+'.pkl'
    fileobj1 = open(namepickle, 'rb')
    f=pickle.load(fileobj1)
    fileobj1.close()
    myfiles=np.vstack((myfiles,f))
myfiles=myfiles[1:]
    
    
    