# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import numpy as np
import ico_propagation as sp
import Network_clean as Network

import pickle
from matplotlib import pyplot as plts

#Kishan'f Data
kishmatrix = np.array([[0.02,0.99981,4.3015e-06],
[0.04,0.99983,3.8088e-06],
[0.06,0.9998,1.0454e-05],
[0.08,0.99968,3.0663e-05],
[0.11,0.99772,0.00044859],
[0.13,0.96099,0.018246],
[0.15,0.60984,0.054379],
[0.17,0.16381,0.041092],
[0.19,0.017807,0.0080603],
[0.21,0.020737,0.016513],
[0.23,4.922e-05,4.8685e-05],
[0.25,0.0001084,8.4968e-05],
[0.27,0,0],
[0.29,0,0],
[0.12,0.99152,0.0027053],
[0.14,0.86184,0.028043],
[0.16,0.29714,0.055185],
[0.18,0.039206,0.013863],
[0.2,0.0056277,0.0028284],
[0.22,4.834e-05,3.6005e-05],
[0.24,0.00082172,0.00081342],
[0.26,0,0],
[0.28,0,0],
[0.3,9.406e-05,9.3115e-05],
[0.1,0.99919,0.00010423]])


kishnu = kishmatrix[:,0]
kishrisk = kishmatrix[:,1]
kisherror = kishmatrix[:,2]
p_transversalconn= kishnu 




    