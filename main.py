import numpy as np
import array as arr
dX=arr.array("i",[1,1,1]) #size of the element
X=np.ndarray([l[0],l[1],l[2],3], dtype=float,order='F') #position

epsilon=np.ndarray([l[0],l[1],l[2],3,3], dtype=float,order='F') #print(sigma.shape)
FV=np.ndarray([l[0],l[1],l[2],3], dtype=float,order='F')

#Boundary Conditions
from BC import *

#EQUATION SETS:
