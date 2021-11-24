import numpy as np
import array as arr
l=arr.array("i",[50,25,10]) #n of elements per side
u=np.ndarray([l[0],l[1],l[2],3], dtype=float,order='F') #displacement
sigma=np.ndarray([l[0],l[1],l[2],3,3], dtype=float,order='F') #cauchy stress tensor

#BCtype=0->Dirichlet 
#BCtype=1->Force

BCtype=("i",[[0,1], #x: [neg face, pos face]
            [1,1],  #y:\\
            [1,1]]) #z:\\
#BCtype(axis,face)=0,1 ; axis in {0,1,2}; face in {0,1}

#Z AXIS - FACE
if BCtype(2,0)==0:   
    for i in range(l[0]):
        for j in range(l[1]):
            u(i,j,0,0)=0  #x #i,j,0 or 1? f**k python indexes LOL
            u(i,j,0,1)=0  #y #check everywhere
            u(i,j,0,2)=0  #z
elif BCtype(2,0)==1:   
    for i in range(l[0]):
        for j in range(l[1]):
            sigma(i,j,0,2,0)=0
            sigma(i,j,0,2,1)=0
            sigma(i,j,0,2,2)=0
#Z AXIS + FACE
if BCtype(2,1)==0:  
    for i in range(l[0]):
        for j in range(l[1]):
            u(i,j,l[2]-1,0)=0
            u(i,j,l[2]-1,1)=0
            u(i,j,l[2]-1,2)=0
elif BCtype(2,1)==1:     
    for i in range(l[0]):
        for j in range(l[1]):
            sigma(i,j,l[2]-1,2,0)=0
            sigma(i,j,l[2]-1,2,1)=0
            sigma(i,j,l[2]-1,2,2)=-100 #[N] axial, from top     

#X AXIS - FACE
if BCtype(0,0)==0: 
    for i in range(l[1])
        for j in range(l[2])
            u(0,i,j,0)=0  #x
            u(0,i,j,1)=0  #y 
            u(0,i,j,2)=0  #z
elif BCtype(0,0)==1:
    for i in range(l[1])
        for j in range(l[2])
            sigma(0,i,j,0,0)=0
            sigma(0,i,j,0,1)=0
            sigma(0,i,j,0,2)=0
#X AXIS + FACE
if BCtype(0,1)==0:  
    for i in range(l[1])
        for j in range(l[2])
            u(l[0]-1,i,j,0)=0  #x
            u(l[0]-1,i,j,1)=0  #y 
            u(l[0]-1,i,j,2)=0  #z
elif BCtype(0,1)==1:  
    for i in range(l[1])
        for j in range(l[2])
            sigma(l[0]-1,i,j,0,0)=0
            sigma(l[0]-1,i,j,0,1)=0
            sigma(l[0]-1,i,j,0,2)=0

#Y AXIS - FACE
if BCtype(1,0)==0:     
    for i in range(l[0])
        for j in range(l[2])
            u(i,0,j,0)=0  #x
            u(i,0,j,1)=0  #y 
            u(i,0,j,2)=0  #z
elif BCtype(1,0)==1:    
    for i in range(l[0])
        for j in range(l[2])
            sigma(i,0,j,0,0)=0
            sigma(i,0,j,0,1)=0
            sigma(i,0,j,0,2)=0
#Y AXIS + FACE
if BCtype(1,1)==0:      
    for i in range(l[0])
        for j in range(l[2])
            u(i,l[1]-1,j,0)=0  #x
            u(i,l[1]-1,j,1)=0  #y 
            u(i,l[1]-1,j,2)=0  #z
elif BCtype(1,1)==1:    
    for i in range(l[0])
        for j in range(l[2])
            sigma(i,l[1]-1,j,0,0)=0
            sigma(i,l[1]-1,j,0,1)=0
            sigma(i,l[1]-1,j,0,2)=0

#HOW TO
#deal with the problem that all the -not imposed values- 
#must not be given at GJ as 0, but as unknowns!!!
#And the -imposed values- must be given at GJ also if 0!!!