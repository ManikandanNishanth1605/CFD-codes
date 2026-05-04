#################################---IMPORTS---########################################
########## Local Imports #########
import laplace as lp
######## 3rd Party Imports #######
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

##################################---INPUTS---########################################
l_x = 1  # Breadth of the mesh
l_y = 2  # Length of the mesh
n_x = 21 # Number on nodes along breadth
n_y = 41 # Number of nodes along length
d_x = l_x/(n_x - 1) # x length
d_y = l_y/(n_y - 1) # y length

############################---INITIAL CONDITIONS---##################################
T_0 = 30 # Temperature all over the plate (but not needed for this problem)

###########################---BOUNDARY CONDITIONS---##################################
T_upper = 30 # Temperature of the upper surface
T_left = 30 # Temperature of the left surface
T_right = 30 # Temperature of the right surface
T_bottom = 100 # Temperature of the bottom surface

#################################---SOLVING---########################################
laplaceInstance = lp.Laplace(l_x,l_y,n_x,n_y,T_0)  # Create a class instance of Laplace
laplaceInstance.generateMesh()  # Create a rectangular mesh
laplaceInstance.applyBoundaryConditions(T_upper,T_left,T_right,T_bottom)  # Apply boundary conditions to the mesh

## 1. Point Gauss Siedel
# resultPGSM,itersPGSM,errorPGSM,timePGSM = laplaceInstance.pointGaussSiedel()

## 2. Line Gauss Siedel
# resultLGSM,itersLGSM,errorLGSM,timeLGSM = laplaceInstance.lineGaussSiedel()

## 3. Point Successive Over Relaxation
resultPSOR,itersPSOR,errorPSOR,timePSOR = laplaceInstance.pSOR(1.2)

## 4. Line Successive Over Relaxation
# resultLSOR,itersLSOR,errorADI,timeLSOR = laplaceInstance.lSOR(1.01)

## 5. Alternating Direction Implicit Method
# resultADI,itersADI,errorADI,timeADI = laplaceInstance.adi()

## 6. Alternating Direction Implicit with SOR
# resultADISOR,itersADISOR,errorADISOR,timeADISOR = laplaceInstance.adiSOR(1.2)

#############################---POST PROCESSING---####################################

# Heat transfer along the bottom boundary
# We need to find the heat transfer along y-axis
result = resultPSOR # Result from solver
k = 50 # Conductivity
diff = np.array([]) # Array of forward differences (heat flux) at each nodes just above the boundary
for i in range(1,n_x-1):
    # Looping thorugh each nodes
    d = (result[n_y-2,i] - result[n_y-3,i])/d_y # Heat flux
    diff = np.append(diff,d)
# Computing heat transfer by average heat flux multiplied by area
Q_y = np.sum(diff)*l_x/(n_x-2)
print(Q_y)

##############################---END OF PROGRAM---####################################