# 0a. Import requiered Python libraries
import numpy as np
import numba as nb
import pandas as pd
import scipy as sp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import mpmath as mp

# 0b. Load FloquetExponents data
dataPO = np.load('FloquetBifDiagra.npz')
dataPO

vector_Iext_e = dataPO["vector_Iext_e"]
vector_Iext_e = vector_Iext_e[0:314]
vector_eps = dataPO["vector_eps"]
vector_eps = vector_eps[0:314]
dataFloquetReal = dataPO["dataFloquetReal"]
dataFloquetReal = dataFloquetReal[0:314,0:661,:]
dataFloquetImaginary = dataPO["dataFloquetImaginary"]
dataFloquetImaginary = dataFloquetImaginary[0:314,0:661,:]
dataStatus = dataPO["dataStatus"]
dataStatus = dataStatus[0:314,0:661,:]

status2 = np.transpose(np.where(dataStatus == 2,dataStatus,np.nan))
status0 = np.transpose(np.where(dataStatus == 0,dataStatus,np.nan))

# 0c. Load Vaps data
dataEqPoints = np.load('EqPointsBifDiagram.npz')
dataEqPoints

dataVapsReal = dataEqPoints["dataVapsReal"]
dataVapsReal = dataVapsReal[0:314,0:661,:]
dataVapsImaginary = dataEqPoints["dataVapsImaginary"]
dataVapsImaginary = dataVapsImaginary[0:314,0:661,:]

# 1a. Treat Periodic Orbits
indices = np.zeros((len(vector_eps),len(vector_Iext_e)))
maxFloquet = np.zeros((len(vector_eps),len(vector_Iext_e)))
imagFloquet = np.zeros((len(vector_eps),len(vector_Iext_e)))
Floquet_2 = np.zeros((len(vector_eps),len(vector_Iext_e)))
Floquet_4 = np.zeros((len(vector_eps),len(vector_Iext_e)))
Floquet_6 = np.zeros((len(vector_eps),len(vector_Iext_e)))
Floquet_20 = np.zeros((len(vector_eps),len(vector_Iext_e)))
Floquet_90 = np.zeros((len(vector_eps),len(vector_Iext_e)))

for i in range(len(vector_eps)):
    for j in range(len(vector_Iext_e)):
        Floquet = dataFloquetReal[j,i,:]
        #print(Floquet[0])
        maxFloquet[i,j] = np.amax(Floquet[1:len(Floquet)])

        if not np.isnan(maxFloquet[i,j]):
            indices[i,j] =  np.argmax(Floquet[1:len(Floquet)])
            imagFloquet[i,j] = dataFloquetImaginary[j,i,int(indices[i,j])]
        else:
            indices[i,j] = np.nan
            imagFloquet[i,j] = np.nan
        
        Floquet_2[i,j] = Floquet[1]
        Floquet_4[i,j] = Floquet[3]
        Floquet_6[i,j] = Floquet[5]
        Floquet_20[i,j] = Floquet[19]
        Floquet_90[i,j] = Floquet[89]

indices = indices+1

Floquet_2_pos = np.where(Floquet_2 >= 0, Floquet_2, np.nan)
Floquet_4_pos = np.where(Floquet_4 >= 0, Floquet_4, np.nan)
Floquet_6_pos = np.where(Floquet_6 >= 0, Floquet_6, np.nan)
Floquet_20_pos = np.where(Floquet_20 >= 0, Floquet_20, np.nan)
Floquet_90_pos = np.where(Floquet_90 >= 0, Floquet_90, np.nan)

# Copy indices
indicesPositive = np.copy(indices)

# Create a mask where maxFloquet is positive
mask = maxFloquet >= 0

# Set the elements in A where the mask is false to NaN
indicesPositive[~mask] = np.nan

positive_directions = np.zeros((len(vector_eps),len(vector_Iext_e)))
for i in range(len(vector_eps)):
    for j in range(len(vector_Iext_e)):
        Floquet = dataFloquetReal[j,i,1:90]

        positive_directions[i,j] = np.sum((Floquet >= 0) & ~np.isnan(Floquet))
positive_directions[~mask] = np.nan

np.savez('TransverseInstabilitiesFloquet.npz',vector_Iext_e=vector_Iext_e,vector_eps=vector_eps, status0 = status0,status2=status2,maxFloquet=maxFloquet,positive_directions=positive_directions,Floquet_2_pos=Floquet_2_pos,Floquet_4_pos=Floquet_4_pos,Floquet_6_pos=Floquet_6_pos,Floquet_20_pos=Floquet_20_pos,Floquet_90_pos=Floquet_90_pos)

# 1b. Treat Equilibrium Points

indicesVaps = np.zeros((len(vector_eps),len(vector_Iext_e)))
maxVaps = np.zeros((len(vector_eps),len(vector_Iext_e)))
imagVaps = np.zeros((len(vector_eps),len(vector_Iext_e)))
Vaps_2 = np.zeros((len(vector_eps),len(vector_Iext_e)))
Vaps_4 = np.zeros((len(vector_eps),len(vector_Iext_e)))
Vaps_6 = np.zeros((len(vector_eps),len(vector_Iext_e)))
Vaps_20 = np.zeros((len(vector_eps),len(vector_Iext_e)))
Vaps_90 = np.zeros((len(vector_eps),len(vector_Iext_e)))


for i in range(len(vector_eps)):
    for j in range(len(vector_Iext_e)):
        current_vaps = dataVapsReal[j,i,:]
        #print(Floquet[0])
        maxVaps[i,j] = np.amax(current_vaps[1:len(current_vaps)])

        if maxVaps[i,j] != 0:
            indicesVaps[i,j] =  np.argmax(current_vaps[1:len(current_vaps)])
            imagVaps[i,j] = dataVapsImaginary[j,i,int(indicesVaps[i,j])]
        else:
            indicesVaps[i,j] = np.nan
            imagVaps[i,j] = np.nan
        
        Vaps_2[i,j] = current_vaps[1]
        Vaps_4[i,j] = current_vaps[3]
        Vaps_6[i,j] = current_vaps[5]
        Vaps_20[i,j] = current_vaps[19]
        Vaps_90[i,j] = current_vaps[89]
        
indicesVaps = indicesVaps+1

Vaps_2_pos = np.where(Vaps_2 > 0, Vaps_2, np.nan)
Vaps_4_pos = np.where(Vaps_4 > 0, Vaps_4, np.nan)
Vaps_6_pos = np.where(Vaps_6 > 0, Vaps_6, np.nan)
Vaps_20_pos = np.where(Vaps_20 > 0, Vaps_20, np.nan)
Vaps_90_pos = np.where(Vaps_90 > 0, Vaps_90, np.nan)

# Copy indices
indicesPositiveVaps = np.copy(indicesVaps)

# Create a mask where maxFloquet is positive
maskVaps = maxVaps >= 0

# Set the elements in A where the mask is false to NaN
indicesPositiveVaps[~mask] = np.nan

positive_directionsVaps = np.zeros((len(vector_eps),len(vector_Iext_e)))
for i in range(len(vector_eps)):
    for j in range(len(vector_Iext_e)):
        Vaps = dataVapsReal[j,i,1:90]

        positive_directionsVaps[i,j] = np.sum((Vaps >= 0) & ~np.isnan(Vaps))
positive_directionsVaps[~maskVaps] = np.nan

np.savez('TransverseInstabilitiesVaps.npz',maxVaps=maxVaps,positive_directionsVaps=positive_directionsVaps,Vaps_2_pos=Vaps_2_pos,Vaps_4_pos=Vaps_4_pos,Vaps_6_pos=Vaps_6_pos,Vaps_20_pos=Vaps_20_pos,Vaps_90_pos=Vaps_90_pos)

