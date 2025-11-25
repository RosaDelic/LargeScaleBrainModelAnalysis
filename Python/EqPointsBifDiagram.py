# 0a. Load required Python libraries
import numpy as np
import numba as nb
import pandas as pd
import scipy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Load required routines and functions
from PlotVapsEqPoints import PlotVaps
from HomogeneousSystem import HomogeneousSystem
from EqPointsHomogeneousSystem import HomogeneousSystemZero,JacobianEqPoints

# 0b. Load Floquet Data obtained from FloquedBifDiagram.py
saved_data = np.load('FloquetBifDiagram.npz')
dataVector_Iext_e = saved_data['vector_Iext_e']
dataVector_Iext_e.shape
dataVector_eps = saved_data['vector_eps']
dataVector_eps.shape
dataStatus = saved_data['dataStatus']
dataStatus.shape

# 0c. Load structural connectivity matrix
data=np.load('NormalizedMatrix.npz')
norm_matrix = data['normalized_matrix']

# Compute eigenvalues and eigenvectors of connectivity matrix W
vapsConn,vepsConn = np.linalg.eig(norm_matrix) 

# 1a. Set parameters of the model
Nvariables = 6
Npop = 90
ModelParams = dict(tau_e = 8,
            tau_i = 8,
            tau_se=1,
            tau_si=5,
            nu_e = -5,
            nu_i = -5,
            Delta_e = 1,
            Delta_i = 1,
            Jee = 5,
            Jei = 13,
            Jii = 5,
            Jie = 13,
            Iext_i=0,
            Iext_e = 0,
            eps = 0)

# Initialize empty (long_Iext_e,long_eps,Npop) array structure to store Floquet exponents
dataVapsReal = np.zeros((dataStatus.shape[0]+1,dataStatus.shape[1]+1,Npop))
dataVapsImaginary = np.zeros((dataStatus.shape[0]+1,dataStatus.shape[1]+1,Npop))

# For loop over fixed points
for idx_Iext_e in range(dataStatus.shape[0]):
    for idx_eps in range(dataStatus.shape[1]):
        # Check that we are in a fixed point region
        if dataStatus[idx_Iext_e,idx_eps]==0:
            # Show progress
            print(' -----------  Iext_e: ',dataVector_Iext_e[idx_Iext_e],'eps: ',dataVector_eps[idx_eps],'  -----------')

            # Update parameters for current point
            ModelParams['Iext_e'] = dataVector_Iext_e[idx_Iext_e]
            ModelParams['eps'] = dataVector_eps[idx_eps]

            #Initialize empty NvariablesxNpop matrix to store Nvariables-Vaps of each dimension alpha 
            matrixVaps = np.zeros((Nvariables,Npop),dtype=np.complex128)

            # Define discretization
            h = 0.01

            # Set time of integration
            t0 = 0
            tf = 5000
            N = int((tf-t0)/h)
            time = np.arange(t0,tf,h)

            # Set initial condition
            x0 = np.zeros(6)

            # Integrate system a long time to tend to the equilibrium point
            sol = solve_ivp(HomogeneousSystem, [t0,tf], x0,t_eval=time, method='RK45', rtol=1e-6, atol=1e-9,args=(Nvariables,ModelParams))
            final_point = sol.y[:,len(sol.t)-1]

            # Compute zeros system to find the exact fixed point 
            eqPoint,infodict,ier,msg = fsolve(HomogeneousSystemZero, final_point,xtol=10**(-8), args=(Nvariables,ModelParams),full_output=True) #maxfev = 20
            print('Exit status: ',ier)
            print('Exit message: ',msg)
            print('Function evaluations: ',infodict['nfev'])
            print('Final residuals: ',infodict['fvec'])

            for idx in range(len(vapsConn)):
                #print('------ Iteration: '+str(idx)+' ------')
                #Set eigenvalue for current iteration
                eig = vapsConn[idx]

                #Compute Jacobian matrix for each eigenvalue
                jacobian = JacobianEqPoints(eqPoint,ModelParams,eig)

                #Compute eigenvalues and eigenvectors for each alpha
                vapsAlpha,vepsAlpha = np.linalg.eig(jacobian)

                #Save current 6 eigenvalues
                matrixVaps[:,idx] = vapsAlpha[:]
            
            # Save 90 maximum vaps per each pair parameter
            # Real part
            dataVapsReal[idx_Iext_e,idx_eps,:] = np.amax(np.real(matrixVaps),axis=0)

            # Imaginary part
            indicesMaxReal = np.argmax(np.real(matrixVaps),axis=0)
            dataVapsImaginary[idx_Iext_e,idx_eps,:] = np.imag(matrixVaps)[indicesMaxReal,np.arange(Npop)]
    # =============== SAVE ALL DATA  ================
    np.savez('EqPointsBifDiagram.npz',vector_Iext_e=dataVector_Iext_e,vector_eps=dataVector_eps,dataFloquetReal=dataVapsReal,dataFloquetImaginary=dataVapsImaginary,dataStatus=dataStatus)