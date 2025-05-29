from matplotlib.transforms import Bbox
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

import Lyapunov
from DynamicalSystem import ContinuousDS
from HomogeneousSystemLyapunov import HomogeneousSystemLyapunov
from JacobianHomogeneous import JacobianHomogeneousSystemLyapunov
from Lyapunov import LCE

def LyapunovExponents(Nvariables,params,save):
    #Define some parameters
    # Initial time
    t0 = 0
    # Initial condition
    x0 = np.zeros(6)
    # Number of variables of the system
    Nvariables = 6
    #Discretization used for integration
    h = 0.001

    # Create system object to compute Lyapunov exponents
    SystemObject = ContinuousDS(x0, t0, HomogeneousSystemLyapunov, JacobianHomogeneousSystemLyapunov, h, Nvariables = Nvariables,param=params)

    Iext_e_vector = np.arange(4,12.1,step=0.1)

    Lyapunov_vector = np.zeros((Nvariables,len(Iext_e_vector)))

    print(Iext_e_vector.shape)


    for idx in range(len(Iext_e_vector)):
        print(idx)
        #Update params with current Iext_e
        params['Iext_e'] = Iext_e_vector[idx]

        # Need to import this again
        from Lyapunov import LCE
        # Computation of LCE
        LCE, history = LCE(SystemObject, 6, 0, 10**6, True)

        #Save Lyapunov exponents for current Iext_e
        Lyapunov_vector[:,idx] = LCE

    return Lyapunov_vector

