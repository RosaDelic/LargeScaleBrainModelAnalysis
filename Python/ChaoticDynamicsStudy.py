import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from functools import partial
from scipy.signal import find_peaks
from HomogeneousSystem import HomogeneousSystem
from matplotlib.transforms import Bbox
import matplotlib.colors as mcolors

import Lyapunov
from DynamicalSystem import ContinuousDS
from HomogeneousSystemLyapunov import HomogeneousSystemLyapunov
from JacobianHomogeneous import JacobianHomogeneousSystemLyapunov
from Lyapunov import LCE


def ChaoticMaxMin(Nvariables,params,save):
    # Set colors dictionary
    color_dict = mcolors.CSS4_COLORS
    # Initialize Iext_e vector
    Iext_e_vector = np.arange(8.5,12,step=0.001)
    # Initialize impacts list
    impacts_list = []

    print(Iext_e_vector.shape)


    for idx in range(len(Iext_e_vector)):
        print(idx)
        # Update to current Iext_e parameter
        params['Iext_e'] = Iext_e_vector[idx]

        # 1. Initial integration to compute mean of s_e that we will use as Poincare section

        #Initial and final times of integration
        t0 = 0
        tf = 1000
        #Discretization used for integration
        h = 0.001
        #Number of points to evaluate the time integration
        N = int((tf-t0)/h)
        time_eval = np.linspace(t0, tf, N)
        #Initial condition
        x0 = np.zeros(6)
        #Integrate
        sol = solve_ivp(HomogeneousSystem, [t0,tf], x0,t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9, args=(Nvariables,params))

        # 2. Compute the mean s_e of the previous integration
        idx_ms = 200
        mean = np.mean(sol.y[2,len(sol.t)-int(idx_ms/h):len(sol.t)])
        max = np.amax(sol.y[2,len(sol.t)-int(idx_ms/h):len(sol.t)])
        #Final time
        tf = 500
        #Number of points to evaluate the time integration
        N = int((tf-t0)/h)
        time_eval = np.linspace(t0, tf, N)
        # 2a. Integrate Homogeneous state along tf time
        sol = solve_ivp(HomogeneousSystem, [t0,tf], sol.y[:,-1],t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9, args=(Nvariables,params))

        peaks,_ = find_peaks(sol.y[2,:])#,height=mean+0.01*max)

        impacts_list.append(np.unique(np.round(sol.y[2,peaks],decimals=5)))

        # 3. Reestructure data
        #Compute maximum dimension
        max_dim_impacts = np.amax([len(arr) for arr in impacts_list])

        #Create padded array of (max dimension,len_Iext_e_vector)
        padded_arrays = np.array([np.pad(arr, (0, max_dim_impacts - len(arr)), constant_values=np.nan) for arr in impacts_list])

    
    fig=plt.figure(figsize=(20,8))
    ax=plt.axes()
    #plt.plot(sol.t,sol.y[5,:],color='yellow',linewidth=6,label='s_i')    legend_handles = []
    for idx in range(max_dim_impacts):
        plt.scatter(Iext_e_vector,padded_arrays[:,idx],color=color_dict['crimson'],facecolors=color_dict['crimson'],edgecolors=color_dict['crimson'], s=2)   
        
    # Add the legend
    plt.xlabel(r'$I_{ext}^e$',fontsize=40)
    plt.ylabel(r'$max(s_E)$',fontsize=40)
    plt.xticks(fontsize=35)#,rotation=45)
    plt.yticks(fontsize=35)
    plt.ylim([0.14,0.32])
    plt.title(r'$\epsilon=$'+str(params['eps']),fontsize=40,fontname='Times New Roman')
    if save:
        plt.savefig('PeriodDoublings/Free/eps_'+str(params['eps'])+'.png',dpi=500,bbox_inches=Bbox([[0,-1],fig.get_size_inches()]))
    plt.show()

    return padded_arrays,fig

