import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from HomogeneousSystem import HomogeneousSystem
from EqPointsHomogeneousSystem import HomogeneousSystemZero,JacobianEqPoints


def RoutineVaps(vapsConn,Nvariables,Npop,params,size):
        
    #Inputs model: 
    #Nvariables: (Integer) Number of variables of the neural mass model
    #Npop: (Integer) Number of populations considered in the network
    #W: (NpopsxNpop) Normalized Structural Connectivity matrix
    #Parameters: (1x15 Float) Array containing the parameters of the model
        #tau_e,tau_i, tau_se, tau_si: Time constants of the model
        #nu_e,nu_i: Baseline constant current for excitatory,inhibitory neurons
        #Delta_e,Delta_i: Mean neuron noise intensity over excitatory,inhibitory
        #Jpq: For p,q in {e,i} Synaptic strength between E-I populations
        #Iext_e,Iext_i: External currents
        #eps: Coupling strength between network populations

    #Initialize empty NvariablesxNpop matrix to store Nvariables-Vaps of each dimension alpha 
    matrixVaps = np.zeros((Nvariables,size),dtype=np.complex128)

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
    sol = solve_ivp(HomogeneousSystem, [t0,tf], x0,t_eval=time, method='RK45', rtol=1e-6, atol=1e-9,args=(Nvariables,params))
    final_point = sol.y[:,len(sol.t)-1]

    # Compute zeros system to find the exact fixed point 
    eqPoint,infodict,ier,msg = fsolve(HomogeneousSystemZero, final_point,xtol=10**(-8), args=(Nvariables,params),full_output=True) #maxfev = 20
    print('Exit status: ',ier)
    print('Exit message: ',msg)
    print('Function evaluations: ',infodict['nfev'])
    # print('Number jacobian calls: ',infodict['njev'])
    print('Final residuals: ',infodict['fvec'])

    for idx in range(len(vapsConn)):
        #print('------ Iteration: '+str(idx)+' ------')
        #Set eigenvalue for current iteration
        eig = vapsConn[idx]

        #Compute Jacobian matrix for each eigenvalue
        jacobian = JacobianEqPoints(eqPoint,params,eig)

        #Compute eigenvalues and eigenvectors for each alpha
        vapsAlpha,vepsAlpha = np.linalg.eig(jacobian)

        #Save current 6 eigenvalues
        matrixVaps[:,idx] = vapsAlpha[:]


    return matrixVaps