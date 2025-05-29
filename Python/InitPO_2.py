import numpy as np
from scipy.integrate import solve_ivp
from HomogeneousSystem import HomogeneousSystem
from HomogeneousVectorField import HomogeneousVectorField
from VariationalsHomogeneous import VariationalsHomogeneous
import matplotlib.pyplot as plt
from functools import partial
from scipy.optimize import fsolve
from scipy.signal import find_peaks

# Define Poincare section
def PoincareEvent(t,x,Nvar,Params,alpha):
    return x[2]-alpha

# Event should be triggered when g(x) == 0 <=> x[1]=0
PoincareEvent.terminal = False  # Stop the integration when the event occurs
PoincareEvent.direction = 1    # Trigger when the second component crosses alpha in ascending direction


# Define the function that we want to find the zero of (phi(T, x0) - x0)
def SystemPeriodicity(initCond, Nvariables,params,alpha):
    #Unfold initCond
    T_approx = initCond[0]
    x0_approx = initCond[1:len(initCond)]

    #print(T_approx)

    # 1. Set parameters 
    #Initial and final times of integration
    t0 = 0
    tf = T_approx
    #Discretization used for integration
    h = 0.001
    #Number of points to evaluate the time integration
    N = int((tf-t0)/h)
    #print(N)
    time_eval = np.linspace(t0, tf, N)
    #Initial condition
    x0 = np.zeros(Nvariables)

    # Integrate the system from initial condition x0_guess for time T_guess
    sol = solve_ivp(HomogeneousSystem, [t0,tf], x0_approx, t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9, args=(Nvariables,params))
    
    # Compute the difference between the final state and the initial state   g(x)=phi(T,x)-alpha=0// phi(T,x)=x
    diff = np.append(sol.y[2, -1]-alpha,sol.y[:, -1]-x0_approx)

    return diff


def InitCondPOHomogeneous_2(Nvariables,params):

    ## ==============================  1.  FIND FIRST APPROXIMATION  ===================================
    # 1. Initial and final times of integration
    t0 = 0
    tf = 1500
    #Discretization used for integration
    h = 0.001
    #Number of points to evaluate the time integration
    N = int((tf-t0)/h)
    time_eval = np.linspace(t0, tf, N)
    #Initial condition
    x0 = np.zeros(Nvariables)

    # 1a. Initial integration to check that the model oscillates
    sol = solve_ivp(HomogeneousSystem, [t0,tf], x0,t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9, args=(Nvariables,params))

    # 1b. Check if there are oscillations
    # Look at the peaks of the v_i voltage eliminating the first part of the integration
    idx_ms = 100
    #peaks,_ = find_peaks(sol.y[4,len(sol.t)-int(idx_ms/h):len(sol.t)],height=1)
    difference = np.amax(sol.y[4,len(sol.t)-int(idx_ms/h):len(sol.t)])-np.amin(sol.y[4,len(sol.t)-int(idx_ms/h):len(sol.t)])
    
    #Here we check that the model oscillates for this parameters --> WE ARE SURE WITH THIS CONDITION
    #print(difference)
    if difference < 0.1:
        # Status 0 encodes no oscillations in the model
        status = 0
        return status, np.zeros(Nvariables), 0

    # 2. Compute the mean s_e of the previous integration
    mean = np.mean(sol.y[2,len(sol.t)-int(idx_ms/h):len(sol.t)])
    #print(mean)
    PoincareEventAlpha = partial(PoincareEvent, alpha=mean)
    PoincareEventAlpha.direction = 1  
    tf = 500
    #Number of points to evaluate the time integration
    N = int((tf-t0)/h)
    time_eval = np.linspace(t0, tf, N)
    # 2a. Integrate Homogeneous state along tf time with poincare events
    sol = solve_ivp(HomogeneousSystem, [t0,tf], sol.y[:,-1],t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9,events = PoincareEventAlpha, args=(Nvariables,params))

    #fig1=plt.figure(figsize=(12,7))
    #ax=plt.axes()
    #plt.xlabel('Time (ms)',fontsize=30,fontname='Times New Roman')
    #plt.ylabel(r'$s_E$',fontsize=30,fontname='Times New Roman')
    #plt.scatter(sol.t_events[0][:],sol.y_events[0][:,2],color='black')
    #plt.xticks(fontsize=30)
    #plt.yticks(fontsize=30)
    # Create custom legend handles for each row
    #plt.plot(sol.t,sol.y[2,:])
    #plt.axhline(mean,color="black", ls="-")


    # 2b. Check if the oscillations are periodic

    #print(sol.t_events[0])
    #period1 = sol.t_events[0][-3:][1]-sol.t_events[0][-3:][0]
    #period2 = sol.t_events[0][-3:][2]-sol.t_events[0][-3:][1]

    vector_periods = sol.t_events[0]
    diff_periods = np.diff(vector_periods)
    equal_diff_periods = np.diff(diff_periods)
    
    check = np.where(equal_diff_periods>5*10**(-3)*np.ones(len(equal_diff_periods)))

    #print(equal_diff_periods)
    #print(check)
    
    if len(check[0])>0:
        # Status 1 encodes chaotic region
        status = 1
        return status, np.zeros(Nvariables), 0

    else:
        # Status 2 encodes periodic oscillation region
        status = 2
    
        # 3. Take final point of the integration --> should be approx to periodic orbit
        x0_approx = sol.y_events[0][-1,:]
        T_approx = np.abs(sol.t_events[0][-2:][1]-sol.t_events[0][-2:][0])

        # Use fsolve to find the solution
        T_solution,infodict,ier,msg = fsolve(SystemPeriodicity, np.append(T_approx,x0_approx),xtol=10**(-8), args=(Nvariables,params,mean),full_output=True) #maxfev = 20
        print('Exit status: ',ier)
        print('Exit message: ',msg)
        print('Function evaluations: ',infodict['nfev'])
    # print('Number jacobian calls: ',infodict['njev'])
        print('Final residuals: ',infodict['fvec'])

        if ier != 1 or T_solution[0] <= 0:
            # Status 4 encodes convergence problems --> Ens quedem amb l'anterior
            status = 4
            return status, x0_approx, T_approx


        return status,T_solution[1:len(T_solution)],T_solution[0]


'''
    fig2=plt.figure(figsize=(20,8))
    ax=plt.axes()
    plt.plot(integrationT.t,integrationT.y[3,:],color='blue',linewidth=6,label='r_i')
    plt.plot(integrationT.t,integrationT.y[0,:],color='red',linewidth=6,label='r_e')
    plt.plot(integrationT.t_events[0],integrationT.y_events[0][:,0],'o',color='black',linewidth=10)
    plt.xlabel("Time (ms)",fontsize=40)
    plt.ylabel("Firing Rate",fontsize=40)
    plt.xticks(fontsize=40)#,rotation=45)
    plt.yticks(fontsize=40)
    ax.legend(loc="upper right",fontsize=30,ncol=2)
    plt.show()
    plt.close()
    

'''

'''
    ## ==============================  2. REFINEMENT NEWTON'S METHOD  ===================================
    print('-----------------  Refinement: Newton Method  ----------------')
    # 1. Set approximated initial condition and period for the Newton's routine
    x0 = integrationT.y_events[0][3,:]
    T0 = integrationT.t_events[0][3]-integrationT.t_events[0][2]

    print('Newton x0: '+ str(x0))
    print('Newton T0: ' +str(T0))

    #Tolerance for the method to converge and maximum number of iterations
    tol = 10^(-12)
    Nmax = 200
    Niter = 1
    while np.linalg.norm(x0-x_approx)>tol and Niter < Nmax:
        x_approx = x0
        T_approx = T0

        # 2. Initial condition of variational equations
        initVariationals = np.eye(Nvariables).flatten()

        # 3. Define parameters for integration
        #Number of points to evaluate the time integration up to period T0
        N = int((T0-t0)/h)
        time_eval = np.linspace(t0, T0, N)
        eig = 1
        #Complete initial condition for integration
        initCond = np.concatenate((x0,initVariationals),axis=0)

        # 5. Integrate the system up to the Poincare section
        integrationT = solve_ivp(VariationalsHomogeneous, [t0,T0], initCond, t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9, args=(Nvariables,eig,params))

        # 6. Monodromy matrix: Get solution of the variational equations at time T
        sol = integrationT.y[0:Nvariables,len(integrationT.t)-1]
        Monodromy = integrationT.y[Nvariables:Nvariables+Nvariables*Nvariables,len(integrationT.t)-1].reshape((Nvariables,Nvariables))

        # 7. Evaluate vector field on solution
        sol_vectorField = HomogeneousVectorField(sol,Nvariables,params)

        # 8. Concatenate solutions in big Jacobian
        J = np.concatenate((sol_vectorField.reshape((len(sol_vectorField),1)),Monodromy-np.eye(Nvariables)),axis=1)
        upper_vector = np.zeros((1,Nvariables+1))
        upper_vector[0,1] = 1
        J = np.concatenate((upper_vector,J),axis=0)

        # 9. Solve linear system to compute dk
        dk = -np.linalg.solve(J, b)



        Niter += 1


'''



