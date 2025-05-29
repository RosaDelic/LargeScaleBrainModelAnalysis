import numpy as np
from scipy.integrate import solve_ivp
from HomogeneousSystem import HomogeneousSystem
from HomogeneousVectorField import HomogeneousVectorField
from VariationalsHomogeneous import VariationalsHomogeneous
import matplotlib.pyplot as plt
from functools import partial

def PoincareEvent(t,x,Nvar,Params,alpha):
    return x[0]-alpha

# Event should be triggered when g(x) == 0 <=> x[1]-alpha=0
#PoincareEvent.terminal = True  # Stop the integration when the event occurs
#PoincareEvent.direction = 0    # Trigger when the second component crosses alpha in any direction

def InitCondPOHomogeneous(Nvariables,params):

    ## ==============================  1.  FIND FIRST APPROXIMATION  ===================================

    # 1. Set parameters 
    #Initial and final times of integration
    t0 = 0
    tflong = 5000
    #Discretization used for integration
    h = 0.001
    #Number of points to evaluate the time integration
    N = int((tflong-t0)/h)
    #print(N)
    time_eval = np.linspace(t0, tflong, N)
    #Initial condition
    x0 = np.zeros(6)
    #Tolerance used to stop Poincare section
    tol=5*10**(-5)
    
    # 2. Integrate Homogeneous state along tf time
    sol = solve_ivp(HomogeneousSystem, [t0,tflong], x0,t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9,args=(Nvariables,params))

    # 3. Take final point --> should be on a periodic orbit and define Poincare section
    #Take as approximate initCond the last time of the integration
    x0_int = sol.y[:,len(sol.t)-1]
    #Take as Poincare section the second to last time of integration so when we start again, we do not stop at time t=0
    Poincare = sol.y[:,len(sol.t)-2]
    alpha_val = Poincare[0]

    # 4. Integrate the homogeneous system from initial condition x_approx up to the first crossing with the Poincare section

    # 4a. Redefine PoincareEvent with updated value of alpha
    PoincareEventAlpha = partial(PoincareEvent, alpha=alpha_val)
    PoincareEventAlpha.terminal = True  # Stop the integration when the event occurs
    #Set the direction of stopping 
    if sol.y[0,len(sol.t)-1]-alpha_val > 0: #We crossed the section with ascending direction 
        PoincareEventAlpha.direction = 1    # Trigger when the second component crosses alpha in ascending direction
        #direction = 1
    if sol.y[0,len(sol.t)-1]-alpha_val < 0: #We crossed the section with descending direction 
        PoincareEventAlpha.direction = -1    # Trigger when the second component crosses alpha in descending direction
        #direction = -1

    # 4b. Integrate system
    integrationT = solve_ivp(HomogeneousSystem, [t0,tflong], x0_int,t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9, events=[PoincareEventAlpha],args=(Nvariables,params))
    
    fig=plt.figure(figsize=(20,8))
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
    
    #Print result
    #print(integrationT)

    #First index 0 to get the element of the list corresp to the first event. 
    #Then for y_events we get a tx6-dim array 
    #For t_events we get a 1xt-dim array

    #print(integrationT.t_events[0].shape)
    #print(integrationT.y_events[0][:,0].shape)

    #print('Last value integration: '+ str(x_approx))
    #print('Zero intersection: '+str(integrationT.y_events[0][0,:]))
    #print('Difference in two consecutive solutions: ' +str(abs(integrationT.y_events[0][3,:]-integrationT.y_events[0][2,:])))
    #print('Approximation Period: ' + str(integrationT.t_events[0][3]-integrationT.t_events[0][2])+' ms.')


    # 5. Set approximated initial condition and period for the Newton's routine
    #x0 = integrationT.y_events[0][3,:]
    #T0 = integrationT.t_events[0][3]-integrationT.t_events[0][2]
    T0 = integrationT.t_events[0][0]
    while np.linalg.norm(Poincare-integrationT.y_events[0][0,:])>tol:
        print('=========  Entering while  ==========')
        print('Previous x_approx: '+str(x0_int))
        print('New x: '+str(integrationT.y_events[0][0,:]))
        print('Norm of difference: ' +str(np.linalg.norm(Poincare-integrationT.y_events[0][0,:])))
        print('EvalCondition: ' +str(np.linalg.norm(Poincare-integrationT.y_events[0][0,:])>tol))
        print('Previous period: ' +str(T0))

        #Update approximated initCondition and period
        x0_int = integrationT.y_events[0][0,:]
        T0 = integrationT.t_events[0][0]

        #Update the direction of stopping 
        #if direction == 1: #We crossed the section with ascending direction 
        #    PoincareEventAlpha.direction = 1    # Trigger when the second component crosses alpha in ascending direction
        #if direction == -1: #We crossed the section with descending direction 
        #    PoincareEventAlpha.direction = -1    # Trigger when the second component crosses alpha in descending direction

        # Refinement: integrate system
        integrationT = solve_ivp(HomogeneousSystem, [t0,tflong], x0_int,t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9, events=[PoincareEventAlpha],args=(Nvariables,params))
        
        fig=plt.figure(figsize=(20,8))
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
    
    #Set initial condition and period on the Periodic Orbit
    #Definitive initial condition
    x0 = integrationT.y_events[0][0,:]

    #Once we hace the condition on the Periodic Orbit, integrate again to find period
    tfshort = 250
    #Discretization used for integration
    h = 0.0001
    #Number of points to evaluate the time integration
    N = int((tfshort-t0)/h)
    #print(N)
    time_eval = np.linspace(t0, tfshort, N)
    PoincareEventAlpha.terminal = False  # Do not stop the integration when the event occurs
    integrationT = solve_ivp(HomogeneousSystem, [t0,tfshort], x0_int,t_eval=time_eval, method='RK45', rtol=1e-6, atol=1e-9, events=[PoincareEventAlpha],args=(Nvariables,params))
    
    fig=plt.figure(figsize=(20,8))
    ax=plt.axes()
    plt.plot(integrationT.t,integrationT.y[3,:],color='blue',linewidth=6,label='r_i')
    plt.plot(integrationT.t,integrationT.y[0,:],color='red',linewidth=6,label='r_e')
    plt.plot(integrationT.t_events[0][:],integrationT.y_events[0][:,0],'o',color='black',linewidth=10)
    plt.xlabel("Time (ms)",fontsize=40)
    plt.ylabel("Firing Rate",fontsize=40)
    plt.xticks(fontsize=40)#,rotation=45)
    plt.yticks(fontsize=40)
    ax.legend(loc="upper right",fontsize=30,ncol=2)
    plt.show()
    plt.close()

    #Definitive Period
    T0 = integrationT.t_events[0][1]-integrationT.t_events[0][0]

    print('=========  Final iteration  ==========')
    print('Norm of difference: ' +str(np.linalg.norm(Poincare-integrationT.y_events[0][0,:])))
    print('EvalCondition: ' +str(np.linalg.norm(Poincare-integrationT.y_events[0][0,:])>tol))
    print('InitialCondition on Periodic orbit: ' +str(x0))
    print('Period of Periodic orbit: ' +str(T0))

    return x0,T0


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



