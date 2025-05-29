import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect
from InitCondPOHomogeous import InitCondPOHomogeneous
from FloquetExponentsVariationals import ComputeFloquetExponents
from RoutineFloquet import RoutineFloquet
from PlotFloquetExponents import PlotFloquetExponents
from InitPO_2 import InitCondPOHomogeneous_2

def MaximumFloquetsContinuation(Iext_e,vapsConn,Nvariables,Npop,ModelParams):
    #Update parameter Iext_e
    ModelParams['Iext_e'] = Iext_e
    #Compute FloquetExponents and Multipliers for current pair of parameters (Iext_e,eps)
    matrixFloquetExp,matrixFloquetMult = RoutineFloquet(vapsConn,Nvariables,Npop,ModelParams)

    #Save maximum of the REAL part among the Nvariables FloquetExponents per each eigenvalue of W
    vectorMaxRealFloquetExp  = np.amax(np.real(matrixFloquetExp),axis=0)[1:Npop]
    
    #Compute maximum of all of them
    max_mu1 = np.amax(vectorMaxRealFloquetExp)

    return max_mu1

def my_bisection(f, a, mu_a, b, mu_b, tol, args): 
    # Approximates a root, R, of f bounded 
    # by a and b to within tolerance 
    # | f(m) | < tol with m the midpoint 
    # between a and b Recursive implementation
    print('******* Inside bisection *******')
    
    # Check if a and b bound a root
    if np.sign(mu_a) == np.sign(mu_b):
        raise Exception(
         "The scalars max_mu1 and max_mu2 do not bound a root")
        
    # Get midpoint
    c = (a + b)/2
    # Compute value of f at midpoint
    mu_c = f(c,args[0],args[1],args[2],args[3])
    
    if np.abs(mu_c) < tol:
        # stopping condition, report m as root
        return c
    elif np.sign(mu_a) == np.sign(mu_c):
        # case where c is an improvement on a. 
        # Make recursive call with a = c
        return my_bisection(f, c, mu_c, b, mu_b, tol, args)
    elif np.sign(mu_b) == np.sign(mu_c):
        # case where c is an improvement on b. 
        # Make recursive call with b = c
        return my_bisection(f, a, mu_a, c, mu_c, tol,args)

def ContinuationRoutine(vapsConn,Nvariables,Npop,ModelParams):
    #Set step to do the continuation
    h = 0.1

    #Initialize epsilon and Iext_e (approx from images) ---> 
    # FirstLeft: eps: 2.5//Iext_e1: 16 FirstUp: eps:21.465//Iext_e1:?  
    # FirstRight: eps: 19// Iext_e1: 8.3473 --> Until 12.5
    # OriginalsSecond: eps: 12.6//Iext_e1:12.32 (Is the same as the first right)
    # Second: eps Right: eps: 20//Iext_e: 9.485  --> Until 24
    # Second eps Left: eps: 19.38//Iext_e:8.309

    eps = 18.12857143
    Iext_e1 = 9.72423437

    #Initialize array of epsilons
    vector_eps = np.linspace(eps,18.85,int((18.85-eps)/h))  #---> For tomorrow: eps=24.081632653061224
    #vector_eps = np.linspace(eps,15,int((eps-15)/h))
    #Initialize array of Iext_e 
    vector_Iext_e = np.zeros(len(vector_eps))

    #Update parameters model with networks params for current simulation
    ModelParams['Iext_e'] = Iext_e1
    ModelParams['eps'] = eps

    #Compute FloquetExponents and Multipliers for current pair of parameters (Iext_e,eps)
    matrixFloquetExp,matrixFloquetMult = RoutineFloquet(vapsConn,Nvariables,Npop,ModelParams)

    #Save maximum of the REAL part among the Nvariables FloquetExponents per each eigenvalue of W
    vectorMaxRealFloquetExp  = np.amax(np.real(matrixFloquetExp),axis=0)[1:Npop]
    #Compute maximum of all of them
    max_mu = np.amax(vectorMaxRealFloquetExp)

    print('Maximum Floquet Exponent starter: '+str(max_mu))

    for idx in range(1,len(vector_eps)):
        #Set current eps
        eps = vector_eps[idx]
        print('================================================')
        print('epsilon: ' +str(eps))
        #Update parameter epsilon
        ModelParams['eps'] = eps

        #Compute maximum Floquet exponent among the Npop
        max_mu1 = MaximumFloquetsContinuation(Iext_e1,vapsConn,Nvariables,Npop,ModelParams)
        print('max_mu1: '+str(max_mu1))

        #Check if the max Floquet exponent for the previous Iext_e is positive or negative
        if max_mu1 > 0:
            #We need to increase Iext_e
            Iext_e2 = Iext_e1 + h
            dir = 1
        elif max_mu1 < 0:
            #We need to decrease Iext_e
            Iext_e2 = Iext_e1 - h
            dir = -1
        else:
            print('Error mu is zero')
            break

        cross = False
        while not cross:
            print('------ Inside cross-------')
            #Compute maximum Floquet exponent among the Npop
            max_mu2 = MaximumFloquetsContinuation(Iext_e2,vapsConn,Nvariables,Npop,ModelParams)
            print('Approx Iext_i: '+str(Iext_e2))
            print('max_mu2: '+str(max_mu2))

            if abs(max_mu2)<1e-5:
                cross = True
                vector_Iext_e[idx] = Iext_e2
                Iext_e1 = Iext_e2
            elif max_mu1*max_mu2 < 0:
                cross = True
                if dir == -1:
                    Output = my_bisection(MaximumFloquetsContinuation, Iext_e2, max_mu2, Iext_e1, max_mu1, tol=1e-5, args=(vapsConn,Nvariables,Npop,ModelParams))
                elif dir == 1:
                    Output = my_bisection(MaximumFloquetsContinuation, Iext_e1, max_mu1, Iext_e2, max_mu2, tol=1e-5, args=(vapsConn,Nvariables,Npop,ModelParams))
                else: 
                    print('Error the direction is not defined')
                    break
                print('Searching direction: '+str(dir))
                vector_Iext_e[idx] = Output
                Iext_e1 = Output
                #Call Bisection method with ComputeFloquetExp function and parameters Iext1,Iext2 --> Output = vector_Iext_e[idx]
                #Iext_e1 = Output
            else:
                #Keep increasing or decreasing Iext_e until we find a change of sign
                Iext_e1 = Iext_e2
                Iext_e2 = Iext_e2 + dir*h

        print(' DEFINITIVE Iext_e: ' +str(vector_Iext_e[idx]))

    return vector_eps,vector_Iext_e


