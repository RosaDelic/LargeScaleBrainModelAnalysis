import numpy as np
#import numba as nb
import mpmath as mp
mp.dps = 50 #Establish the number of precision digits to represent pi (in our case)
#from numba import jit, njit, prange


def NextGenPopModel(t0, x, W, Nvariables, Npop, param):
        #Function to build Large Brain Scale model

        #InputsModel:
            #t0: (Float) Time
            #x: (1x(NpopxNvariables) Float) Array defining the vector field x=np.array([r_e,v_e,s_e,r_i,v_i,s_i])
            #W: (NpopxNpop Integer) Normalised structural connectivity matrix
            #Nvariables: (Integer) Number of variables of the model
            #Npop: (Integer) Number of populations considered in the network
            #Parameters: (1x15 Float) Array containing the parameters of the model
                #tau_e,tau_i, tau_se, tau_si: Time constants of the model
                #nu_e,nu_i: Baseline constant current for excitatory,inhibitory neurons
                #Delta_e,Delta_i: Mean neuron noise intensity over excitatory,inhibitory
                #Jpq: For p,q in {e,i} Synaptic strength between E-I populations
                #Iext_e,Iext_i: External currents
                #eps: Coupling strength between network populations
        
        #Vector to store the field that describes the model
        dx = np.zeros(Npop*Nvariables)


        #---------  VARIABLES MODEL  ---------
        
        #Unzip the positions of the variables in separated vectors
        idx_ve = 1
        idx_se = 2
        idx_ri = 3
        idx_vi = 4
        idx_si = 5
        r_e_vector = x[0:90]
        v_e_vector = x[90*idx_ve:90*idx_se]
        s_e_vector = x[90*idx_se:90*idx_ri]
        r_i_vector = x[90*idx_ri:90*idx_vi]
        v_i_vector = x[90*idx_vi:90*idx_si]
        s_i_vector = x[90*idx_si:90*6]

        #---------  PARAMETERS MODEL  ---------
        tau_e = param['tau_e']
        tau_i = param['tau_i']
        tau_se = param['tau_se']
        tau_si = param['tau_si']
        Delta_e = param['Delta_e']
        Delta_i = param['Delta_i']
        nu_e = param['nu_e']
        nu_i = param['nu_i']
        Jee = param['Jee']
        Jei = param['Jei']
        Jii = param['Jii']
        Jie = param['Jie']
        Iext_e = param['Iext_e']
        Iext_i = param['Iext_i']
        eps = param['eps']

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  DIFFERENTIAL FIELD MODEL EQUATIONS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #np.sum(W*s_e_vector,axis=1)
        coupling_term = W@s_e_vector
        dx[0:90] = (Delta_e/(tau_e*np.pi)+2*r_e_vector*v_e_vector)/tau_e
        dx[90*idx_ve:90*idx_se] = (v_e_vector**2+nu_e-(np.pi*r_e_vector*tau_e)**2+Iext_e+tau_e*Jee*s_e_vector+tau_e*eps*coupling_term-tau_e*Jei*s_i_vector)/tau_e
        dx[90*idx_se:90*idx_ri] = (-s_e_vector+r_e_vector)/tau_se
        dx[90*idx_ri:90*idx_vi] = (Delta_i/(tau_i*np.pi)+2*r_i_vector*v_i_vector)/tau_i
        dx[90*idx_vi:90*idx_si] = (v_i_vector**2+nu_i-(np.pi*r_i_vector*tau_i)**2+Iext_i+tau_i*Jie*s_e_vector+tau_i*eps*coupling_term-tau_i*Jii*s_i_vector)/tau_i
        dx[90*idx_si:90*6] = (-s_i_vector+r_i_vector)/tau_si

        return dx