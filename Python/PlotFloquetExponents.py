from RoutineFloquet import RoutineFloquet
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox

def PlotFloquetExponents(W,Nvariables,Npop,ModelParams,NetworkParams,save):

    #Inputs model: 
    #Nvariables: (Integer) Number of variables of the neural mass model
    #Npop: (Integer) Number of populations considered in the network
    #W: (NpopsxNpop) Normalized Structural Connectivity matrix
    #Parameters: (1x15 Float) Array containing the parameters of the model
        #tau_e,tau_i, tau_se, tau_si: Time constants of the model
        #nu_e,nu_i: Baseline constant current for excitatory,inhibitory neurons
        #Delta_e,Delta_i: Mean neuron noise intensity over excitatory,inhibitory
        #Jpq: For p,q in {e,i} Synaptic strength between E-I populations
        #Iext_i: External inhibitory current
        
    #ModelParams: (Tuple of float pairs) Contains pairs of the form (Iext_e,eps) for which compute the FloquetExponents
        #eps: Coupling strength between network populations

    #Compute eigenvalues and eigenvectors of connectivity matrix W
    vapsConn,vepsConn = np.linalg.eig(W)

    #Define tuple of invented vaps connectivity matrix
    vapsCurve = np.arange(-1,1+0.005,0.005)

    #Define colormap and vector of colors of the length of the tuple ModelParams
    colormap = cm.Set1
    colors_vector = colormap(np.linspace(0, 0.5, len(NetworkParams)))

    #Initialize matrices to store FloquetExponents
    matrixMaxRealFloquetExp = np.zeros((len(NetworkParams),Npop))
    matrixMaxImagFloquetExp = np.zeros((len(NetworkParams),Npop))
    InventedMatrixMaxRealFloquetExp = np.zeros((len(NetworkParams),len(vapsCurve)))

    for idx in range(len(NetworkParams)):
        #Get current pair of parameters for the simulation
        tuple = NetworkParams[idx]

        #Set color for plot
        #current_color = colors_vector[idx]

        #Get pairs of (Iext_e,eps) parameters for network simulation
        Iext_e = tuple[0]
        eps = tuple[1]
                
        #Update parameters model with networks params for current simulation
        ModelParams['Iext_e'] = Iext_e
        ModelParams['eps'] = eps

        #Compute FloquetExponents and Multipliers for current pair of parameters (Iext_e,eps)
        matrixFloquetExp,matrixFloquetMult = RoutineFloquet(vapsConn,Nvariables,Npop,ModelParams,len(vapsConn))

        #Compute Continuous curve of FloquetExponents and Multipliers for current pair of parameters (Iext_e,eps)
        InventedMatrixFloquetExp,InventedMatrixFloquetMult = RoutineFloquet(vapsCurve,Nvariables,Npop,ModelParams,len(vapsCurve))

        #Save maximum of the REAL part among the Nvariables FloquetExponents per each eigenvalue of W
        matrixMaxRealFloquetExp[idx,:]  = np.amax(np.real(matrixFloquetExp),axis=0)
        InventedMatrixMaxRealFloquetExp[idx,:]  = np.amax(np.real(InventedMatrixFloquetExp),axis=0)
        #Save indices corresponding to the previous maximum values
        indicesMaxReal = np.argmax(np.real(matrixFloquetExp),axis=0)
        #Save imaginary part of the previous maximum values
        matrixMaxImagFloquetExp[idx,:] = np.imag(matrixFloquetExp)[indicesMaxReal,np.arange(Npop)]

    fig1=plt.figure(figsize=(12,7))
    ax=plt.axes()
    plt.title(r'Real part Floquet Exponents $\epsilon = '+ f'{NetworkParams[idx][1]}$',fontsize=30,fontname='Times New Roman')
    plt.xlabel(r'$\mathcal{\Lambda}_{\alpha}$',fontsize=30,fontname='Times New Roman')
    plt.ylabel(r'$\mu$',fontsize=30,fontname='Times New Roman')
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    # Create custom legend handles for each row
    legend_handles = []
    for idx in range(len(NetworkParams)):
        plt.plot(vapsCurve,InventedMatrixMaxRealFloquetExp[idx,:],linewidth = 1,color=colors_vector[idx])
        plt.scatter(vapsConn,matrixMaxRealFloquetExp[idx,:],color=colors_vector[idx],facecolors='none',edgecolors=colors_vector[idx], s=30)  
   #     plt.plot(vapsConn,matrixMaxRealFloquetExp[idx,:],color=colors_vector[idx],linewidth=1)
        #legend_handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='none', markeredgecolor=colors_vector[idx], markersize=6, label=r'$I_{ext}^e = ' + f'{NetworkParams[idx][0]}$'))
    # Add the legend
    #plt.legend(handles=legend_handles,ncol=4, fontsize=15, title='Legend', title_fontsize=15,loc='lower left')
    plt.axhline(0,color="black", ls="-")
    plt.xlim([-1.01,1.01])
    plt.xticks([-1,-0.5,0,0.5,1])
    #plt.ylim([0,130])
    if save:
        plt.savefig('Definitive/CurvesRealFloquetExponents__'+str(eps)+'.png', dpi=500,bbox_inches=Bbox([[-1,-1],fig1.get_size_inches()]))
    plt.show()

    fig2=plt.figure(figsize=(12,7))
    ax=plt.axes()
    plt.title(r'Imaginary part Floquet Exponents $\epsilon = '+ f'{NetworkParams[idx][1]}$',fontsize=30,fontname='Times New Roman')
    plt.xlabel(r'$\mathcal{\Lambda}_{\alpha}$',fontsize=30,fontname='Times New Roman')
    plt.ylabel(r'$\beta$',fontsize=30,fontname='Times New Roman')
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    legend_handles = []
    for idx in range(len(NetworkParams)):
        plt.scatter(vapsConn,matrixMaxImagFloquetExp[idx,:],color=colors_vector[idx],facecolors='none',edgecolors=colors_vector[idx], s=30)   
        plt.plot(vapsConn,matrixMaxImagFloquetExp[idx,:],color=colors_vector[idx],linewidth=1)
        legend_handles.append(Line2D([0], [0], marker='o', color='w',markerfacecolor='none', markeredgecolor=colors_vector[idx], markersize=6, label=r'$I_{ext}^e = ' + f'{NetworkParams[idx][0]}$'))
    # Add the legend
    
    plt.legend(handles=legend_handles,ncol=4, fontsize=15, title='Legend', title_fontsize=15,loc='lower left')
    plt.axhline(0,color="black", ls="-")
    plt.xlim([-1.01,1.01])
    plt.xticks([-1,-0.5,0,0.5,1])
    #plt.ylim([0,130])
    #if save:
        #plt.savefig('Definitive/ImaginaryFloquetExponents__'+str(eps)+'.png',dpi=500,bbox_inches=Bbox([[-1,-1],fig1.get_size_inches()]))
    plt.show()

    return fig1,fig2