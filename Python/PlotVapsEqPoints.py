from RoutineEqPoints import RoutineVaps
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox

def PlotVaps(W,Nvariables,Npop,ModelParams,NetworkParams,save):

    # Enable LaTeX rendering
    plt.rcParams.update({
        "text.usetex": True,  # Use LaTeX for all text
        "font.family": "serif",
        "font.serif": ["Computer Modern"],
    })

    #Compute eigenvalues and eigenvectors of connectivity matrix W
    vapsConn,vepsConn = np.linalg.eig(W)

    #Define tuple of invented vaps connectivity matrix
    vapsCurve = np.arange(-1,1+0.005,0.005)

    #Define colormap and vector of colors of the length of the tuple ModelParams
    colormap = cm.Set1
    colors_vector = colormap(np.linspace(0, 0.5, len(NetworkParams)))

    #Initialize matrices to store FloquetExponents
    matrixMaxRealVaps = np.zeros((len(NetworkParams),Npop))
    matrixMaxImagVaps = np.zeros((len(NetworkParams),Npop))
    InventedMatrixMaxRealVaps = np.zeros((len(NetworkParams),len(vapsCurve)))
    InventedMatrixMaxImagVaps = np.zeros((len(NetworkParams),len(vapsCurve)))

    for idx in range(len(NetworkParams)):
        #Get current pair of parameters for the simulation
        tuple = NetworkParams[idx]

        #Get pairs of (Iext_e,eps) parameters for network simulation
        Iext_e = tuple[0]
        eps = tuple[1]
                
        #Update parameters model with networks params for current simulation
        ModelParams['Iext_e'] = Iext_e
        ModelParams['eps'] = eps

        #Compute Vaps for current pair of parameters (Iext_e,eps) ---> Returns a 6xlenVapsConn matrix
        matrixVaps = RoutineVaps(vapsConn,Nvariables,Npop,ModelParams,len(vapsConn))

        #Compute Continuous curve of Vaps for current pair of parameters (Iext_e,eps)
        InventedMatrixVaps = RoutineVaps(vapsCurve,Nvariables,Npop,ModelParams,len(vapsCurve))

        #Save maximum of the REAL part among the Nvariables Vaps per each eigenvalue of W
        matrixMaxRealVaps[idx,:]  = np.amax(np.real(matrixVaps),axis=0)
        InventedMatrixMaxRealVaps[idx,:]  = np.amax(np.real(InventedMatrixVaps),axis=0)

        #Save indices corresponding to the previous maximum values
        indicesMaxReal = np.argmax(np.real(matrixVaps),axis=0)
        InventedIndicesMaxReal = np.argmax(np.real(InventedMatrixVaps),axis=0)

        #Save imaginary part of the previous maximum values
        matrixMaxImagVaps[idx,:] = np.imag(matrixVaps)[indicesMaxReal,np.arange(Npop)]
        InventedMatrixMaxImagVaps[idx,:] = np.imag(InventedMatrixVaps)[InventedIndicesMaxReal,np.arange(len(vapsCurve))]
    
    fig1=plt.figure(figsize=(12,7))
    ax=plt.axes()
    plt.title(r'$\varepsilon = '+ f'{NetworkParams[idx][1]}$',fontsize=40,fontname='Times New Roman',loc="left",pad=20)
    plt.xlabel(r'$\Lambda_{\alpha}$',fontsize=40,fontname='Times New Roman')
    plt.ylabel(r'$\mu_{max}^{(\alpha)}$',fontsize=40,fontname='Times New Roman')
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    # Create custom legend handles for each row
    legend_handles = []
    for idx in range(len(NetworkParams)):
        plt.plot(vapsCurve,InventedMatrixMaxRealVaps[idx,:],linewidth = 1,color=colors_vector[idx])
        plt.scatter(vapsConn,matrixMaxRealVaps[idx,:],color=colors_vector[idx],facecolors='none',edgecolors=colors_vector[idx], s=30)
        legend_handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='none', markeredgecolor=colors_vector[idx], markersize=6, label=r'$I_{ext}^E = ' + f'{NetworkParams[idx][0]}$'))
    # Add the legend
    plt.legend(handles=legend_handles,ncol=4, fontsize=22,loc='upper right', bbox_to_anchor=(1.04, 1.16),columnspacing=0.5,handletextpad=0.5,frameon=True)
    plt.axhline(0,color="black", ls="-")
    plt.xlim([-1.01,1.01])
    plt.xticks([-1,-0.5,0,0.5,1])
    plt.yticks([-0.3,-0.2,-0.1,0,0.1])
    plt.ylim(-0.3,0.15)
    if save:
        plt.savefig('Vaps__'+str(eps)+'.png', dpi=500,bbox_inches="tight")
    plt.show()

    fig2=plt.figure(figsize=(12,7))
    ax=plt.axes()
    plt.title(r'Imaginary part Vaps $\varepsilon = '+ f'{NetworkParams[idx][1]}$',fontsize=30,fontname='Times New Roman')
    plt.xlabel(r'$\mu_{max}^{(\alpha)}$',fontsize=30,fontname='Times New Roman')
    plt.ylabel(r'$\beta$',fontsize=30,fontname='Times New Roman')
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    legend_handles = []
    for idx in range(len(NetworkParams)):
        plt.plot(vapsCurve,InventedMatrixMaxImagVaps[idx,:],linewidth = 1,color=colors_vector[idx])
        plt.scatter(vapsConn,matrixMaxImagVaps[idx,:],color=colors_vector[idx],facecolors='none',edgecolors=colors_vector[idx], s=30)   
        legend_handles.append(Line2D([0], [0], marker='o', color='w',markerfacecolor='none', markeredgecolor=colors_vector[idx], markersize=6, label=r'$I_{ext}^E = ' + f'{NetworkParams[idx][0]}$'))
    # Add the legend
    plt.legend(handles=legend_handles,ncol=2, fontsize=15, title='Legend', title_fontsize=15,loc='lower left')
    plt.axhline(0,color="black", ls="-")
    plt.xlim([-1.01,1.01])
    plt.xticks([-1,-0.5,0,0.5,1])
    if save:
        plt.savefig('ImaginaryVaps__'+str(eps)+'.png',dpi=500,bbox_inches=Bbox([[-1,-1],bbox_inches="tight")
    plt.show()

    return fig1,fig2
