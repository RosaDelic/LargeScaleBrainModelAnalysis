# Large-scale Brain Model Analysis

Understanding the dynamics that emerge from large-scale brain models is challenging due to the high complexity of the system. In this project, we build upon the study in [1] to investigate the emergence of complex spatiotemporal patterns in a network of $90$ interconnected brain regions, each modeled as an excitatory-inhibitory (E-I) network whose dynamics are described by next generation neural mass models. We analyze the homogeneous oscillatory state of the system and study its stability under uniform perturbations. To assess its stability against non-uniform perturbations, we apply the Master Stability Function formalism, which allows us to characterize the emergence of complex spatiotemporal patterns from the unstable directions of the homogeneous state.

[1] Pau Clusella et al. “Complex spatiotemporal oscillations emerge from transverse instabilities in large-scale brain networks”. In: PLOS Computational Biology 

We provide access to three folders with the codes that we implemented to carry out this work.

    1. Auto: Contains the codes that we implemented with the AUTO-07p software to study the homogeneous system. In this folder we provide all the necessary files for AUTO to run the homogeneous system of the large-scale brain where the dynamics of each node is modeled by means of next generation neural mass models. Also contains the files that we implemented to study the 1-dimensional bifurcation diagram as a function of the $I_{ext}^E$ parameter and the 2-dimensional bifurcation diagram as a function of the two parameters $(I_{ext}^E,\varepsilon)$.
    
    2. Julia: Contains the codes that we implemented to perform simulations of both the homogeneous system and the large-scale brain model. Also contains files with examples of how these systems should be integrated and how did we compute the Lyapunov exponents. We also include one file with the procedure that we used to compute the Power Spectrum of the signals.
    
    3. Python: Contains the codes that we implemented to perform simulations of the uncoupled next generation neural mass model, the homogeneous system and the large-scale brain model. We also provide here the routine that we implemented to study the transverse instabilities that arise in the model.


In all the Julia and Python folders, we also include the NormalizedMatrix.npz file that contains the normalized structural connectivity matrix that we used for our simulations.
