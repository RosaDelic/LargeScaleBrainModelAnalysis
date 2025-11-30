# Large-scale Brain Model Analysis

Understanding the dynamics of large-scale brain models remains a central challenge due to the inherent complexity of these systems.
In this work, we explore the emergence of complex spatiotemporal patterns in a large scale-brain model composed of 90 interconnected brain regions coupled through empirically derived anatomical connectivity. An important aspect of our formulation is that the  local dynamics of each brain region are described by a next-generation neural mass model, which explicitly captures the macroscopic gamma activity of coupled excitatory and inhibitory neural populations (PING mechanism). 
We first identify the system’s homogeneous states—both resting and oscillatory—and analyze their stability under uniform perturbations.
Then, we determine the stability against non-uniform perturbations by obtaining dispersion relations for the perturbation growth rate. This analysis enables us to link unstable directions of the homogeneous solutions to the emergence of rich spatiotemporal patterns, that we characterize by means of Lyapunov exponents and frequency spectrum analysis.
Our results show that, compared to previous studies with classical neural mass models, next-generation neural mass models provide a broader dynamical repertoire, both within homogeneous states and in the heterogeneous regime. Additionally, we identify a key role for anatomical connectivity in  cross-frequency coupling, allowing for the emergence of gamma oscillations with amplitude modulated by slower rhythms. These findings suggest that such models are not only more biophysically grounded but also particularly well-suited to capture the full complexity of large-scale brain dynamics. Overall, our study advances the analytical understanding of emerging spatiotemporal patterns in whole-brain models. 

We provide access to three folders with the codes that we implemented to carry out this work.

    1. Auto: Contains the codes that we implemented with the AUTO-07p software to study the homogeneous system. In this folder we provide all the necessary files for AUTO to run the homogeneous system of the large-scale brain where the dynamics of each node is modeled by means of next generation neural mass models. Also contains the files that we implemented to study the 1-dimensional bifurcation diagram as a function of the $I_{ext}^E$ parameter and the 2-dimensional bifurcation diagram as a function of the two parameters $(I_{ext}^E,\varepsilon)$.
    
    2. Julia: Contains the codes that we implemented to perform simulations of both the homogeneous system and the large-scale brain model. Also contains files with examples of how these systems should be integrated and how did we compute the Lyapunov exponents. We also include one file with the procedure that we used to compute the Power Spectrum of the signals.
    
    3. Python: Contains the codes that we implemented to perform simulations of the uncoupled next generation neural mass model, the homogeneous system and the large-scale brain model. We also provide here the routine that we implemented to study the transverse instabilities that arise in the model.


In all the Julia and Python folders, we also include the NormalizedMatrix.npz file that contains the normalized structural connectivity matrix that we used for our simulations.
