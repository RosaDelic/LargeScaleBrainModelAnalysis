# Julia codes

1. HomogenousSystem.jl: Homogeneous system describing the dynamics on the homogeneous invariant manifold of the network

2. Network.jl: Large-scale brain model, where each node is an E-I network that we model with next generation neural mass models.

3. HomogeneousSystemSims.jl: Contains the routine that we used to integrate the homogeneous system and plot the results.

4. NetworkSims.jl: Contains the routine that we used to integrate the large-scale brain model and plot the results.

5. FourierAnalysis: Contains the routine that we used to compute the Power Spectrum of the signals with which we were working. We include the Power Spectrum computation of each node and the computation for the mean.

6. LyapunovRoutineHomogeneous.jl: Routine to compute the Lyapunov exponents of the homogeneous system across the parameter space within the chaotic region.

7. LyapunovRoutineNetwork.jl: Routine to compute the Lyapunov exponents of the large-scale brain model across all the parameter space.

8. PlotTransverseInstabilities.jl: Routine to plot the transverse instabilities data generated with TreatTransverseInstabilities.py

9. NormalizedMatrix.npz: Npz file containing the normalized structural connectivity matrix.
