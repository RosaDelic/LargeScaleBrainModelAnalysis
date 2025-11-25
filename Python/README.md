# Python codes

1. NextGenDynSyn.py: Next-Generation model.

2. NextGenPopModel.py: Large-scale brain model, where each node is an E-I network that we model with next generation neural mass models and they are coupled.

3. HomogeneousSystem.py: Homogeneous system equations.

4. HomogeneousVectorField.py: Vector field of the homogeneous system, used to compute equilibrium points.

5. JacobianHomogeneous.py: Jacobian of the homogeneous system, used to compute stability of equilibrium points.

6. InitPO_2.py: Routine to find initial condition on a periodic orbit.

7. FloquetExponentsVariationals.py: Routine to integrate the variational equations of the NextGenDynSyn model to compute the Floquet exponents and Floquet multipliers for a given parameter combination $(I_{ext}^E,\varepsilon)$.

8. VariationalsHomogeneous.py: System of the variational equations of the NextGenDynSyn model.

9. ChaoticDynamicsStudy.py: Routine to generate the chaotic cascade figure.

10. RoutineEqPoints.py: Routine to compute the 90 equilibrium points with maximum real part of the large-scale brain model for a given parameter combination $(I_{ext}^E,\varepsilon)$. Take the real and imaginary parts.

11. PlotVapsEqPoints.py: Routine to plot the real and imaginary parts of the 90 equilibrium points computed with RoutineEqPoints.py.

12. RoutineFloquet.py: Routine to compute the 90 Floquet exponents with maximum real part of the large-scale brain model for a given parameter combination $(I_{ext}^E,\varepsilon)$. Take the real and imaginary parts.

13. PlotFloquetExponents.py: Routine to plot the real and imaginary parts of the 90 Floquet exponents computed with RoutineFloquet.py.

14. FloquetBifDiagram.py: Routine to compute the 90 Floquet exponents with maximum real part among the considered range of parameters $(I_{ext}^E,\varepsilon)$ detecting where the homogeneous system tends to a periodic orbit and discarding the remaining regimes.

15. EqPointsBifDiagram.py: Routine to compute the 90 Floquet exponents with maximum real part among the considered range of parameters $(I_{ext}^E,\varepsilon)$ using the data obtained in FloquetBifDiagram.py corresponding to equilibrium points of the homogeneous system.

16. TreatTransverseInstabilities.py: Routine to treat the data from transverse instabilities generated with FloquetBifDiagram.py and EqPointsBifDiagram.py.

17. Videos.py: Routine to generate the videos from the simulations obtained from FourierAnalysis.jl.

18. NormalizedMatrix.npz: Npz file containing the normalized structural connectivity matrix.

19. aal.dat: File containing the spatial distribution of the nodes in the large-scale brain model.
