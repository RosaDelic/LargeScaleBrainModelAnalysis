# LyapunovRoutineHomogeneous.jl

module LyapunovRoutineHomogeneous

    # Some imports
    using ChaosTools
    using DifferentialEquations
    using DynamicalSystems
    using StaticArrays
    using JLD2
    include("HomogeneousSystem.jl")

    # Define the system
    function LyapunovHomogeneousSystem!()
        # Initial eps
        eps0 = 9
        # Final Iext_e
        epsf = 15
        # Initial Iext_e
        I0 = 8
        # Final Iext_e
        If = 12
        # Discretization for the mesh
        h = 0.01
        # Create a vector for eps with step h from I0 to If (inclusive)
        vectorIext_e = range(I0, stop=If, step=h);
        # Create a vector for Iext_e with step h from I0 to If (inclusive)
        vector_eps = range(eps0, stop=epsf, step=h);

        # Parameters for Lyapunov exp function
        dt = 0.01
        ns = 0.1
        tmax = 1e4
        trans = 8e2
        Nsteps = Int(tmax/ns)

        # Integrate homogeneous system to set initial condition for Lyapunov Exponentsu0 = zeros(6)

        # Time span
        # Initial time
        t0 = 0
        # Final time
        tf = 10000
        tspan = (t0, tf)

        # Initial conditions (size should match the system, here we have a 2D system)
        uHom = zeros(6);

        # Define the parameters as a dictionary
        p = Float64[8,8,1,5,-5,-5,1,1,5,13,5,13,0,I0,eps0];

        # Initialize an array of SVectors, where each row is an SVector{6}
        matrix = Matrix{MVector{6, Float64}}(undef, length(vector_eps), length(vectorIext_e))

        # Loop through the range 1:length(vector_eps) and 1:length(vectorIext_e) to fill the matrix with 6 Lyapunov Exponents in each position
        for col in axes(matrix,2)
            for row in axes(matrix,1)
                # Set current Iext_e
                current_I = vectorIext_e[col]
                
                # Set current eps
                current_eps = vector_eps[row]
                
                println("eps = $current_eps, Iext_e=$current_I")

                # Update current Iext_e in parameters
                p[14] = current_I

                # Update current eps in parameters
                p[15] = current_eps

                # Compute LyapunovExponents
                ds = ContinuousDynamicalSystem(HomogeneousSystem.HomogeneousSystem!, u0, p; diffeq=(alg=RK4(),dt=dt,adaptive=true))
                tands = TangentDynamicalSystem(ds; k=6,J=HomogeneousSystem.JacobianHomogeneous!)
                λ = lyapunovspectrum(tands, Nsteps; Δt = ns, Ttr=trans)

                # Store the SVector in the idx-th position of the matrix
                matrix[row,col] = λ
            end
        end
        # Save all data
        # Save them to a JLD2 file
        @save "ImprovedLyapunovHomogeneous.jld2" matrix vector_eps vectorIext_e

        return vector_eps,vectorIext_e,matrix
    end
end