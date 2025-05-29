# LyapunovRoutineNetwork.jl

module LyapunovRoutineNetwork

    # Some imports
    using ChaosTools
    using DifferentialEquations
    using DynamicalSystems
    using StaticArrays
    using JLD2
    using Distributions
    include("Network.jl")
    include("HomogeneousSystem.jl")

    # Define a parameter struct that contains the scalar parameters and the matrix
    struct Params
        scalar_params::Vector{Float64}  # Vector of scalar parameters (e.g., p[1], p[2], ...)
        matrix_params::Matrix{Float64}  # The normalized connectivity matrix
    end

    # Define the system
    function LyapunovNetwork!()
        # Initial eps
        eps0 = 9
        # Final Iext_e
        epsf = 9
        # Initial Iext_e
        I0 = 0
        # Final Iext_e
        If = 16
        # Discretization for the mesh
        h = 0.2
        # Create a vector for eps with step h from I0 to If (inclusive)
        vectorIext_e = range(I0, stop=If, step=h);
        # Create a vector for Iext_e with step h from I0 to If (inclusive)
        vector_eps = range(eps0, stop=epsf, step=h);

        # Parameters for Lyapunov exp function
        dt = 0.1
        ns = 1
        tmax = 7e3
        trans = 3e2 
        Nsteps = Int(tmax/ns)

        Nvariables = 6
        Npop = 90

        # Define the parameters as a dictionary
        scalar_params = [8,8,1,5,-5,-5,1,1,5,13,5,13,0,I0,eps0];
        
        # Integrate homogeneous system to set initial condition for Lyapunov Exponentsu0 = zeros(6)

        # Time span
        # Initial time
        t0 = 0
        # Final time
        tf = 10000
        tspan = (t0, tf)
        
        # Load the .npz file containing the connectivity matrix into a dictionary-like structure
        data = load("NormalizedMatrix.npz")
        W = data["normalized_matrix"]
        
        # Create the Params struct
        p = Params(scalar_params, W)

        # Initialize an array of SVectors, where each row is an SVector{6}
        matrix = Matrix{MVector{90, Float64}}(undef, length(vector_eps), length(vectorIext_e))

        # Loop through the range 1:length(vector_eps) and 1:length(vectorIext_e) to fill the matrix with 90x6 Lyapunov Exponents in each position
        for col in axes(matrix,2)
            for row in axes(matrix,1)
                
                # Set current Iext_e
                current_I = vectorIext_e[col]
                
                # Set current eps
                current_eps = vector_eps[row]
                
                println("eps = $current_eps, Iext_e=$current_I")

                # Update current Iext_e in parameters
                p.scalar_params[14] = current_I

                # Update current eps in parameters
                p.scalar_params[15] = current_eps

                # Initial conditions 
                uHom = zeros(6);
                # Define the problem using the system and parameters
                prob = ODEProblem(HomogeneousSystem.HomogeneousSystem!, uHom, tspan, p.scalar_params)
                # Solve the ODE problem 
                sol = solve(prob,reltol = 1e-6,abstol = 1e-6,saveat=0.01);
                
                # Initial conditions
                fin = sol[:,end]
                u0 = zeros(Npop*Nvariables)
                u0[1:Npop] = fill(fin[1],Npop)
                u0[1+Npop:2*Npop] = fill(fin[2],Npop)
                u0[1+2*Npop:3*Npop] = fill(fin[3],Npop)
                u0[1+3*Npop:4*Npop] = fill(fin[4],Npop)
                u0[1+4*Npop:5*Npop] = fill(fin[5],Npop)
                u0[1+5*Npop:end] = fill(fin[6],Npop)

                # Define random vector and add to the initial condition of the network
                random_vector = 0.01*rand(Npop*Nvariables)
                u0 .+= random_vector

                # Initialize Jacobian of the network
                M = zeros(90*6,90*6)
                Network.JacobianInit!(M,u0,p,0);

                ds=ContinuousDynamicalSystem(Network.Network!, u0, p;diffeq=(dt=dt,adaptive=false))
                tands = TangentDynamicalSystem(ds; k=90,J=Network.JacobianNetwork!,J0=M)
                @time λ = lyapunovspectrum(tands, Nsteps; Δt = ns, Ttr=trans)


                # Store the SVector in the idx-th position of the matrix
                matrix[row,col] = λ

            end
        end
        # Save all data
        # Save them to a JLD2 file
        @save "eps_9_LyapunovNetwork.jld2" matrix vector_eps vectorIext_e

        return vector_eps,vectorIext_e,matrix
    end
end
