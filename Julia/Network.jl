# Network.jl

# Set module
module Network

# Import packages
using StaticArrays
using LinearAlgebra

    # Define the system
    function Network!(du, u, p, t)

        # Extract parameters from p
        scalar_params = p.scalar_params
        W = p.matrix_params

        # Number of populations
        Npop = size(W)[1]

        # Define necessary indices
        idx_ve = 1
        idx_se = 2
        idx_ri = 3
        idx_vi = 4
        idx_si = 5

        # Unpack state variables
        @inbounds r_e_vector = @view u[1:Npop]
        @inbounds v_e_vector = @view u[Npop*idx_ve+1:Npop*idx_se]
        @inbounds s_e_vector = @view u[Npop*idx_se+1:Npop*idx_ri]
        @inbounds r_i_vector = @view u[Npop*idx_ri+1:Npop*idx_vi]
        @inbounds v_i_vector = @view u[Npop*idx_vi+1:Npop*idx_si]
        @inbounds s_i_vector = @view u[Npop*idx_si+1:end]
    
        # Unpack parameters 
        tau_e, tau_i, tau_se, tau_si, nu_e, nu_i, Delta_e, Delta_i, Jee, Jei, Jii, Jie, Iext_i, Iext_e, eps = scalar_params

        # Differential field model equations
        @inbounds d_re = @view du[1:Npop]
        @inbounds d_ve = @view du[Npop*idx_ve+1:Npop*idx_se]
        @inbounds d_se = @view du[Npop*idx_se+1:Npop*idx_ri]
        @inbounds d_ri = @view du[Npop*idx_ri+1:Npop*idx_vi]
        @inbounds d_vi = @view du[Npop*idx_vi+1:Npop*idx_si]
        @inbounds d_si = @view du[Npop*idx_si+1:end]

        # Compute coupling sum
        coupling_sum = W * s_e_vector

        @. d_re = (Delta_e/(tau_e * π) + 2 * r_e_vector * v_e_vector) / tau_e
        @. d_ve = (v_e_vector^2 + nu_e - (π * r_e_vector * tau_e)^2 + Iext_e + tau_e * Jee * s_e_vector + tau_e * eps * coupling_sum - tau_e * Jei * s_i_vector) / tau_e
        @. d_se = (-s_e_vector + r_e_vector) / tau_se
        @. d_ri = (Delta_i/(tau_i * π) + 2 * r_i_vector * v_i_vector) / tau_i
        @. d_vi = (v_i_vector^2 + nu_i - (π * r_i_vector * tau_i)^2 + Iext_i + tau_i * Jie * s_e_vector + tau_i *  eps * coupling_sum - tau_i * Jii * s_i_vector) / tau_i
        @. d_si = (-s_i_vector + r_i_vector) / tau_si
        
    end

    function JacobianNetwork!(M, u, p, t)
        # M: Matrix where the Jacobian is stored
        # u: vector of state variables
        # p: Vector of parameters
        # t: time variables

        # Extract parameters from p
        scalar_params = p.scalar_params
        W = p.matrix_params

        # Number of populations
        Npop = 90

        # Number of variables
        Nvar = 6

        # Define necessary indices        
        idx_ve = 1
        idx_se = 2
        idx_ri = 3
        idx_vi = 4
        idx_si = 5

        # Unpack state variables
        @inbounds r_e_vector = @view u[1:Npop]
        @inbounds v_e_vector = @view u[Npop*idx_ve+1:Npop*idx_se]

        @inbounds r_i_vector = @view u[Npop*idx_ri+1:Npop*idx_vi]
        @inbounds v_i_vector = @view u[Npop*idx_vi+1:Npop*idx_si]


        # Unpack parameters 
        tau_e, tau_i, tau_se, tau_si, nu_e, nu_i, Delta_e, Delta_i, Jee, Jei, Jii, Jie, Iext_i, Iext_e, eps = scalar_params

        # For the main diagonal blocks
        d0 = diagind(M);
        @inbounds r_e_diag = view(M,d0[1:Npop])
        @inbounds v_e_diag = view(M,d0[1+Npop:idx_se*Npop])
        @inbounds s_e_diag = view(M,d0[1+idx_se*Npop:idx_ri*Npop])
        @inbounds r_i_diag = view(M,d0[1+idx_ri*Npop:idx_vi*Npop])
        @inbounds v_i_diag = view(M,d0[1+idx_vi*Npop:idx_si*Npop])
        @inbounds s_i_diag = view(M,d0[1+idx_si*Npop:end])

        # For the first upper diagonal blocks
        d1 = diagind(M,Npop);
        @inbounds r_e_UpDiag = view(M,d1[1:Npop])
        @inbounds v_e_UpDiag = view(M,d1[1+Npop:idx_se*Npop])
        @inbounds r_i_Updiag = view(M,d1[1+idx_ri*Npop:idx_vi*Npop])
        @inbounds v_i_Updiag = view(M,d1[1+idx_vi*Npop:end])

        # For the first lower diagonal blocks
        ds1 = diagind(M,-Npop);
        @inbounds v_e_DownDiag = view(M,ds1[1:Npop])
        @inbounds v_i_DownDiag = view(M,ds1[1+idx_ri*Npop:idx_vi*Npop])

        # For the second lower diagonal blocks
        ds2 = diagind(M,-idx_se*Npop);
        @inbounds s_e_SecondDownDiag = view(M,ds2[1:Npop])
        @inbounds v_i_SecondDownDiag = view(M,ds2[1+idx_se*Npop:idx_ri*Npop])
        @inbounds s_i_SecondDownDiag = view(M,ds2[1+idx_ri*Npop:end])

        # For the left most diagonal block
        ds_left = diagind(M,idx_vi*Npop);
        @inbounds v_e_left_most = view(M,ds_left[1+Npop:end])

        # ----------------  Partial derivatives of r_e  ----------------
        # Partial r_e
        @. r_e_diag = (2/tau_e)*v_e_vector
        # Partial v_e
        @. r_e_UpDiag = (2/tau_e)*r_e_vector

        # ----------------  Partial derivatives of v_e  ----------------
        # Partial r_e
        @. v_e_DownDiag = -2*r_e_vector*tau_e*(π)^2
        # Partial v_e
        @. v_e_diag = (2/tau_e)*v_e_vector

        # ----------------  Partial derivatives of r_i  ----------------
        # Partial r_i
        @. r_i_diag = (2/tau_i)*v_i_vector
        # Partial v_i
        @. r_i_Updiag = (2/tau_i)*r_i_vector
        # Partial s_e, r_i, v_i, s_i --> All zeros

        # ----------------  Partial derivatives of v_i  ----------------

        # Partial r_i
        @. v_i_DownDiag = -2*r_i_vector*tau_i*(π)^2
        # Partial v_i
        @. v_i_diag = (2/tau_i)*v_i_vector




    end

    function JacobianInit!(M, u, p, t)
        # M: Matrix where the Jacobian is stored
        # u: vector of state variables
        # p: Vector of parameters
        # t: time variables

        # Extract parameters from p
        scalar_params = p.scalar_params
        W = p.matrix_params

        # Number of populations
        Npop = 90

        # Number of variables
        Nvar = 6

        # Define necessary indices        
        idx_ve = 1
        idx_se = 2
        idx_ri = 3
        idx_vi = 4
        idx_si = 5

        # Unpack state variables
        @inbounds r_e_vector = @view u[1:Npop]
        @inbounds v_e_vector = @view u[Npop*idx_ve+1:Npop*idx_se]

        @inbounds r_i_vector = @view u[Npop*idx_ri+1:Npop*idx_vi]
        @inbounds v_i_vector = @view u[Npop*idx_vi+1:Npop*idx_si]


        # Unpack parameters 
        tau_e, tau_i, tau_se, tau_si, nu_e, nu_i, Delta_e, Delta_i, Jee, Jei, Jii, Jie, Iext_i, Iext_e, eps = scalar_params

        # For the main diagonal blocks
        d0 = diagind(M);
        @inbounds r_e_diag = view(M,d0[1:Npop])
        @inbounds v_e_diag = view(M,d0[1+Npop:idx_se*Npop])
        @inbounds s_e_diag = view(M,d0[1+idx_se*Npop:idx_ri*Npop])
        @inbounds r_i_diag = view(M,d0[1+idx_ri*Npop:idx_vi*Npop])
        @inbounds v_i_diag = view(M,d0[1+idx_vi*Npop:idx_si*Npop])
        @inbounds s_i_diag = view(M,d0[1+idx_si*Npop:end])

        # For the first upper diagonal blocks
        d1 = diagind(M,Npop);
        @inbounds r_e_UpDiag = view(M,d1[1:Npop])
        @inbounds v_e_UpDiag = view(M,d1[1+Npop:idx_se*Npop])
        @inbounds r_i_Updiag = view(M,d1[1+idx_ri*Npop:idx_vi*Npop])
        @inbounds v_i_Updiag = view(M,d1[1+idx_vi*Npop:end])

        # For the first lower diagonal blocks
        ds1 = diagind(M,-Npop);
        @inbounds v_e_DownDiag = view(M,ds1[1:Npop])
        @inbounds v_i_DownDiag = view(M,ds1[1+idx_ri*Npop:idx_vi*Npop])

        # For the second lower diagonal blocks
        ds2 = diagind(M,-idx_se*Npop);
        @inbounds s_e_SecondDownDiag = view(M,ds2[1:Npop])
        @inbounds v_i_SecondDownDiag = view(M,ds2[1+idx_se*Npop:idx_ri*Npop])
        @inbounds s_i_SecondDownDiag = view(M,ds2[1+idx_ri*Npop:end])

        # For the left most diagonal block
        ds_left = diagind(M,idx_vi*Npop);
        @inbounds v_e_left_most = view(M,ds_left[1+Npop:end])

        # ----------------  Partial derivatives of r_e  ----------------
        # Partial r_e
        @. r_e_diag = (2/tau_e)*v_e_vector
        # Partial v_e
        @. r_e_UpDiag = (2/tau_e)*r_e_vector
        # Partial s_e, r_i, v_i, s_i --> All zeros

        # ----------------  Partial derivatives of v_e  ----------------
        # Partial r_e
        @. v_e_DownDiag = -2*r_e_vector*tau_e*(π)^2
        # Partial v_e
        @. v_e_diag = (2/tau_e)*v_e_vector
        # Partial s_e
        v_e_UpDiag .= fill(Jee,Npop)
        @inbounds @. M[1+Npop:idx_se*Npop,1+idx_se*Npop:idx_ri*Npop] += eps * W
        # Partial s_i
        v_e_left_most .= fill(-Jei,Npop)
        # Partial r_i, v_i --> All zeros

        # ----------------  Partial derivatives of s_e  ----------------
        # Partial r_e
        s_e_SecondDownDiag .= fill(1/tau_se,Npop)
        # Partial s_e
        s_e_diag .= fill(-1/tau_se,Npop)
        # Partial v_e, r_i, v_i, s_i --> All zeros

        # ----------------  Partial derivatives of r_i  ----------------
        # Partial r_i
        @. r_i_diag = (2/tau_i)*v_i_vector
        # Partial v_i
        @. r_i_Updiag = (2/tau_i)*r_i_vector
        # Partial s_e, r_i, v_i, s_i --> All zeros

        # ----------------  Partial derivatives of v_i  ----------------
        # Partial s_e
        v_i_SecondDownDiag .= fill(Jie,Npop)
        @inbounds @. M[1+idx_vi*Npop:idx_si*Npop,1+idx_se*Npop:idx_ri*Npop] += eps * W
        
        # Partial r_i
        @. v_i_DownDiag = -2*r_i_vector*tau_i*(π)^2
        # Partial v_i
        @. v_i_diag = (2/tau_i)*v_i_vector
        # Partial s_i
        v_i_Updiag .= fill(-Jii,Npop)
        # Partial r_e, v_e --> All zeros

        # ----------------  Partial derivatives of s_i  ----------------
        # Partial r_i
        s_i_SecondDownDiag .= fill(1/tau_si,Npop)
        # Partial s_i
        s_i_diag .= fill(-1/tau_si,Npop)
        # Partial v_e, r_i, v_i, s_i --> All zeros

    end
end
