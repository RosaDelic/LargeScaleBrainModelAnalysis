# HomogeneousSystem.jl

# Set module
module HomogeneousSystem

# Import packages
using DifferentialEquations
using StaticArrays

# Define the system
function HomogeneousSystem!(du, u, p, t)
    # du, u: Velocity and state variables vectors
    # p: Vector of parameters
    # t: Time

    # Unpack state variables
    r_e, v_e, s_e, r_i, v_i, s_i = u

    # Unpack parameters 
    tau_e, tau_i, tau_se, tau_si, nu_e, nu_i, Delta_e, Delta_i, Jee, Jei, Jii, Jie, Iext_i, Iext_e, eps = p

    #println(typeof(du))

    @inbounds du[1] = (Delta_e/(tau_e * π) + 2 * r_e * v_e) / tau_e
    @inbounds du[2] = (v_e^2 + nu_e - (π * r_e * tau_e)^2 + Iext_e + tau_e * (Jee + eps) * s_e - tau_e * Jei * s_i) / tau_e
    @inbounds du[3] = (-s_e + r_e) / tau_se
    @inbounds du[4] = (Delta_i/(tau_i * π) + 2 * r_i * v_i) / tau_i
    @inbounds du[5] = (v_i^2 + nu_i - (π * r_i * tau_i)^2 + Iext_i + tau_i * (Jie + eps) * s_e - tau_i * Jii * s_i) / tau_i
    @inbounds du[6] = (-s_i + r_i) / tau_si
    
end

# Define the jacobian of the system
function JacobianHomogeneous!(M, u, p, t)
    # M: Matrix where the Jacobian is stored
    # u: vector of state variables
    # p: Vector of parameters
    # t: time variables

    # Unpack state variables
    r_e, v_e, s_e, r_i, v_i, s_i = u

    # Unpack parameters 
    tau_e, tau_i, tau_se, tau_si, nu_e, nu_i, Delta_e, Delta_i, Jee, Jei, Jii, Jie, Iext_i, Iext_e, eps = p

    # Jacobian equations
    # r_e equation
    @inbounds M[1,1] = (2*v_e)/tau_e
    @inbounds M[1,2] = (2*r_e)/tau_e

    # v_e equation
    @inbounds M[2,1] = -2*r_e*tau_e*(π)^2
    @inbounds M[2,2] = (2*v_e)/tau_e
    @inbounds M[2,3] = Jee+eps
    @inbounds M[2,6] = -Jei

    # s_e equation
    @inbounds M[3,1] = 1/tau_se
    @inbounds M[3,3] = -1/tau_se

    # r_i equation
    @inbounds M[4,4] = (2*v_i)/tau_i
    @inbounds M[4,5] =  (2*r_i)/tau_i

    # v_i equation
    @inbounds M[5,3] = Jie+eps
    @inbounds M[5,4] = -2*r_i*tau_i*(π)^2
    @inbounds M[5,5] = (2*v_i)/tau_i
    @inbounds M[5,6] = -Jii

    # s_i equation
    @inbounds M[6,4] = 1/tau_si
    @inbounds M[6,6] = -1/tau_si

end

end