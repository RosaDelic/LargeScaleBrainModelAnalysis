# Load libraries
using ChaosTools
using DifferentialEquations
using Plots
using LinearAlgebra
using Plots.PlotMeasures
using LaTeXStrings
using JLD2
using NPZ
using Distributions
using DynamicalSystems

# Load files
include("Network.jl")
include("HomogeneousSystem.jl")

# Define a parameter struct that contains the scalar parameters and the matrix
struct Params
    scalar_params::Vector{Float64}  # Vector of scalar parameters (e.g., p[1], p[2], ...)
    matrix_params::Matrix{Float64}  # The normalized connectivity matrix
end

# Initial conditions (size should match the system, here we have a 2D system)
u0 = zeros(6)

# Time span
# Initial time
t0 = 0
# Final time
tf = 2000
tspan = (t0, tf)

# Define the parameters as a dictionary
p = [8,8,1,5,-5,-5,1,1,5,13,5,13,0,15,9];

# Define the problem using the system and parameters
prob = ODEProblem(HomogeneousSystem.HomogeneousSystem!, u0, tspan, p)
# Solve the ODE problem using RK45
sol = solve(prob,reltol = 1e-6,abstol = 1e-6,saveat=0.001);

Nvariables = 6
Npop = 90

# Set initial condition on network
u0 = zeros(Npop*Nvariables)
u0[1:Npop] = fill(fin[1],Npop)
u0[1+Npop:2*Npop] = fill(fin[2],Npop)
u0[1+2*Npop:3*Npop] = fill(fin[3],Npop)
u0[1+3*Npop:4*Npop] = fill(fin[4],Npop)
u0[1+4*Npop:5*Npop] = fill(fin[5],Npop)
u0[1+5*Npop:end] = fill(fin[6],Npop)
u0

Nvariables = 6
Npop = 90

# Define the normal distribution (mean=0, standard deviation=1)
dist = Normal(0, 0.01)

# Generate random vector
random_vector = rand(dist, Npop*Nvariables)

# Initial conditions 
u0 .+= random_vector

# Time span
# Initial time
t0 = 0
# Final time
tf = 1000
tspan = (t0, tf)

# Define the parameters as a dictionary
scalar_params = [8,8,1,5,-5,-5,1,1,5,13,5,13,0,15,9];

# Load the .npz file containing the connectivity matrix into a dictionary-like structure
data = load("NormalizedMatrix.npz")
W = data["normalized_matrix"]

# Create the Params struct
p = Params(scalar_params, W)

# Define the problem using the system and parameters
prob = ODEProblem(Network.Network!, u0, tspan, p)

# Solve the ODE problem using RK45

#@profview sol = solve(prob,reltol = 1e-6,abstol = 1e-6,saveat=0.001);
sol = solve(prob,reltol = 1e-6,abstol = 1e-6,RK4(); dt=0.01, adaptive=false);

# Unpack state variables
idx_ve = 1
idx_se = 2
idx_ri = 3
idx_vi = 4
idx_si = 5
@inbounds r_e_vector = @view sol[1:90,:]
@inbounds v_e_vector = @view sol[90*idx_ve+1:90*idx_se,:]
@inbounds s_e_vector = @view sol[90*idx_se+1:90*idx_ri,:]
@inbounds r_i_vector = @view sol[90*idx_ri+1:90*idx_vi,:]
@inbounds v_i_vector = @view sol[90*idx_vi+1:90*idx_si,:]
@inbounds s_i_vector = @view sol[90*idx_si+1:90*6,:]

# =========================   PLOT RESULTS   =======================
# Set the backend to GR
gr()  
vmin = -2
vmax = 2
idx_ms = 100
step = 0.01
# Loop over each column and plot it in a different color
plt= heatmap()
heatmap(1:Npop,sol.t[end-Int(floor(idx_ms/step)):end],transpose(v_e_vector[:,end-Int(floor(idx_ms/step)):end]), color=:RdBu,clims=(vmin, vmax),left_margin=2cm,right_margin=2cm,colorbar = true, colorbar_position = :right)

# Set labels
xlabel!("Population i ",xguidefontsize=20,bottom_margin = 10mm)
ylabel!("Time (ms)",yguidefontsize=20,left_margin = 10mm)
plot!(xtickfont=14,ytickfont=14)
# Set the legend to have two columns
#plot!(legend=:outertop, legendfontsize=5,legendcolumns=3)
# Set the figure size (make it wider)
plot!(title=L"\epsilon = " * string(p.scalar_params[15]) * L", I_{ext}^e = " * string(p.scalar_params[14]),titlefont = font(20),top_margin = 10mm)
#savefig("CheckSimsNetwork/eps=" * string(p.scalar_params[15]) * "_Iext_e=" * string(p.scalar_params[14]) * ".png")
plot!(size=(1000, 500)) 