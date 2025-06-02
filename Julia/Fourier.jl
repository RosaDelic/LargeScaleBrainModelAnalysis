# Load libraries
using DifferentialEquations
using Plots
using LinearAlgebra
using Plots.PlotMeasures
using LaTeXStrings
using JLD2
using NPZ
using Distributions
using DynamicalSystems
using FFTW

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
p = [8,8,1,5,-5,-5,1,1,5,13,5,13,0,11.5,5];

# Define the problem using the system and parameters
prob = ODEProblem(HomogeneousSystem.HomogeneousSystem!, u0, tspan, p)
# Solve the ODE problem using RK45
sol = solve(prob,reltol = 1e-6,abstol = 1e-6,saveat=0.01);

# Set initial condition for the network
Nvariables = 6
Npop = 90

# Fill large-scale brain model intitial condition
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

# Generate a random vector of size 10 with values from the normal distribution
random_vector = rand(dist, Npop*Nvariables)

# Initial conditions (size should match the system, here we have a 2D system)
u0 .+= random_vector

# Time span
# Initial time
t0 = 0
# Final time
tf = 10000
tspan = (t0, tf)

# Define the parameters as a dictionary
scalar_params = [8,8,1,5,-5,-5,1,1,5,13,5,13,0,9.5,9];

# Load the .npz file containing the connectivity matrix into a dictionary-like structure
data = load("NormalizedMatrix.npz")
W = data["normalized_matrix"]

# Create the Params struct
p = Params(scalar_params, W)

# Define the problem using the system and parameters
prob = ODEProblem(Network.Network!, u0, tspan, p)

# Solve the ODE problem using RK45

#@profview sol = solve(prob,reltol = 1e-6,abstol = 1e-6,saveat=0.001);
sol = solve(prob,reltol = 1e-6,abstol = 1e-6,saveat=0.01,RK4();dt = 0.01,adaptive=false);

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

# ===========================   FIRST PLOT    ==========================
# Set the backend to GR
gr()  
vmin = -2
vmax = 2
idx_ms = 100
step = 0.01
# Loop over each column and plot it in a different color
plt= heatmap()
heatmap(1:Npop,sol.t[end-Int(floor(idx_ms/step)):end],transpose(v_e_vector[:,end-Int(floor(idx_ms/step)):end]), color=:RdBu,clims=(-2, 2),left_margin=2cm,right_margin=2cm,colorbar = true, colorbar_position = :right,dpi=500)

# Set labels
xlabel!("Population i ",xguidefontsize=30,bottom_margin = 10mm)
ylabel!("Time (ms)",yguidefontsize=30,left_margin = 10mm)

# Tick positions and numeric labels
xtick_positions = [1, 11, 21, 31, 41, 51, 61, 71, 81]
xtick_labels = [0, 10, 20, 30, 40, 50, 60, 70, 80]
ytick_positions = [9900, 9925, 9950, 9975,10000]
ytick_labels = [0,25,50,75,100]

plot!(xticks=(xtick_positions, xtick_labels),yticks=(ytick_positions, ytick_labels),xtickfont=30,ytickfont=30,tick_direction = :out,)
# Set the legend to have two columns
#plot!(legend=:outertop, legendfontsize=5,legendcolumns=3)

# Set the figure size (make it wider)
plot!(title=L"\epsilon = " * string(p.scalar_params[15]) * L", I_{ext}^e = " * string(p.scalar_params[14]),titlefont = font(30),top_margin = 10mm)
plot!(size=(1000, 575)) 
#savefig("Fourier10/NonPerturbedSimeps=" * string(p.scalar_params[15]) * "_Iext_e=" * string(p.scalar_params[14]) * ".png")
plot!(size=(1000, 575)) 


# =============================   Compute the Power Spectrum across all nodes   ==============================
NonPerturbedSimepsdelta = 0.01/1000 #Sampling rate
samFrec = 1/delta #Sampling frequency
# Eliminate offset
no_offset = v_e_vector[:,Int(1000/0.01):end].-mean(v_e_vector[:,Int(1000/0.01):end],dims=1)

# Compute Power Spectrum
matrix_PSD = M = Matrix{Float64}(undef, Npop, Int(floor(length(sol.t[Int(1000/0.01):end])/2)))

for idx in axes(matrix_PSD,1)
    #Set current signal
    current_signal = no_offset[idx,:]

    # Compute FFT
    F = fft(current_signal)
    freqs =  fftfreq(length(sol.t[Int(1000/0.01):end]), samFrec)

    # Select only half frequency spectrum
    half_idx = Int(floor(length(freqs)/2+1))
    F = F[2:half_idx]
    freqs = freqs[2:half_idx]

    # Normalize abs value FFT
    normalized = ((abs.(F)).^2)./maximum((abs.(F)).^2)

    # Save current data
    matrix_PSD[idx,:] = 10*log10.(normalized)
end

freqs =  fftfreq(length(sol.t), samFrec)
half_idx = Int(floor(length(freqs)/2+1))
freqs = freqs[2:half_idx]

# ===========================   SECOND PLOT    ==========================
# Set the backend to GR
gr()  
vmin = -2
vmax = 2
idx_ms = 100
step = 0.01
# Loop over each column and plot it in a different color
plt= heatmap()
heatmap(1:90,freqs[2:2001],transpose(matrix_PSD[:,2:2001]), color=:inferno,clims=(-60, 0),left_margin=2cm,right_margin=2cm,colorbar = true, colorbar_position = :right,dpi=500)
# Set title
#plot!(title=L"\epsilon = " * string(p[15]) * L", I_{ext}^e = " * string(p[14]),titlefont = font(20),top_margin = 10mm)
# Set labels
xlabel!("Population i ",xguidefontsize=30,bottom_margin = 10mm)
ylabel!("Frequency (Hz)",yguidefontsize=30,left_margin = 10mm)

# Tick positions and numeric labels
xtick_positions = [1, 11, 21, 31, 41, 51, 61, 71, 81]
xtick_labels = [0, 10, 20, 30, 40, 50, 60, 70, 80]
#ytick_positions = [9900, 9925, 9950, 9975,10000]
#ytick_labels = [0,25,50,75,100]
plot!(xticks=(xtick_positions, xtick_labels),xtickfont=30,ytickfont=30,tick_direction = :out)
#plot!(xtickfont=30,ytickfont=30,tick_direction = :out)
# Set the legend to have two columns
#plot!(legend=:outertop, legendfontsize=5,legendcolumns=3)
# Set the figure size (make it wider)
plot!(title=L"\epsilon = " * string(p.scalar_params[15]) * L", I_{ext}^e = " * string(p.scalar_params[14]),titlefont = font(30),top_margin = 10mm)
plot!(size=(1000, 575)) 
savefig("Fourier10/Sim_eps=" * string(p.scalar_params[15]) * "_Iext_e=" * string(p.scalar_params[14]) * ".png")
plot!(size=(1000, 575)) 


# =============================   Compute mean and Power Spectrum mean   ==============================
mean_ve = mean(v_e_vector,dims=1)

idx_ms = 2000
step = 0.01

gr()  # Set the backend to GR
plot= scatter()
# Loop over each column and plot it in a different color
plot!(sol.t[end-Int(floor(idx_ms/step)):end]/1000,mean_ve[1,end-Int(floor(idx_ms/step)):end],linewidth=2,markerstrokewidth=0,grid=false,legend=false,xticks=8:0.5:10,color="darkorange1",dpi=500)

# Set title
plot!(title=L"\epsilon = " * string(p.scalar_params[15]) * L", I_{ext}^e = " * string(p.scalar_params[14]) ,titlefont = font(30),top_margin = 10mm,left_margin=20mm,xlims=(8, 10))
# Set labels
xlabel!("Time (s)",xguidefontsize=30,bottom_margin = 10mm)
ylabel!(L"\bar{v}_E",yguidefontsize=30,left_margin = 10mm)
plot!(xtickfont=30,ytickfont=30)
# Set the legend to have two columns
#plot!(legend=:outertop, legendfontsize=5,legendcolumns=1)
# Set the figure size (make it wider)
plot!(size=(1000, 600)) 
#savefig(plot,"Fourier10/Modulation_eps" * string(p.scalar_params[15])*".png")
display(plot) 


no_offset = mean_ve[1,Int(1000/0.01):end].-mean(mean_ve[1,Int(1000/0.01):end])
F = fft(no_offset)
freqs =  fftfreq(length(sol.t), samFrec)
half_idx = Int(floor(length(freqs)/2+1))
F = F[1:half_idx]
freqs = freqs[1:half_idx]
normalized = ((abs.(F)).^2)./maximum((abs.(F)).^2)

# ===========================   THIRD PLOT    ==========================
gr()
plot = scatter()
plot!(freqs,10*log10.(normalized),linewidth=2,grid=false,legend=false,color="darkorange1",dpi=500)
plot!(xtickfont=30,ytickfont=30,xlims=(freqs[2],200),ylims=(-85,0))
plot!(title=L"\epsilon = " * string(p.scalar_params[15]) * L", I_{ext}^e = " * string(p.scalar_params[14]) ,titlefont = font(30),top_margin = 10mm,left_margin=20mm)
xlabel!("Frequency (Hz)",xguidefontsize=30,bottom_margin = 10mm)
ylabel!(L"\bar{v}_E"*" PSD (dB)",yguidefontsize=30,left_margin = 10mm)
plot!(xtickfont=30,ytickfont=30)
plot!(size=(1000, 600))
#savefig(plot,"Fourier10/PSDMean_eps" * string(p.scalar_params[15])*".png")
display(plot)


