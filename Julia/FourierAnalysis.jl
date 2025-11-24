# 0a. Import required Julia libraries

using ChaosTools
using DifferentialEquations
using LinearAlgebra
using LaTeXStrings
using JLD2
using NPZ
using Distributions
using DynamicalSystems
using FFTW
using PyPlot

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["font.family"] = "serif"
rcParams["font.serif"] = ["Computer Modern"]

# Ob. Import required codes
include("Network.jl")
include("HomogeneousSystem.jl")

# 0c. Define a parameter struct that contains the scalar parameters and the matrix
struct Params
    scalar_params::Vector{Float64}  # Vector of scalar parameters (e.g., p[1], p[2], ...)
    matrix_params::Matrix{Float64}  # The normalized connectivity matrix
end

#-------  1a. Integrate homogeneous system to find point on asymptotic regime  ------

# Initial conditions (size should match the system, here we have a 6D system)
u0 = zeros(6)

# Time span
# Initial time
t0 = 0
# Final time
tf = 10000
tspan = (t0, tf)

# Define the parameters as a dictionary
eps = 5
Iext_e = 13
p = [8,8,1,5,-5,-5,1,1,5,13,5,13,0,Iext_e,eps]

# Define the problem using the system and parameters
prob = ODEProblem(HomogeneousSystem.HomogeneousSystem!, u0, tspan, p)
# Solve the ODE problem using RK45
sol = solve(prob,reltol = 1e-6,abstol = 1e-6,saveat=0.01);
fin = sol[:,end]

# Uncouple variables
Nvariables = 6
Npop = 90

u0 = zeros(Npop*Nvariables)
u0[1:Npop] = fill(fin[1],Npop)
u0[1+Npop:2*Npop] = fill(fin[2],Npop)
u0[1+2*Npop:3*Npop] = fill(fin[3],Npop)
u0[1+3*Npop:4*Npop] = fill(fin[4],Npop)
u0[1+4*Npop:5*Npop] = fill(fin[5],Npop)
u0[1+5*Npop:end] = fill(fin[6],Npop)
u0

#-------  1b. Integrate large-scale brain model starting from homogeneous state ------

# Define the normal distribution (mean=0, standard deviation=1)
dist = Normal(0, 0.01)

# Generate a random vector of size 10 with values from the normal distribution
random_vector = rand(dist, Npop*Nvariables)

# Initial conditions (size should match the system, here we have a 2D system)
u0 .+= random_vector

# Time span
# Initial time
t0 = 0

# Define the parameters as a dictionary
scalar_params = [8,8,1,5,-5,-5,1,1,5,13,5,13,0,Iext_e,eps];

# Load the .npz file containing the connectivity matrix into a dictionary-like structure
data = load("NormalizedMatrix.npz")
W = data["normalized_matrix"]

# Create the Params struct
p = Params(scalar_params, W)

# Final time
tf = 10000
tspan = (t0, tf)
# Define the problem using the system and parameters
prob = ODEProblem(Network.Network!, u0, tspan, p)

# Solve the ODE problem using RK45
sol1 = solve(prob,reltol = 1e-6,abstol = 1e-8,saveat=0.01,RK4();dt = 0.01,adaptive=false,maxiters=1e7);

# Integrate again to eliminate transients
u0 = sol1[:,end]
# Final time
tf = 50000
tspan = (t0, tf)
# Define the problem using the system and parameters
prob = ODEProblem(Network.Network!, u0, tspan, p)
sol1 = solve(prob,reltol = 1e-6,abstol = 1e-8,saveat=0.01,RK4();dt = 0.01,adaptive=false,maxiters=1e7);
# Unpack state variables
#r_e, v_e, s_e, r_i, v_i, s_i = u
idx_ve = 1
idx_se = 2
idx_ri = 3
idx_vi = 4
idx_si = 5
temps = sol1.t
@inbounds v_e_vector = @view sol1[90*idx_ve+1:90*idx_se,:]


#-------  2a. Compute PS all nodes ------
#Sampling rate
delta = 0.01/1000 
#Sampling frequency
samFrec = 1/delta 
# Eliminate offset
no_offset = v_e_vector.-mean(v_e_vector,dims=1)

matrix_PSD = Matrix{Float64}(undef, Npop, Int(floor(length(temps)/2)))
freqs =  fftfreq(length(temps), samFrec)
half_idx = Int(floor(length(freqs)/2))
freqs = freqs[1:half_idx]

# Compute PS node by node
for idx in axes(matrix_PSD,1)
    print(idx)

    print("\n")

    current_signal = no_offset[idx,:]

    # Compute FFT
    F = fft(current_signal)
    F = F[1:half_idx]
    freqs =  fftfreq(length(temps), samFrec)
    # Normalize abs value FFT
    normalized = ((abs.(F)).^2)./maximum((abs.(F)).^2)

    # Save current data
    matrix_PSD[idx,:] = 10*log10.(normalized)
    
end

transp = transpose(matrix_PSD[:,1:Int(100/0.01)+2])


# ------ Figure PS all nodes  ------
using PyPlot

# Population indices
x = 1:90  

# Create figure and axis
fig, ax = subplots(figsize=(10, 5.75))  # 1000x575 px

# Create the heatmap
cax = ax.imshow(transp,
    aspect="auto",
    origin="lower",
    extent=[x[1], x[end], freqs[1], freqs[Int(100/0.01)+2]],
    cmap="inferno",
    vmin=-60, vmax=0
)

# Custom x-ticks
xtick_positions = [10, 30, 50, 70, 90]
xtick_labels = [10,30,50,70, 90]
ax.set_xticks(xtick_positions)
ax.set_xticklabels(xtick_labels, fontsize=50)

# Custom y-ticks
ytick_positions = [0, freqs[Int(25/0.01)], freqs[Int(50/0.01)], freqs[Int(75/0.01)], freqs[Int(100/0.01)]]
ytick_labels = [0, 50, 100, 150, 200]
ax.set_yticks(ytick_positions)
ax.set_yticklabels(ytick_labels, fontsize=50)

# y-tick font size
ax.tick_params(axis="y", labelsize=50)

# Tick direction
ax.tick_params(direction="out")

# Title (with LaTeX)
ϵ = p.scalar_params[15]
Iext = p.scalar_params[14]
# Labels and title
ax.set_xlabel("Population j", fontsize=50)
ax.set_ylabel("Frequency (Hz)", fontsize=50)
ϵ = p.scalar_params[15]
Iext = p.scalar_params[14]

title_str = L"(I_{ext}^E, \varepsilon) = (" * "$(Iext)" * L", " * "$(ϵ)" * L")"
ax.set_title(title_str, fontsize=50, pad=30)

# Colorbar
cb = fig.colorbar(cax, ax=ax)
cb.set_label(L"\mathrm{Power}\ v_E\ (\mathrm{dB})", fontsize=50)
cb.ax.tick_params(labelsize=50)
cb.ax.set_yticks([0, -30, -60])  # Assuming vertical colorbar

# Layout
tight_layout()

# Save to file
#PyPlot.savefig("SimPSD_eps=$(ϵ)_Iext_e=$(Iext)2.png", dpi=500)

# Show the plot
display(fig)

# ------ Figure vE all nodes ------
using PyPlot

# Parameters
vmin = -2
vmax = 2
idx_ms = 100
step = 0.01

# Prepare time and data slices
t_slice = sol.t[end - Int(floor(idx_ms/step)) + 1:end]
v_slice = transpose(v_e_vector[:, end - Int(floor(idx_ms/step)) + 1:end])  # shape: time x Npop

fig, ax = subplots(figsize=(10, 5.75))  # approx 1000x575 pixels

# Plot the heatmap
cax = ax.imshow(v_slice, aspect="auto", origin="lower", cmap="RdBu_r",
                extent=[1, Npop, t_slice[1], t_slice[end]], vmin=vmin, vmax=vmax)

# Colorbar
cb = fig.colorbar(cax, ax=ax)
cb.set_label(L"v_E", fontsize=50)
cb.ax.tick_params(labelsize=50)
cb.ax.set_yticks([-2,-1, 0, 1, 2])  # For vertical colorbar

xtick_positions = [10, 30, 50, 70, 90]
xtick_labels = [10,30,50,70, 90]
ytick_positions = [9900,9925,9950,9975,10000]
ytick_labels = [0,25,50,75,100]

#plot!(xticks=(xtick_positions,xtick_labels),yticks=(ytick_positions,ytick_labels),labelsize=30,xtickfont = 40,ytickfont = 40, tick_direction = :out)
# Set custom tick positions and labels
ax.set_xticks(xtick_positions)
ax.set_xticklabels(xtick_labels, fontsize=50)

ax.set_yticks(ytick_positions)
ax.set_yticklabels(ytick_labels, fontsize=50)

# Labels and title
ax.set_xlabel("Population j", fontsize=50)
ax.set_ylabel("Time (ms)", fontsize=50)
ax.set_title(L"(I_{ext}^E, \varepsilon) = ("*string(p.scalar_params[14])*L", "*string(p.scalar_params[15])*L")",
             fontsize=50, pad=30)

# Tick label size
#ax.tick_params(axis="both", labelsize=30)

# Save if needed
#PyPlot.savefig("SimPerturbed_eps=$(p.scalar_params[15])_Iext_e=$(p.scalar_params[14])2.png", dpi=500, bbox_inches="tight")

display(fig)

# ------  2b. Plot mean std of PS nodes ------
meanFourier = mean(transp,dims=2)[:]
stdFourier = std(transp,dims=2)[:]
infBoundFourier = meanFourier-stdFourier
supBoundFourier = meanFourier+stdFourier

using PyPlot

# Plot setup
fig, ax = subplots(figsize=(10, 6))
ax.plot(freqs[1:Int(100/0.01)+2], meanFourier, color="darkorange", linewidth=2)
ax.fill_between(freqs[1:Int(100/0.01)+2],infBoundFourier,supBoundFourier,color="darkorange",alpha=0.3)

# Axis limits
ax.set_xticks([0,50,100,150,200])
ax.set_xlim(0, 200)
ax.set_yticks([-100,-50,0])
ax.set_ylim(-100, 0)

# Tick font sizes
ax.tick_params(axis="both", labelsize=40)

# Labels
ax.set_xlabel("Frequency (Hz)", fontsize=40)
ax.set_ylabel(L"\mathrm{Power}\ \bar{v}_E\ (\mathrm{dB})", fontsize=40)

# Layout and save
tight_layout()
#PyPlot.savefig("PeriodicOrbintPSDMean_PO.png", dpi=500)
display(fig)

# ------  3a. Compute mean across nodes  ------

# Compute mean across nodes. --> vE_mean(t)
mean_ve = mean(v_e_vector,dims=1)


# ------ Figure time series mean vE across nodes  ------
using PyPlot

idx_ms = 2000
step = 0.01

# Compute time and data slice
t_slice = temps[end - Int(floor(idx_ms / step)) : end] ./ 1000
ve_slice = mean_ve[1, end - Int(floor(idx_ms / step)) : end]

# Create figure and plot
fig, ax = subplots(figsize=(10, 6))
ax.plot(t_slice, ve_slice, color="darkorange", linewidth=2)

# Title and labels
ϵ = p.scalar_params[15]
Iext = p.scalar_params[14]
#title_str = "\$(I_{ext}^E, \\epsilon) = ($Iext, $ϵ)\$"
#ax.set_title(title_str, fontsize=30,pad=20)
ax.set_xlabel("Time (s)", fontsize=50)
ax.set_ylabel(L"\bar{v}_E", fontsize=50)

xtick_positions = [48.0,48.5,49, 49.5,50.0]
xtick_labels = [0,0.5,1,1.5,2]
ax.set_xticks(xtick_positions)
ax.set_xticklabels(xtick_labels, fontsize=50)

# Axis limits and ticks
ax.set_xlim(48,50)
ax.tick_params(axis="both", labelsize=50)

# Layout and save
tight_layout()
#PyPlot.savefig("PeriodicOrbitModulation_eps$(ϵ)_2.png", dpi=500)


# ------  3b. Compute PS of mean vE across nodes  ------

# Eliminate offset from mean
no_offset_mean = mean_ve.-mean(mean_ve)
F = fft(no_offset_mean)
freqs =  fftfreq(length(temps), samFrec)
half_idx = Int(floor(length(freqs)/2))
F = F[1:half_idx]
freqs = freqs[1:half_idx]
normalized = ((abs.(F)).^2)./maximum((abs.(F)).^2)

# ------ Figure PS of mean vE ------
using PyPlot

# Prepare data
ϵ = p.scalar_params[15]
Iext = p.scalar_params[14]

# Compute the PSD in dB
psd_dB = 10 .* log10.(normalized)

# Plot setup
fig, ax = subplots(figsize=(10, 6))
ax.plot(freqs, psd_dB, color="darkorange", linewidth=2)

# Axis limits
ax.set_xticks([0,50,100,150,200])
ax.set_xlim(0, 200)
ax.set_yticks([-100,-50,0])
ax.set_ylim(-100, 0)
# Tick font sizes
ax.tick_params(axis="both", labelsize=50)

# Labels
ax.set_xlabel("Frequency (Hz)", fontsize=50)
ax.set_ylabel(L"\mathrm{Power}\ \bar{v}_E\ (\mathrm{dB})", fontsize=50)

# Layout and save
tight_layout()
#PyPlot.savefig("PeriodicOrbintPSDMean_eps$(ϵ)_2.png", dpi=500)
display(fig)

# ------  4. Save all  ------
npzwrite("Fourier_eps_"*string(p.scalar_params[15])*L"_Iext_e"*string(p.scalar_params[14])*L".npz",Dict("temps" => t_slice,"v_e_vector" => v_slice,"freqs" => freqs, "transp" => transp, "long_temps" => temps,"mean_ve" => mean_ve))