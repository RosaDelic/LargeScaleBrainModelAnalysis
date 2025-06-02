# Load libraries
using ChaosTools
using DifferentialEquations
using Plots
using LinearAlgebra
using Plots.PlotMeasures
using LaTeXStrings
using DynamicalSystems

# Load files
include("HomogeneousSystem.jl")

# Initial conditions (size should match the system, here we have a 2D system)
u0 = zeros(6)

# Time span
# Initial time
t0 = 0
# Final time
tf = 200
tspan = (t0, tf)

# Define the parameters as a dictionary
p = [8,8,1,5,-5,-5,1,1,5,13,5,13,0,12,8];

# Define the problem using the system and parameters
prob = ODEProblem(HomogeneousSystem.HomogeneousSystem!, u0, tspan, p)

# Solve the ODE problem using RK45
sol = solve(prob,reltol = 1e-6,abstol = 1e-6,saveat=0.001);


fig = plot(sol.t, [sol[4,:] sol[1,:]], label=["r_i" "r_e"],linewidth=3, color=[:green :blue])

# Customize labels, ticks, and size
xlabel!("Time (ms)",xguidefontsize=20,bottom_margin = 10mm)
ylabel!("Voltage",yguidefontsize=20,left_margin = 10mm)
plot!(xtickfont=14,ytickfont=14)
# Set the legend to have two columns
plot!(legend=:topright, legendfontsize=12,legendcolumns=2)
# Set the figure size (make it wider)
plot!(size=(800, 400)) 
# Adjust the margins to prevent labels from being cut off
#plot!(margin=100)  # This will increase the margins around the plot
# Set the plot title with LaTeX formatting
plot!(title=L"\epsilon = " * string(p[15]) * L", I_{ext}^e = " * string(p[14]),titlefont = font(20),top_margin = 10mm)


fig = plot(sol.t, [sol[5,:] sol[2,:]], label=["v_i" "v_e"],linewidth=3, color=[:magenta :red])
# Customize labels, ticks, and size
xlabel!("Time (ms)",xguidefontsize=20,bottom_margin = 10mm)
ylabel!("Voltage",yguidefontsize=20,left_margin = 10mm)
plot!(xtickfont=14,ytickfont=14)
# Set the legend to have two columns
plot!(legend=:topright, legendfontsize=12,legendcolumns=2)
# Set the figure size (make it wider)
plot!(size=(800, 400)) 
# Adjust the margins to prevent labels from being cut off
#plot!(margin=100)  # This will increase the margins around the plot
# Set the plot title with LaTeX formatting
plot!(title=L"\epsilon = " * string(p[15]) * L", I_{ext}^e = " * string(p[14]),titlefont = font(20),top_margin = 10mm)


fig = plot(sol.t, [sol[6,:] sol[3,:]], label=["s_i" "s_e"],linewidth=3, color=[:yellow :orange])
# Customize labels, ticks, and size
xlabel!("Time (ms)",xguidefontsize=20,bottom_margin = 10mm)
ylabel!("Synapsis",yguidefontsize=20,left_margin = 10mm)
plot!(xtickfont=14,ytickfont=14)
# Set the legend to have two columns
plot!(legend=:topright, legendfontsize=12,legendcolumns=2)
# Set the figure size (make it wider)
plot!(size=(800, 400)) 
# Adjust the margins to prevent labels from being cut off
#plot!(margin=100)  # This will increase the margins around the plot
# Set the plot title with LaTeX formatting
plot!(title=L"\epsilon = " * string(p[15]) * L", I_{ext}^e = " * string(p[14]),titlefont = font(20),top_margin = 10mm)
