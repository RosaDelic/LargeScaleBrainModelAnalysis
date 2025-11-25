# 0a. Load requiered libraries in Julia
using ChaosTools
using DifferentialEquations
using PyPlot
using LinearAlgebra
using Plots.PlotMeasures
using LaTeXStrings
using JLD2
using IJulia
using DynamicalSystems
using NPZ

# Set LaTeX rendering and Computer Modern font
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["font.family"] = "serif"
rcParams["font.serif"] = ["Computer Modern"]

# 0b. Load requiered data of TransverseEqPoints
data_eqPoints = npzread("TransverseInstabilitiesVaps.npz")
dataVapsReal = data_eqPoints["dataVapsReal"]
dataVapsImaginary = data_eqPoints["dataVapsImaginary"]
maxFloquetEqPoint = data_eqPoints["maxFloquet"]
positive_directionsVaps = data_eqPoints["positive_directionsVaps"]
Vaps_2_pos = data_eqPoints["Vaps_2_pos"]
Vaps_4_pos = data_eqPoints["Vaps_4_pos"]
Vaps_6_pos = data_eqPoints["Vaps_6_pos"]
Vaps_20_pos = data_eqPoints["Vaps_20_pos"]
Vaps_90_pos = data_eqPoints["Vaps_90_pos"]

# 0c. Load requiered data of TransversePeriodicOrbits
data = npzread("TransverseInstabilitiesFloquet.npz")
vector_Iext_e = data["vector_Iext_e"]
vector_eps = data["vector_eps"]
status0 = data["status0"]
status2 = data["status2"]
maxFloquet = data["maxFloquet"]
positive_directions = data["positive_directions"]
Floquet_2_pos = data["Floquet_2_pos"]
Floquet_4_pos = data["Floquet_4_pos"]
Floquet_6_pos = data["Floquet_6_pos"]
Floquet_20_pos = data["Floquet_20_pos"]
Floquet_90_pos = data["Floquet_90_pos"]

# 0d. Load required data from BifurcationDiagramAuto
dataAuto = npzread("BifDiagramAuto.npz")
lc1_2d = dataAuto["lc1_2d"]
sd_lc_bif = dataAuto["sd_lc_bif"]
sd_lc_bif_18 = dataAuto["sd_lc_bif_18"]
sd_lc_bif_25 = dataAuto["sd_lc_bif_25"]
sn_bif_29 = dataAuto["sn_bif_29"]
pd1_bif_12 = dataAuto["pd1_bif_12"]

# pd1_bif_12
Iext_e_pd1_bif_12 = vec(pd1_bif_12[1, :, 1])
eps_pd1_bif_12 = vec(pd1_bif_12[1, :, 5])
fig = plot(Iext_e_pd1_bif_12, eps_pd1_bif_12,linewidth=3, color=:black)

# lc1_2d
Iext_e_lc1_2d = vec(lc1_2d[1, :, 1])
eps_lc1_2d = vec(lc1_2d[1, :, 5])

# sd_lc_bif
Iext_e_sd_lc_bif = vec(sd_lc_bif[1, :, 1])
eps_sd_lc_bif = vec(sd_lc_bif[1, :, 5])

# sd_lc_bif_18
Iext_e_sd_lc_bif_18 = vec(sd_lc_bif_18[1, :, 1])
eps_sd_lc_bif_18 = vec(sd_lc_bif_18[1, :, 5])

# sd_lc_bif_25
Iext_e_sd_lc_bif_25 = vec(sd_lc_bif_25[1, :, 1])
eps_sd_lc_bif_25 = vec(sd_lc_bif_25[1, :, 5])

# sn_bif_29
Iext_e_sn_bif_29 = vec(sn_bif_29[1, :, 1])
eps_sn_bif_29 = vec(sn_bif_29[1, :, 5])


using PyPlot
using PyCall
np = pyimport("numpy")  # Import NumPy

# Import matplotlib.colors
mcolors = pyimport("matplotlib.colors")
dict_color = mcolors.CSS4_COLORS 

# Create meshgrid that matches the extent
X, Y = np.meshgrid(vector_Iext_e, vector_eps)

fig, ax = subplots(figsize=(15, 7))  # Optional figure size

# Plot the heatmap for Equilibrium Points
cax2 = ax.imshow(maxFloquetEqPoint,
    origin="lower",
    extent=[minimum(vector_Iext_e), maximum(vector_Iext_e),
            minimum(vector_eps), maximum(vector_eps)],
    cmap="PiYG_r",
    vmin=-0.025, vmax=0.025,
    aspect="auto")

# Plot the heatmap for Periodic Orbits
cax = ax.imshow(maxFloquet,
                origin="lower",
                extent=[minimum(vector_Iext_e), maximum(vector_Iext_e),
                        minimum(vector_eps), maximum(vector_eps)],
                cmap="RdBu_r",
                vmin=-0.025, vmax=0.025,
                aspect="auto")

# Add colorbar
cb = fig.colorbar(cax2, ax=ax)
cb.set_label(L"\hat{\mu}_{\small{EQ}}", fontsize=30)
cb.ax.tick_params(labelsize=30)

# Axis labels
ax.set_xlabel(L"I_{ext}^E", fontsize=30)
ax.set_ylabel(L"\varepsilon", fontsize=30)

# Ticks and styling (optional)
ax.tick_params(labelsize=30)

# Adjust layout to prevent label cut-off
tight_layout()

CS2 = ax.contour(X,Y,maxFloquetEqPoint, levels=[0], colors="black", linewidths=2)

# Add contour line where maxFloquet == 0 (Periodi Orbits)
CS1 = ax.contour(X,Y,maxFloquet, levels=[0], colors="black", linewidths=2)

ax.plot(Iext_e_pd1_bif_12, eps_pd1_bif_12, linewidth=3, color=dict_color["red"], label="PD1")
ax.plot(Iext_e_lc1_2d, eps_lc1_2d, linewidth=3, color=dict_color["indigo"], label="Hopf")
ax.plot(Iext_e_sd_lc_bif_25, eps_sd_lc_bif_25, linewidth=3, color=dict_color["seagreen"], label="SN_LC_25")
ax.plot(Iext_e_sd_lc_bif, eps_sd_lc_bif, linewidth=3, color=dict_color["seagreen"], label="SN_LC")
ax.plot(Iext_e_sd_lc_bif_18, eps_sd_lc_bif_18, linewidth=3, color=dict_color["seagreen"], label="SN_LC_18")
ax.plot(Iext_e_sn_bif_29, eps_sn_bif_29, linewidth=3, color=dict_color["dodgerblue"], label="SN")

# Set axis limits
ax.set_xlim(0, 16)
ax.set_ylim(0, 30)

# Save the figure (with high DPI)
#savefig("MaxFloquetMixed2.png", dpi=500, bbox_inches="tight")
display(fig)


using PyPlot
using PyCall
np = pyimport("numpy")  

# Import matplotlib.colors
mcolors = pyimport("matplotlib.colors")
dict_color = mcolors.CSS4_COLORS 

# Create a ListedColormap with one color
ListedColormap = mcolors.ListedColormap
green_cmap = ListedColormap([dict_color["darkorange"]])

# Create meshgrid that matches the extent
X, Y = np.meshgrid(vector_Iext_e, vector_eps)

fig, ax = subplots(figsize=(15, 7))  # Optional figure size

# Plot the heatmap
cax = ax.imshow(Floquet_2_pos,
                origin="lower",
                extent=[minimum(vector_Iext_e), maximum(vector_Iext_e),
                        minimum(vector_eps), maximum(vector_eps)],
                cmap=green_cmap,
                vmin=-0.025, vmax=0.025,
                aspect="auto")

cax2 = ax.imshow(Vaps_2_pos,
                origin="lower",
                extent=[minimum(vector_Iext_e), maximum(vector_Iext_e),
                        minimum(vector_eps), maximum(vector_eps)],
                cmap=green_cmap,
                vmin=-0.025, vmax=0.025,
                aspect="auto")

# Add colorbar
cb = fig.colorbar(cax, ax=ax)
cb.set_label(L"\hat{\mu}", fontsize=30)
cb.ax.tick_params(labelsize=30)

# Axis labels
# Axis labels
ax.set_xlabel(L"I_{ext}^E", fontsize=30)
ax.set_ylabel(L"\varepsilon", fontsize=30)

# Ticks and styling (optional)
ax.tick_params(labelsize=30)

# Adjust layout to prevent label cut-off
tight_layout()

# Add contour line where maxFloquet == 0
CS = ax.contour(X,Y,maxFloquet, levels=[0], colors="black", linewidths=2)

ax.plot(Iext_e_pd1_bif_12, eps_pd1_bif_12, linewidth=3, color=dict_color["red"], label="PD1")
ax.plot(Iext_e_lc1_2d, eps_lc1_2d, linewidth=3, color=dict_color["indigo"], label="Hopf")
ax.plot(Iext_e_sd_lc_bif_25, eps_sd_lc_bif_25, linewidth=3, color=dict_color["seagreen"], label="SN_LC_25")
ax.plot(Iext_e_sd_lc_bif, eps_sd_lc_bif, linewidth=3, color=dict_color["seagreen"], label="SN_LC")
ax.plot(Iext_e_sd_lc_bif_18, eps_sd_lc_bif_18, linewidth=3, color=dict_color["seagreen"], label="SN_LC_18")
ax.plot(Iext_e_sn_bif_29, eps_sn_bif_29, linewidth=3, color=dict_color["dodgerblue"], label="SN")

# Set axis limits
ax.set_xlim(0, 16)
ax.set_ylim(0, 30)

# Save the figure (with high DPI)
#savefig("SecondAlpha.png", dpi=500, bbox_inches="tight")

using PyPlot
using PyCall
np = pyimport("numpy")  # Import NumPy

# Import matplotlib.colors
mcolors = pyimport("matplotlib.colors")
dict_color = mcolors.CSS4_COLORS 

# Create a ListedColormap with one color
ListedColormap = mcolors.ListedColormap
green_cmap = ListedColormap([dict_color["darkorange"]])

# Assume vector_Iext_e and vector_eps are 1D arrays (x and y axes)
# and maxFloquet is a 2D matrix of shape (length(vector_eps), length(vector_Iext_e))

# Create meshgrid that matches the extent
X, Y = np.meshgrid(vector_Iext_e, vector_eps)

fig, ax = subplots(figsize=(15, 7))  # Optional figure size

# Plot the heatmap
cax = ax.imshow(positive_directions,
                origin="lower",
                extent=[minimum(vector_Iext_e), maximum(vector_Iext_e),
                        minimum(vector_eps), maximum(vector_eps)],
                cmap=:plasma,
                vmin=0, vmax=89,
                aspect="auto")

cax2 = ax.imshow(positive_directions,
        origin="lower",
        extent=[minimum(vector_Iext_e), maximum(vector_Iext_e),
                minimum(vector_eps), maximum(vector_eps)],
        cmap=:plasma,
        vmin=0, vmax=89,
        aspect="auto")

# Add colorbar
cb = fig.colorbar(cax, ax=ax)
cb.set_label(L"\#"*" Unstable directions", fontsize=30)
cb.ax.tick_params(labelsize=30)

# Axis labels
ax.set_xlabel(L"I_{ext}^E", fontsize=30)
ax.set_ylabel(L"\varepsilon", fontsize=30)

# Ticks and styling (optional)
ax.tick_params(labelsize=30)

# Adjust layout to prevent label cut-off
tight_layout()

# Add contour line where maxFloquet == 0
CS = ax.contour(X,Y,maxFloquet, levels=[0], colors="black", linewidths=2)

ax.plot(Iext_e_pd1_bif_12, eps_pd1_bif_12, linewidth=3, color=dict_color["red"], label="PD1")
ax.plot(Iext_e_lc1_2d, eps_lc1_2d, linewidth=3, color=dict_color["indigo"], label="Hopf")
ax.plot(Iext_e_sd_lc_bif_25, eps_sd_lc_bif_25, linewidth=3, color=dict_color["seagreen"], label="SN_LC_25")
ax.plot(Iext_e_sd_lc_bif, eps_sd_lc_bif, linewidth=3, color=dict_color["seagreen"], label="SN_LC")
ax.plot(Iext_e_sd_lc_bif_18, eps_sd_lc_bif_18, linewidth=3, color=dict_color["seagreen"], label="SN_LC_18")
ax.plot(Iext_e_sn_bif_29, eps_sn_bif_29, linewidth=3, color=dict_color["dodgerblue"], label="SN")

# Set axis limits
ax.set_xlim(0, 16)
ax.set_ylim(0, 30)

# Save the figure (with high DPI)
#savefig("PosDirections.png", dpi=500, bbox_inches="tight")

display(fig)