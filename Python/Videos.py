# 0a. Import required Python libraries
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation

# Enable LaTeX rendering
plt.rcParams.update({
    "text.usetex": True,  # Use LaTeX for all text
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
})

# 0b. Sort the spatial distribution of the nodes
# Load .dat file — adjust the separator (e.g. whitespace, comma, etc.)
df = pd.read_csv('aal.dat', delim_whitespace=True, header=None)  # or use sep=',' or sep='\t'

# Assume the label column is the last one
label_col = 6

# Split the dataframe into 'L' and 'R' groups
df_R = df[df[label_col] == 'R'].iloc[::-1]  # reverse the order
df_L = df[df[label_col] == 'L']            # keep as is

# Concatenate back: L first (reversed), then R
df_sorted = pd.concat([df_L, df_R], ignore_index=True)

# Relabel index column
df_sorted.iloc[:, 0] = range(len(df_sorted))

# Save and load Sorted spatial distribution
df_sorted.to_csv('Sorted_aal.dat', sep='\t', index=False, header=False)
SpatialDistribSort = np.genfromtxt('Sorted_aal.dat', delimiter='\t', encoding=None,usecols=(0, 2, 3, 4,7))
# Load Normalized Structural Connectivity matrix
data=np.load('NormalizedMatrix.npz')
norm_matrix = data['normalized_matrix']
vaps,veps = np.linalg.eig(norm_matrix)

# 1. Make movie

#Initial and final times of integration
t0 = 0
tf = 1000
#Discretization used for integration
h = 0.01
#Number of points to evaluate the time integration
N = int((tf-t0)/h)
time_eval = np.linspace(t0, tf, N)
#Number of variables and equations of the model
Nvariables = 6 #Number of variables of each population model
Npop = 90
Neq = Npop*Nvariables #Number of variables of the whole network

#Dictionary containing parameters of the model
ParametersPop = dict(tau_e = 8,
              tau_i = 8,
              tau_se=1,
              tau_si=5,
              nu_e = -5,
              nu_i = -5,
              Delta_e = 1,
              Delta_i = 1,
              Jee = 5,
              Jei = 13,
              Jii = 5,
              Jie = 13,
              Iext_i=0,
              Iext_e=13,
              eps=5)

Sim = np.load("Fourier_eps_13_Iext_e_5.npz")
temps = Sim["temps"]
v_e_matrix = Sim["v_e_vector"]

# Setup figure and axes
fig = plt.figure(figsize=(15, 5))
gs = GridSpec(1, 4, width_ratios=[1, 1, 1, 0.05], figure=fig)
view_angles = [(90, 90), (30, 60), (0, 90)]
point_size = 600

colors = np.zeros(Npop)  # Initial colors
sc_list = []             # Store scatter objects for animation

for i in range(3):
    ax = fig.add_subplot(gs[i], projection='3d')
    sc = ax.scatter(
        SpatialDistribSort[:, 1], SpatialDistribSort[:, 2], SpatialDistribSort[:, 3],
        c=colors, s=point_size, cmap='RdBu_r', vmin=-2, vmax=2
    )
    ax.view_init(*view_angles[i])

    # Turn off axes cleanly
    ax.grid(False)
    ax.set_axis_off()

    if i==1:
        # Set title
        ax.set_title(r'$(I_{ext}^E,\varepsilon) = $(' + str(ParametersPop['Iext_e']) + r', $' + str(ParametersPop['eps'])+ r')$',fontsize=30,fontname='Times New Roman')

    # Add text at the bottom center
    time_text = fig.text(0.5, 0.1, ' ',ha='center', va='center', fontsize=30, fontname='Times New Roman')
    sc_list.append(sc)

# Add colorbar
cbar_ax = fig.add_subplot(gs[3])
cbar = fig.colorbar(sc_list[0], cax=cbar_ax)
cbar.ax.tick_params(labelsize=30)
cbar.set_label(r'$v_E$', fontsize=30, fontname='Times New Roman')

# Init function
def init():
    time_text.set_text('t = 0') 
    for sc in sc_list:
        sc.set_array(np.zeros(Npop))  # Set to initial zero state
    return sc_list + [time_text]

# Animation update function
def animate(i):
    ms_idx = 140  # time window
    time_idx = len(temps) - int(ms_idx / h) + i-1
    color_frame = v_e_matrix[time_idx,:]

    for sc in sc_list:
        sc.set_array(color_frame)

    # Update the time text
    time_text.set_text(f't = {i*h:.0f} ms')

    return sc_list+ [time_text]

# Create animation
ms_idx = 140
frame_count = int(ms_idx / h)
anim = FuncAnimation(fig, animate, init_func=init,
                     frames=np.arange(int(140/h)+10,step=10), interval=1, blit=False)

# Save
anim.save('Animation_eps_'+str(ParametersPop['eps'])+
          '_Iext_e_'+str(ParametersPop['Iext_e'])+'.mp4',
          fps=10, extra_args=['-vcodec', 'libx264'])

plt.show()