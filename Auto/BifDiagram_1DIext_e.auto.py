# Initializing and exploring one-parameter line bifurcations from a fixed point:
# 1. Use Euler integrator up to steady state to initialize the diagram
eps = 0
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':-5.0, 'eps' : eps})
#2. Continue along p an uncover the bifurcations (2 saddle-node and 3 Hopf)         
ic=init(201)
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#save(fp,'fp') 

# 2. Limit cycles:
# 2a. Continue the limit cycle solution (IPS=2) emerging from the first Hopf.
# This branch dies on a SNIC bifurcation, characterized by a infinite period.
# Therefore, we stop the bifurcation at PERIOD = 100.0 (it can be continued further if needed)
hb=fp('HB1')
lc1=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})
#   =======================   PROBLEMS   =============================
lc1=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={'PERIOD' : 300})
hom1 = run(lc1('UZ'),ISW=1,ICP=[1,2,11],ISP=0)
hom2 = run(lc1('UZ'),ISW=1,ICP=[1,2,11],ISP=0,DS="-")
hom_35 = merge(hom1_35+hom2_35)
hom_35=relabel(hom_35)

if eps<=50:
    hb=fp('HB2')
    lc2=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})
if eps >=80 and eps <85:
    hb=fp('HB3')
    lc3=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

# 2b. The second Hopf produces a limit cycle that vanishes at another Hopf (which auto labels as 'BP' when continuing limit-cycles)
#hb=fp('HB2')
#lc2=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, STOP=['BP1'])
#save(lc2,'lc2')


# =================================   PLOTS   =================================

eps = 12
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':-5.0, 'eps' : eps})
#2. Continue along p an uncover the bifurcations (2 saddle-node and 3 Hopf)
ic=init(201)
fp_12=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#3. Limit cycle
hb_12=fp_12('HB1')
lc1_12=run(hb_12,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

eps = 15
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':-5.0, 'eps' : eps})
#2. Continue along p an uncover the bifurcations (2 saddle-node and 3 Hopf)
ic=init(201)
fp_15=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#3. Limit cycle
hb_15=fp_15('HB1')
lc1_15=run(hb_15,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

eps = 19
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':-5.0, 'eps' : eps})
#2. Continue along p an uncover the bifurcations (2 saddle-node and 3 Hopf)
ic=init(201)
fp_19=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#3. Limit cycle
hb_19=fp_19('HB1')
lc1_19=run(hb_19,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

# =========================================================================

eps = 22
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':-5.0, 'eps' : eps})
#2. Continue along p an uncover the bifurcations (2 saddle-node and 3 Hopf)
ic=init(201)
fp_22=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#3. Limit cycle
hb_22=fp_22('HB1')
lc1_22=run(hb_22,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={'PERIOD' : -10})

eps = 26
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':-5.0, 'eps' : eps})
#2. Continue along p an uncover the bifurcations (2 saddle-node and 3 Hopf)
ic=init(201)
fp_26=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#3. Limit cycle
hb_26=fp_26('HB1')
lc1_26=run(hb_26,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

eps = 29.5
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':-5.0, 'eps' : eps})
#2. Continue along p an uncover the bifurcations (2 saddle-node and 3 Hopf)
ic=init(201)
fp_29=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#3. Limit cycle
hb_29=fp_29('HB1')
lc1_29=run(hb_29,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

eps = 30
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':-5.0, 'eps' : eps})
#2. Continue along p an uncover the bifurcations (2 saddle-node and 3 Hopf)
ic=init(201)
fp_30=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#3. Limit cycle
hb_30=fp_30('HB1')
lc1_30=run(hb_30,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})



# 3. Plot results
# 3a. Some imports
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import BoundaryNorm
from matplotlib.transforms import Bbox
import numpy as np

# 3b. Auxiliary functions
def pt_vals(f):
	return np.array([f[0][i]['PT'] for i in range(len(f[0]))])

def bifs(f,par):
	exceptions = ['No Label', 'RG', 'EP', 'UZ','LP']  # List of exceptions
	return [[f[0][i]['TY name'],f[0][i][par]] for i in range(len(f[0])) if f[0][i]['TY name'] not in exceptions]

col_map = colormaps['Purples']
boundary = BoundaryNorm([0,col_map.N],col_map.N)
Ntons = col_map.N

# Figure
fig1=plt.figure(figsize=(15,7))
ax=plt.axes()
plt.title('$\epsilon =$ '+str(eps),fontsize=50,fontname='Times New Roman')
plt.xlabel('Parameter $I_{ext}^E$',fontsize=50,fontname='Times New Roman')
plt.ylabel('$v_E$ (mV)',fontsize=50,fontname='Times New Roman')
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
# General case

plt.scatter(lc1['Iext_e'],lc1['y_max'], c=Ntons * (pt_vals(lc1) < 0) + 2 * (pt_vals(lc1) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(lc1['Iext_e'],lc1['y_min'], c=Ntons * (pt_vals(lc1) < 0) + 2 * (pt_vals(lc1) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(fp['Iext_e'],fp['v_e'], c=Ntons * (pt_vals(fp) < 0) + 2 * (pt_vals(fp) > 0), cmap='PuRd', s=12,norm=boundary)
# eps = 12
plt.scatter(lc1_12['Iext_e'],lc1_12['y_max'], c=Ntons * (pt_vals(lc1_12) < 0) + 2 * (pt_vals(lc1_12) > 0), cmap='Blues', s=12,norm=boundary)
plt.scatter(lc1_12['Iext_e'],lc1_12['y_min'], c=Ntons * (pt_vals(lc1_12) < 0) + 2 * (pt_vals(lc1_12) > 0), cmap='Blues', s=12,norm=boundary)
# eps = 15
plt.scatter(lc1_15['Iext_e'],lc1_15['y_max'], c=Ntons * (pt_vals(lc1_15) < 0) + 2 * (pt_vals(lc1_15) > 0), cmap='Greens', s=12,norm=boundary)
plt.scatter(lc1_15['Iext_e'],lc1_15['y_min'], c=Ntons * (pt_vals(lc1_15) < 0) + 2 * (pt_vals(lc1_15) > 0), cmap='Greens', s=12,norm=boundary)
# eps = 19
plt.scatter(lc1_19['Iext_e'],lc1_19['y_max'], c=Ntons * (pt_vals(lc1_19) < 0) + 2 * (pt_vals(lc1_19) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(lc1_19['Iext_e'],lc1_19['y_min'], c=Ntons * (pt_vals(lc1_19) < 0) + 2 * (pt_vals(lc1_19) > 0), cmap='Purples', s=12,norm=boundary)

# eps = 22
plt.scatter(lc1_22['Iext_e'],lc1_22['y_max'], c=Ntons * (pt_vals(lc1_22) < 0) + 2 * (pt_vals(lc1_22) > 0), cmap='Blues', s=12,norm=boundary)
plt.scatter(lc1_22['Iext_e'],lc1_22['y_min'], c=Ntons * (pt_vals(lc1_22) < 0) + 2 * (pt_vals(lc1_22) > 0), cmap='Blues', s=12,norm=boundary)
#eps = 26
plt.scatter(lc1_26['Iext_e'],lc1_26['y_max'], c=Ntons * (pt_vals(lc1_26) < 0) + 2 * (pt_vals(lc1_26) > 0), cmap='Greens', s=12,norm=boundary)
plt.scatter(lc1_26['Iext_e'],lc1_26['y_min'], c=Ntons * (pt_vals(lc1_26) < 0) + 2 * (pt_vals(lc1_26) > 0), cmap='Greens', s=12,norm=boundary)
# eps = 29
plt.scatter(lc1_29['Iext_e'],lc1_29['y_max'], c=Ntons * (pt_vals(lc1_29) < 0) + 2 * (pt_vals(lc1_29) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(lc1_29['Iext_e'],lc1_29['y_min'], c=Ntons * (pt_vals(lc1_29) < 0) + 2 * (pt_vals(lc1_29) > 0), cmap='Purples', s=12,norm=boundary)
# eps = 30
plt.scatter(lc1_30['Iext_e'],lc1_30['y_max'], c=Ntons * (pt_vals(lc1_30) < 0) + 2 * (pt_vals(lc1_30) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(lc1_30['Iext_e'],lc1_30['y_min'], c=Ntons * (pt_vals(lc1_30) < 0) + 2 * (pt_vals(lc1_30) > 0), cmap='Purples', s=12,norm=boundary)

if eps<=50:
    plt.scatter(lc2['Iext_e'],lc2['y_min'], c=Ntons * (pt_vals(lc2) < 0) + 2 * (pt_vals(lc2) > 0), cmap='Purples', s=12,norm=boundary)
    plt.scatter(lc2['Iext_e'],lc2['y_max'], c=Ntons * (pt_vals(lc2) < 0) + 2 * (pt_vals(lc2) > 0), cmap='Purples', s=12,norm=boundary)
if eps>=80 and eps<85:
    plt.scatter(lc3['Iext_e'],lc3['y_min'], c=Ntons * (pt_vals(lc3) < 0) + 2 * (pt_vals(lc3) > 0), cmap='Purples', s=12,norm=boundary)
    plt.scatter(lc3['Iext_e'],lc3['y_max'], c=Ntons * (pt_vals(lc3) < 0) + 2 * (pt_vals(lc3) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(fp_29['Iext_e'],fp_29['v_e'], c=Ntons * (pt_vals(fp_29) < 0) + 2 * (pt_vals(fp_29) > 0), cmap='PuRd', s=12,norm=boundary)
#plt.scatter(lc2['Iext_e'],lc2['y_max'], c=Ntons * (pt_vals(lc2) < 0) + 2 * (pt_vals(lc2) > 0), cmap='Purples', s=12,norm=boundary)
#plt.scatter(lc2['Iext_e'],lc2['y_min'], c=Ntons * (pt_vals(lc2) < 0) + 2 * (pt_vals(lc2) > 0), cmap='Purples', s=12,norm=boundary)

#bfp = bifs(fp,'Iext_e')
#for b in bfp:
#	plt.axvline(b[1],color="black", ls="--",alpha=0.7)
#bfp = bifs(lc1,'Iext_e')
#for b in bfp:
#    plt.axvline(b[1],color="red", ls="--",alpha=0.7)
plt.xlim([0,70])
plt.savefig('Definitive/ZoomedBifDiagramEps__'+str(eps)+'.png', dpi=600,bbox_inches=Bbox([[0,-1],fig1.get_size_inches()]))
#plt.show()

fig1=plt.figure(figsize=(15,7))
ax=plt.axes()
plt.title('$\epsilon =$ '+str(eps),fontsize=30,fontname='Times New Roman')
plt.xlabel('Parameter $I_{ext}^e$',fontsize=30,fontname='Times New Roman')
plt.ylabel('$Frequency$ (Hz)',fontsize=30,fontname='Times New Roman')
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.scatter(lc1_12['Iext_e'],1000/lc1_12['PERIOD'], c=Ntons * (pt_vals(lc1_12) < 0) + 2 * (pt_vals(lc1_12) > 0), cmap='Blues', s=12,norm=boundary)
plt.scatter(lc1_15['Iext_e'],1000/lc1_15['PERIOD'], c=Ntons * (pt_vals(lc1_15) < 0) + 2 * (pt_vals(lc1_15) > 0), cmap='Greens', s=12,norm=boundary)
plt.scatter(lc1_19['Iext_e'],1000/lc1_19['PERIOD'], c=Ntons * (pt_vals(lc1_19) < 0) + 2 * (pt_vals(lc1_19) > 0), cmap='Purples', s=12,norm=boundary)

plt.scatter(lc1_22['Iext_e'],1000/lc1_22['PERIOD'], c=Ntons * (pt_vals(lc1_22) < 0) + 2 * (pt_vals(lc1_22) > 0), cmap='Blues', s=12,norm=boundary)
plt.scatter(lc1_26['Iext_e'],1000/lc1_26['PERIOD'], c=Ntons * (pt_vals(lc1_26) < 0) + 2 * (pt_vals(lc1_26) > 0), cmap='Greens', s=12,norm=boundary)
plt.scatter(lc1_29['Iext_e'],1000/lc1_29['PERIOD'], c=Ntons * (pt_vals(lc1_29) < 0) + 2 * (pt_vals(lc1_29) > 0), cmap='Purples', s=12,norm=boundary)
if eps<=50:
    plt.scatter(lc2['Iext_e'],1000/lc2['PERIOD'], c=Ntons * (pt_vals(lc2) < 0) + 2 * (pt_vals(lc2) > 0), cmap='Greys', s=12,norm=boundary)
if eps>=80 and eps<85:
    plt.scatter(lc3['Iext_e'],1000/lc3['PERIOD'], c=Ntons * (pt_vals(lc3) < 0) + 2 * (pt_vals(lc3) > 0), cmap='Greys', s=12,norm=boundary)
#plt.scatter(lc2['Iext_e'],1000/lc2['PERIOD'], c=Ntons * (pt_vals(lc2) < 0) + 2 * (pt_vals(lc2) > 0), cmap='Greys', s=12,norm=boundary)
#bfp = bifs(fp,'Iext_e')
#for b in bfp:
#    plt.axvline(b[1],color="black", ls="--",alpha=0.7)
#bfp = bifs(lc1,'Iext_e')
#for b in bfp:
#    plt.axvline(b[1],color="red", ls="--",alpha=0.7)
plt.xlim([0,15])
plt.ylim([0,70])
#plt.savefig('Definitive/FrequencyEps2__'+str(eps)+'.png', dpi=600,bbox_inches=Bbox([[0,-1],fig1.get_size_inches()]))
plt.show()


