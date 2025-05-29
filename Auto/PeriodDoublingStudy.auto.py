# -----------------------------  eps = 12  ---------------------------
eps = 12
# 1. Use Euler integrator up to steady state to initialize the diagram
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':0.0, 'eps' : eps})

#2. Continue along p an uncover the bifurcations (2 Hopf)
ic=init(201)
fp_12=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})

# 3. Continue the limit cycle solutions (ISW=1,IPS=2) emerging from the first Hopf.
hb_12=fp_12('HB1')
lc1_12=run(hb_12,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

# 4. Continue PD emerging from the first Hopf
pd1_12=lc1_12('PD1')
pd1_bif_12=run(pd1_12,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
pd1_bif_12=run(pd1_bif_12('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})

# 5. Branch switch at Second Period Doubling (the beginning) label
# When we do this no more Period Doublings appear
bsw_12=lc1_12('PD1')
bsw_bif_12=run(bsw_12,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=-1, UZSTOP={})

# 6. Continue First PD emerging from the Hopf in 2D
pd2_12=bsw_bif_12('PD1')
pd2_bif_12=run(pd2_12,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
pd2_bif_12=run(pd2_bif_12('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})

# 7. Branch switch at Second Period Doubling (the end) label
# When we do this more Period Doubling appear
SecondBsw_12=bsw_bif_12('PD1')
SecondBsw_bif_12=run(SecondBsw_12,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=-1, UZSTOP={})

# 8. Continue Second PD emerging from the Hopf in 2D
pd3_12=SecondBsw_bif_12('PD1')
pd3_bif_12=run(pd3_12,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
pd3_bif_12=run(pd3_bif_12('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})

# 9. Branch switch at Third Period Doubling (the end) label
# When we do this more Period Doubling appear
ThirdBsw_12=SecondBsw_bif_12('PD1')
ThirdBsw_bif_12=run(ThirdBsw_12,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=-1, UZSTOP={})

# 10. Continue Third PD emerging from the Hopf in 2D
pd4_12=ThirdBsw_bif_12('PD1')
pd4_bif_12=run(pd4_12,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
pd4_bif_12=run(pd4_bif_12('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})

# -----------------------------  eps = 14  ---------------------------
eps = 14
# 1. Use Euler integrator up to steady state to initialize the diagram
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':0.0, 'eps' : eps})

#2. Continue along p an uncover the bifurcations (2 Hopf)
ic=init(201)
fp_14=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})

# 3. Continue the limit cycle solutions (ISW=1,IPS=2) emerging from the first Hopf.
hb_14=fp_14('HB1')
lc1_14=run(hb_14,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

# 4. Continue PD emerging from the first Hopf
pd1_14=lc1_14('PD1')
pd1_bif_14=run(pd1_14,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
pd1_bif_14=run(pd1_bif_14('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})

# 5b. Branch switch at Second Period Doubling (the beginning) label
bsw1_14=lc1_14('PD1')
bsw1_bif_14=run(bsw1_14,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=-1, UZSTOP={})

# 6. Continue Second PD emerging from the Hopf in 2D (the beginning)
pd2_14=bsw1_bif_14('PD1')
pd2_bif_14=run(pd2_14,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
pd2_bif_bsw1_14=run(pd2_bif_14('EP'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})

# 7. Branch switch at Second Period Doubling (the beginning) label
SecondBsw_14=bsw1_bif_14('PD1')
SecondBsw_bif_14=run(SecondBsw_14,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=-1, UZSTOP={})


# 8. Continue Second PD emerging from the Hopf in 2D
pd3_14=SecondBsw_bif_14('PD1')
pd3_bif_14=run(pd3_14,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
pd3_bif_14=run(pd3_bif_14('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})


# 3. Plot results
# 3a. Some imports
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import BoundaryNorm
from matplotlib.transforms import Bbox
import numpy as np
import matplotlib.colors as mcolors
dict_color = mcolors.CSS4_COLORS


# 3b. Auxiliary functions
def pt_vals(f):
	return np.array([f[0][i]['PT'] for i in range(len(f[0]))])

def bifs(f,par):
	exceptions = ['No Label', 'RG', 'EP', 'UZ', 'LP']  # List of exceptions
	return [[f[0][i]['TY name'],f[0][i][par]] for i in range(len(f[0])) if f[0][i]['TY name'] not in exceptions]

col_map = colormaps['Purples']
boundary = BoundaryNorm([0,col_map.N],col_map.N)
Ntons = col_map.N

# Figure 1D: eps 12
fig1=plt.figure(figsize=(15,7))
#Customize axis
eps_title = 12
ax=plt.axes()
plt.title('$\epsilon =$ '+str(eps_title),fontsize=50,fontname='Times New Roman')
plt.xlabel('$I_{ext}^E$',fontsize=50,fontname='Times New Roman')
plt.ylabel('$v_E$ (mV)',fontsize=50,fontname='Times New Roman')
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
#Third PD left
plt.scatter(ThirdBsw_bif_12['Iext_e'],ThirdBsw_bif_12['y_max'], c=Ntons * (pt_vals(ThirdBsw_bif_12) < 0) + 2 * (pt_vals(ThirdBsw_bif_12) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(ThirdBsw_bif_12['Iext_e'],ThirdBsw_bif_12['y_min'], c=Ntons * (pt_vals(ThirdBsw_bif_12) < 0) + 2 * (pt_vals(ThirdBsw_bif_12) > 0), cmap='Purples', s=12,norm=boundary)

#Second PD right
plt.scatter(SecondBsw_bif_12['Iext_e'],SecondBsw_bif_12['y_max'], c=Ntons * (pt_vals(SecondBsw_bif_12) < 0) + 2 * (pt_vals(SecondBsw_bif_12) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(SecondBsw_bif_12['Iext_e'],SecondBsw_bif_12['y_min'], c=Ntons * (pt_vals(SecondBsw_bif_12) < 0) + 2 * (pt_vals(SecondBsw_bif_12) > 0), cmap='Purples', s=12,norm=boundary)
#Second PD left
#plt.scatter(SecondLeftBsw1_bif_12['Iext_e'],SecondLeftBsw1_bif_12['y_max'], c=Ntons * (pt_vals(SecondLeftBsw1_bif_12) < 0) + 2 * (pt_vals(SecondLeftBsw1_bif_12) > 0), cmap='Blues', s=12,norm=boundary)
#plt.scatter(SecondLeftBsw1_bif_12['Iext_e'],SecondLeftBsw1_bif_12['y_min'], c=Ntons * (pt_vals(SecondLeftBsw1_bif_12) < 0) + 2 * (pt_vals(SecondLeftBsw1_bif_12) > 0), cmap='Blues', s=12,norm=boundary)

#First PD left
plt.scatter(bsw_bif_12['Iext_e'],bsw_bif_12['y_max'], c=Ntons * (pt_vals(bsw_bif_12) < 0) + 2 * (pt_vals(bsw_bif_12) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(bsw_bif_12['Iext_e'],bsw_bif_12['y_min'], c=Ntons * (pt_vals(bsw_bif_12) < 0) + 2 * (pt_vals(bsw_bif_12) > 0), cmap='Purples', s=12,norm=boundary)

#First PD right
#plt.scatter(bsw2_bif_12['Iext_e'],bsw2_bif_12['y_max'], c=Ntons * (pt_vals(bsw1_bif_12) < 0) + 2 * (pt_vals(bsw1_bif_12) > 0), cmap='Greens', s=12,norm=boundary)
#plt.scatter(bsw2_bif_12['Iext_e'],bsw2_bif_12['y_min'], c=Ntons * (pt_vals(bsw1_bif_12) < 0) + 2 * (pt_vals(bsw1_bif_12) > 0), cmap='Greens', s=12,norm=boundary)
plt.scatter(lc1_12['Iext_e'],lc1_12['y_max'], c=Ntons * (pt_vals(lc1_12) < 0) + 2 * (pt_vals(lc1_12) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(lc1_12['Iext_e'],lc1_12['y_min'], c=Ntons * (pt_vals(lc1_12) < 0) + 2 * (pt_vals(lc1_12) > 0), cmap='Purples', s=12,norm=boundary)

#Equilibrium Point
plt.scatter(fp_12['Iext_e'],fp_12['v_e'], c=Ntons * (pt_vals(fp_12) < 0) + 2 * (pt_vals(fp_12) > 0), cmap='PuRd', s=12,norm=boundary)
#plt.axhline(0,color="black", ls="-")
plt.xlim([0,15])
#plt.ylim([7,14])
plt.savefig('Zoomed1D_PD12_Purple.png', dpi=600,bbox_inches=Bbox([[0,-1],fig1.get_size_inches()]))
plt.show()

# Figure 1D: eps 14
fig1=plt.figure(figsize=(15,7))
#Customize axis
ax=plt.axes()
plt.title('NextGeneration DynSynapses',fontsize=30,fontname='Times New Roman')
plt.xlabel('Parameter $I_{ext}^e$',fontsize=30,fontname='Times New Roman')
plt.ylabel('Parameter $\epsilon$',fontsize=30,fontname='Times New Roman')
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.scatter(lc1_14['Iext_e'],lc1_14['y_max'], c=Ntons * (pt_vals(lc1_14) < 0) + 2 * (pt_vals(lc1_14) > 0), cmap='Purples', s=12,norm=boundary)
plt.scatter(lc1_14['Iext_e'],lc1_14['y_min'], c=Ntons * (pt_vals(lc1_14) < 0) + 2 * (pt_vals(lc1_14) > 0), cmap='Purples', s=12,norm=boundary)
#First PD left
plt.scatter(bsw1_bif_14['Iext_e'],bsw1_bif_14['y_max'], c=Ntons * (pt_vals(bsw1_bif_14) < 0) + 2 * (pt_vals(bsw1_bif_14) > 0), cmap='Greens', s=12,norm=boundary)
plt.scatter(bsw1_bif_14['Iext_e'],bsw1_bif_14['y_min'], c=Ntons * (pt_vals(bsw1_bif_14) < 0) + 2 * (pt_vals(bsw1_bif_14) > 0), cmap='Greens', s=12,norm=boundary)
#First PD right
plt.scatter(bsw2_bif_14['Iext_e'],bsw2_bif_14['y_max'], c=Ntons * (pt_vals(bsw1_bif_14) < 0) + 2 * (pt_vals(bsw1_bif_14) > 0), cmap='Greens', s=12,norm=boundary)
plt.scatter(bsw2_bif_14['Iext_e'],bsw2_bif_14['y_min'], c=Ntons * (pt_vals(bsw1_bif_14) < 0) + 2 * (pt_vals(bsw1_bif_14) > 0), cmap='Greens', s=12,norm=boundary)
#Second PD left
plt.scatter(SecondBsw1_bif_14['Iext_e'],SecondBsw1_bif_14['y_max'], c=Ntons * (pt_vals(SecondBsw1_bif_14) < 0) + 2 * (pt_vals(SecondBsw1_bif_14) > 0), cmap='Blues', s=12,norm=boundary)
plt.scatter(SecondBsw1_bif_14['Iext_e'],SecondBsw1_bif_14['y_min'], c=Ntons * (pt_vals(SecondBsw1_bif_14) < 0) + 2 * (pt_vals(SecondBsw1_bif_14) > 0), cmap='Blues', s=12,norm=boundary)
#Equilibrium Point
plt.scatter(fp_14['Iext_e'],fp_14['v_e'], c=Ntons * (pt_vals(fp_14) < 0) + 2 * (pt_vals(fp_14) > 0), cmap='PuRd', s=12,norm=boundary)
#plt.axhline(0,color="black", ls="-")
plt.xlim([5,15])
#plt.ylim([7,14])
plt.savefig('1D_14.png', dpi=600,bbox_inches=Bbox([[0,-1],fig1.get_size_inches()]))
plt.show()


# Figure 2D
fig2=plt.figure(figsize=(15,7))
#Customize axis
ax=plt.axes()
ax.grid()
#plt.title('NextGeneration DynSynapses',fontsize=30,fontname='Times New Roman')
plt.xlabel('$I_{ext}^e$',fontsize=30,fontname='Times New Roman')
plt.ylabel('$\epsilon$',fontsize=30,fontname='Times New Roman')
#plt.xticks([2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],fontsize=30)
#plt.yticks([2,4,6,8,10,12,14,16,18,20,22,24,26],fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
#----------------------------  Plot Period Doubling region (eps=12)  ----------------------------
# eps = 12
#First Period doubling region --> Orange
plt.plot(pd1_bif_12['Iext_e'],pd1_bif_12['eps'], linewidth=2, color=dict_color["darkorange"])
#Second Period doubling region --> Orange
plt.plot(pd2_bif_12['Iext_e'],pd2_bif_12['eps'], linewidth=2, color=dict_color["red"])
#Third Period doubling region --> Orange
plt.plot(pd3_bif_12['Iext_e'],pd3_bif_12['eps'], linewidth=2, color=dict_color["darkred"])
#Third Period doubling region --> Orange
plt.plot(pd4_bif_12['Iext_e'],pd4_bif_12['eps'], linewidth=2, color=dict_color["darkred"])
#plt.plot(pd2_other_bif_12['Iext_e'],pd2_other_bif_12['eps'], linewidth=2, color=dict_color["blue"])

# eps = 14
#First Period doubling region --> Orange
#plt.plot(pd1_bif_14['Iext_e'],pd1_bif_14['eps'], linewidth=2, color=dict_color["darkorange"])
#Second Period doubling region --> Orange
plt.plot(pd2_bif_bsw1_14['Iext_e'],pd2_bif_bsw1_14['eps'], linewidth=2, color=dict_color["red"])
#Third Period doubling region --> Orange
plt.plot(pd3_bif_14['Iext_e'],pd3_bif_14['eps'], linewidth=2, color=dict_color["darkred"])


#---------------------------  Filled regions  -------------------------
# eps = 12
plt.fill(pd1_bif_12['Iext_e'],pd1_bif_12['eps'], color=dict_color["darkorange"],alpha=0.3)
plt.fill(pd2_bif_12['Iext_e'],pd2_bif_12['eps'], color=dict_color["red"],alpha=0.3)
plt.fill(pd3_bif_12['Iext_e'],pd3_bif_12['eps'], color=dict_color["darkred"],alpha=0.3)
plt.fill(pd4_bif_12['Iext_e'],pd4_bif_12['eps'], color=dict_color["darkred"],alpha=0.3)

# eps = 14
#plt.fill(pd1_bif_14['Iext_e'],pd1_bif_14['eps'], color=dict_color["darkorange"],alpha=0.3)
plt.fill(pd2_bif_bsw1_14['Iext_e'],pd2_bif_bsw1_14['eps'], color=dict_color["red"],alpha=0.3)
plt.fill(pd3_bif_14['Iext_e'],pd3_bif_14['eps'], color=dict_color["darkred"],alpha=0.3)



plt.xlim([8,12])
plt.ylim([9,15])


#plt.savefig('2D_PDComplet.png', dpi=600,bbox_inches=Bbox([[0,-1],fig1.get_size_inches()]))
plt.show()

