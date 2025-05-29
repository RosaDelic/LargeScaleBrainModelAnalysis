
#=====================  PART HOPF + Bistability eps = 15  ==========================
eps = 15
# 1. Use Euler integrator up to steady state to initialize the diagram
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':0.0, 'eps' : eps})

#2. Continue along p an uncover the bifurcations (2 Hopf)
ic=init(201)
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#save(fp,'fp')

# 3. Continue the limit cycle solutions (ISW=1,IPS=2) emerging from the first Hopf.
hb=fp('HB1')
lc1=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

# 3. Continue LP emerging from the first Hopf (bistability region ~ Green) Both LPs give the same closed curve
lp1=lc1('LP1')
sd_lc=run(lp1,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
sd_lc_bif=run(sd_lc('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#This is not necessary bc previous step ends in EP not in MX
#lp2=lc1('LP2')
#sd_lc2=run(lp2,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#sd_lc_bif_2=run(sd_lc2('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#sd_lc_bif_2=run(lp2,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={},DS="-")

# 4. Continue Hopf region in 2D (two branches needed) --> Several GH and BT appear
lc1a_2d=run(hb,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
lc1b_2d=run(hb,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={},DS="-")
lc1_2d=merge(lc1a_2d+lc1b_2d)
lc1_2d_rel=relabel(lc1_2d)


#=====================  PART HOPF + Bistability eps = 18  ==========================
eps = 18
# 1. Use Euler integrator up to steady state to initialize the diagram
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':0.0, 'eps' : eps})

#2. Continue along p an uncover the bifurcations (2 Hopf)
ic=init(201)
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})
#save(fp,'fp')

# 3. Continue the limit cycle solutions (ISW=1,IPS=2) emerging from the first Hopf.
hb_18=fp('HB1')
lc1_18=run(hb_18,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

# 3. Continue LP emerging from the first Hopf (bistability region ~ Green) Both LPs give the same closed curve
lp1_18=lc1_18('LP1')
sd_lc_18=run(lp1_18,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
sd_lc_bif_18=run(sd_lc_18('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#This is not necessary bc previous step ends in EP not in MX
#lp2=lc1('LP2')
#sd_lc2=run(lp2,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#sd_lc_bif_2=run(sd_lc2('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#sd_lc_bif_2=run(lp2,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={},DS="-")


#=====================  PART Period Doubling eps = 12  ==========================
eps = 12
# 1. Use Euler integrator up to steady state to initialize the diagram
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':0.0, 'eps' : eps})

#2. Continue along p an uncover the bifurcations (2 Hopf)
ic=init(201)
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})

# 3. Continue the limit cycle solutions (ISW=1,IPS=2) emerging from the first Hopf.
hb_12=fp('HB1')
lc1_12=run(hb_12,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

# 3.5. Branch switch at Period Doubling
#pd1_12_second=lc1_12('PD2')
#pd1_bif_12_second=run(pd1_12_second,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=-1, UZSTOP={})

# 4. Continue PDs emerging from the first Hopf (Period Doubling Region ~ Orange)
pd1_12=lc1_12('PD1')
pd1_bif_12=run(pd1_12,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
pd1_bif_12=run(pd1_bif_12('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#Both PD labels give the same curve
#pd2_12=lc1_12('PD2')
#pd2_bif_12=run(pd2_12,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#pd2_bif_12=run(pd2_bif_12('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})


#===================  PART Bistability + Resonances eps = 25  ==========================
eps = 25
# 1. Use Euler integrator up to steady state to initialize the diagram
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':0.0, 'eps' : eps})

#2. Continue along p an uncover the bifurcations (2 Hopf)
ic=init(201)
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})

# 3. Continue the limit cycle solutions (ISW=1,IPS=2) emerging from the first Hopf
hb_25=fp('HB1')
lc1_25=run(hb_25,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

# 4. Continue LPs observed in the continuation of the first Hopf  (Resonance Region ~ Red) --> There appears 1 R1 labels
lp1_25=lc1_25('LP1')
sd_lc_25=run(lp1_25,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
sd_lc_bif_25=run(sd_lc_25('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
lp1_25=sd_lc_bif_25('LP2')
sd_lc_bif_25=run(lp1_25,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={},DS="-")

#=====================  PART SNIC + Period Doubling + Resonance eps = 29  ==========================
eps = 29
# 1. Use Euler integrator up to steady state to initialize the diagram
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':0.0, 'eps' : eps})

#2. Continue along p an uncover the bifurcations (2 Hopf + 2 LP)
ic=init(201)
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})

# 3. Continue the Saddle-Node solutions (ISW=2,IPS=1) emerging from the LP --> Lables CP and BT appear
sn_29 = fp('LP1')
sn1_bif_29=run(sn_29,IPS=1,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
sn2_bif_29=run(sn_29,IPS=1,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={},DS="-")
sn_bif_29=merge(sn1_bif_29+sn2_bif_29)
sn_bif_29=relabel(sn_bif_29)
#sn_25 = fp('LP2')
#sn2_bif_25=run(sn_25,IPS=1,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})

# 4. Continue the limit cycle solutions (ISW=1,IPS=2) emerging from the first Hopf (Period Doublings)
#hb_25=fp('HB1')
#lc1_25=run(hb_25,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})

# 3. Continue the Period Doubling solutions (ISW=2,IPS=1) emerging from the first Hopf --> Lables LP appear
#pd1_25=lc1_25('PD1')
#pd1_bif_25=run(pd1_25,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#pd1_bif_25=run(pd1_bif_25('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})

#pd2_25=lc1_25('PD2')
#pd2_bif_25=run(pd2_25,IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})
#pd2_bif_25=run(pd2_bif_25('EP1'),IPS=2,ISP=2,ICP=[1,2,11,3,4],NMX=20000,ISW=2, UZSTOP={})

#=====================  PART Homoclinic eps = 35  ==========================
eps = 35
# 1. Use Euler integrator up to steady state to initialize the diagram
init=run('ngsyn',IPS=-2,NMX=100000,PAR={'Iext_e':0.0, 'eps' : eps})

#2. Continue along p an uncover the bifurcations (2 Hopf + 2 LP)
ic=init(201)
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'Iext_e' : 0, 'Iext_e' : 100})

# 4. Continue the limit cycle solutions (ISW=1,IPS=2) emerging from the first Hopf stop when period is huge to find homoclinics
hb_35=fp('HB1')
#lc1_30=run(hb_30,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={})
lc1_35=run(hb_35,IPS=2,ICP=[1,11,2,3,4],NMX=300000,ISW=1,UZSTOP = {'PERIOD' : 300}) #1500 ---> hom_bad
hom1_35 = run(lc1_35('UZ'),ISW=1,ICP=[1,2,11],ISP=0)
hom2_35 = run(lc1_35('UZ'),ISW=1,ICP=[1,2,11],ISP=0,DS="-")
hom_35 = merge(hom1_35+hom2_35)
hom_35=relabel(hom_35)


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

# Figure
fig1=plt.figure(figsize=(15,7))
#Customize axis
ax=plt.axes()
plt.title('NextGeneration DynSynapses',fontsize=30,fontname='Times New Roman')
plt.xlabel('Parameter $I_{ext}^e$',fontsize=30,fontname='Times New Roman')
plt.ylabel('Parameter $\epsilon$',fontsize=30,fontname='Times New Roman')
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
#----------------------------  Plot Period Doubling region (eps=12)  ----------------------------
#Period doubling region --> Orange
plt.scatter(pd1_bif_12['Iext_e'],pd1_bif_12['eps'], c=50 * (pt_vals(pd1_bif_12)<0), cmap='Oranges', s=12,norm=boundary)
#plt.scatter(pd2_bif_12['Iext_e'],pd2_bif_12['eps'], c=50 * (pt_vals(pd2_bif_12)<0), cmap='Oranges', s=12,norm=boundary)

#----------------------------  Plot Hopf and bistability region (eps=15)  ----------------------------
#Hopf --> Purple
plt.scatter(lc1_2d['Iext_e'],lc1_2d['eps'], c=Ntons * (pt_vals(lc1_2d)), cmap='Purples', s=12,norm=boundary)
#Bistability region --> Green
plt.scatter(sd_lc_bif['Iext_e'],sd_lc_bif['eps'], c=50 * (pt_vals(sd_lc_bif)<0), cmap='Greens', s=12,norm=boundary)

#----------------------------  Plot Hopf and bistability region (eps=18)  ----------------------------
#Bistability region --> Green
plt.scatter(sd_lc_bif_18['Iext_e'],sd_lc_bif_18['eps'], c=50 * (pt_vals(sd_lc_bif_18)<0), cmap='Greens', s=12,norm=boundary)

#----------------------------  Plot Period doubling and Resonance region (eps=25)  ----------------------------
#Period doubling region --> Orange
plt.scatter(sd_lc_bif_25['Iext_e'],sd_lc_bif_25['eps'], c=50 * (pt_vals(sd_lc_bif_25)<0), cmap='Oranges', s=12,norm=boundary)

#----------------------------  Plot SNIC + Period Doubling + Resonances region (eps=25)  ----------------------------
#SN/SNIC --> Blue
plt.scatter(sn_bif_29['Iext_e'],sn_bif_29['eps'], c=50 * (pt_vals(sn_bif_29)<0), cmap='Blues', s=12,norm=boundary)
#plt.scatter(sn1_bif_25['Iext_e'],sn1_bif_25['eps'], c=50 * (pt_vals(sn1_bif_25)<0), cmap='Blues', s=12,norm=boundary)
#plt.scatter(sn2_bif_25['Iext_e'],sn2_bif_25['eps'], c=50 * (pt_vals(sn2_bif_25)<0), cmap='Blues', s=12,norm=boundary)

#-------------------------------------  Plot Homoclinic region (eps=35)  ---------------------------------------
#Homoclinic --> Brown
plt.scatter(hom_35['Iext_e'],hom_35['eps'], c=Ntons * (pt_vals(hom_35)), cmap='Reds', s=12,norm=boundary)
#plt.scatter(hom1_35['Iext_e'],hom1_35['eps'], c=Ntons * (pt_vals(hom1_35)), cmap='Reds', s=12,norm=boundary)
#plt.scatter(hom2_35['Iext_e'],hom2_35['eps'], c=Ntons * (pt_vals(hom2_35)), cmap='Reds', s=12,norm=boundary)


#bfp = bifs(fp,'Iext_e')
#for b in bfp:
#	plt.axvline(b[1],color="black", ls="--",alpha=0.7)
#bfp = bifs(lc1_2d,'Iext_e')
#for b in bfp:
#    plt.axvline(b[1],color="red", ls="--",alpha=0.7)
plt.axhline(0,color="black", ls="-")
plt.xlim([-5,20])
plt.ylim([0,130])
plt.savefig('2D_1.png', dpi=600,bbox_inches=Bbox([[0,-1],fig1.get_size_inches()]))
plt.show()

# Create a mask for the condition y < alpha
mask = lc1_2d['eps'] < 32.4245

#Load data Floquet Bifurcation diagram
#First right curve
#dataCompleteRight = np.load('CompleteRightMost.npz')
#vector_epsCompleteRight = dataCompleteRight['vector_eps']
#vector_Iext_eCompleteRight = dataCompleteRight['vector_Iext_e']
#First left curve
#dataCompleteLeft = np.load('CompleteLeftMost.npz')
#vector_epsCompleteLeft = dataCompleteLeft['vector_eps']
#vector_Iext_eCompleteLeft = dataCompleteLeft['vector_Iext_e']
#Second right curve
#dataCompleteSecondRight = np.load('CompleteSecondRightMost.npz')
#vector_epsCompleteSecondRight = dataCompleteSecondRight['vector_eps']
#vector_Iext_eCompleteSecondRight = dataCompleteSecondRight['vector_Iext_e']
#Second left curve
#dataCompleteSecondLeft = np.load('CompleteRoundSecondLeftMost.npz')
#vector_epsCompleteSecondLeft = dataCompleteSecondLeft['vector_eps']
#vector_Iext_eCompleteSecondLeft = dataCompleteSecondLeft['vector_Iext_e']

# Figure
fig2=plt.figure(figsize=(15,7))
fig2=plt.figure(figsize=(10,8))
#Customize axis
ax=plt.axes()
#ax.grid()
#plt.title('NextGeneration DynSynapses',fontsize=30,fontname='Times New Roman')
plt.xlabel('$I_{ext}^e$',fontsize=30,fontname='Times New Roman')
plt.ylabel('$\epsilon$',fontsize=30,fontname='Times New Roman')
#plt.xticks([2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],fontsize=30)
#plt.yticks([2,4,6,8,10,12,14,16,18,20,22,24,26],fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#----------------------------  Plot Period Doubling region (eps=12)  ----------------------------
#Period doubling region --> Orange
plt.plot(pd1_bif_12['Iext_e'],pd1_bif_12['eps'], linewidth=2, color=dict_color["darkorange"])
#plt.scatter(pd2_bif_12['Iext_e'],pd2_bif_12['eps'], c=50 * (pt_vals(pd2_bif_12)<0), cmap='Oranges', s=12,norm=boundary)

#----------------------------  Plot Hopf and bistability region (eps=15)  ----------------------------
#Hopf --> Purple
plt.plot(lc1_2d['Iext_e'],lc1_2d['eps'], linewidth=2, color=dict_color["indigo"])

#Bistability region --> Green
plt.plot(sd_lc_bif['Iext_e'],sd_lc_bif['eps'], linewidth=2, color=dict_color["seagreen"])

#----------------------------  Plot Hopf and bistability region (eps=18)  ----------------------------
#Bistability region --> Green
plt.plot(sd_lc_bif_18['Iext_e'],sd_lc_bif_18['eps'], linewidth=2, color=dict_color["seagreen"])

#----------------------------  Plot Period doubling and Resonance region (eps=25)  ----------------------------
#Period doubling region --> Orange
plt.plot(sd_lc_bif_25['Iext_e'],sd_lc_bif_25['eps'], linewidth=2, color=dict_color["seagreen"])

#----------------------------  Plot SNIC + Period Doubling + Resonances region (eps=25)  ----------------------------
#SN/SNIC --> Blue
plt.plot(sn_bif_29['Iext_e'],sn_bif_29['eps'], linewidth=2, color=dict_color["dodgerblue"])
#plt.scatter(sn1_bif_25['Iext_e'],sn1_bif_25['eps'], c=50 * (pt_vals(sn1_bif_25)<0), cmap='Blues', s=12,norm=boundary)
#plt.scatter(sn2_bif_25['Iext_e'],sn2_bif_25['eps'], c=50 * (pt_vals(sn2_bif_25)<0), cmap='Blues', s=12,norm=boundary)

#-------------------------------------  Plot Homoclinic region (eps=35)  ---------------------------------------
#Homoclinic --> Brown
#plt.plot(hom_35['Iext_e'],hom_35['eps'], linewidth=2, color=dict_color["firebrick"])
#plt.scatter(hom1_35['Iext_e'],hom1_35['eps'], c=Ntons * (pt_vals(hom1_35)), cmap='Reds', s=12,norm=boundary)
#plt.scatter(hom2_35['Iext_e'],hom2_35['eps'], c=Ntons * (pt_vals(hom2_35)), cmap='Reds', s=12,norm=boundary)


#---------------------------  Filled regions  -------------------------
#Hopf --> Purple
plt.fill(lc1_2d['Iext_e'][mask],lc1_2d['eps'][mask], color=dict_color["indigo"],alpha=0.3)
#Period doubling region --> Orange
plt.fill(pd1_bif_12['Iext_e'],pd1_bif_12['eps'], color=dict_color["white"],alpha=0.3)
plt.fill(pd1_bif_12['Iext_e'],pd1_bif_12['eps'], color=dict_color["darkorange"],alpha=0.3)
#Bistability region --> Green
#eps 15
plt.fill(sd_lc_bif['Iext_e'],sd_lc_bif['eps'], color=dict_color["white"],alpha=0.3)
plt.fill(sd_lc_bif['Iext_e'],sd_lc_bif['eps'], color=dict_color["seagreen"],alpha=0.3)
#eps 18
plt.fill(sd_lc_bif_18['Iext_e'],sd_lc_bif_18['eps'], color=dict_color["white"],alpha=0.3)
plt.fill(sd_lc_bif_18['Iext_e'],sd_lc_bif_18['eps'], color=dict_color["seagreen"],alpha=0.3)

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

#----------------------------  Plot Hopf and bistability region (eps=15)  ----------------------------

#Bistability region --> Green
plt.plot(sd_lc_bif['Iext_e'],sd_lc_bif['eps'], linewidth=2, color=dict_color["seagreen"])

#Bistability region --> Green
#eps 15
plt.fill(sd_lc_bif['Iext_e'],sd_lc_bif['eps'], color=dict_color["white"],alpha=0.3)
plt.fill(sd_lc_bif['Iext_e'],sd_lc_bif['eps'], color=dict_color["seagreen"],alpha=0.3)


#bfp = bifs(fp,'Iext_e')
#for b in bfp:
#    plt.axvline(b[1],color="black", ls="--",alpha=0.7)
#bfp = bifs(lc1_2d,'Iext_e')
#for b in bfp:
#    plt.axvline(b[1],color="red", ls="--",alpha=0.7)

#----------------------- Examples -------------------------
#epsilonsDestab = np.array([5,5,7.5,7.5,11,11,16,16,17.5,20,20])
#Iext_eDestab = np.array([11.5,14,13,15,6,8,7.5,9.5,9.5,6.5,9])
#epsilonsStab = np.array([5,7.5,11,16,17.5,17.5,20])
#Iext_eStab = np.array([10.5,12,5,6.5,8.5,10.5,7.5])

#plt.scatter(Iext_eStab,epsilonsStab,s=20,color='green')
#plt.scatter(Iext_eDestab,epsilonsDestab,s=20,color='red')

#-----------------------  Floquet bifurcation line  -------------------------
#plt.plot(vector_Iext_eCompleteLeft,vector_epsCompleteLeft, linewidth=2, color=dict_color["slategray"])
#plt.plot(vector_Iext_eCompleteRight,vector_epsCompleteRight, linewidth=2, color=dict_color["slategray"])
#plt.plot(vector_Iext_eCompleteSecondRight,vector_epsCompleteSecondRight, linewidth=2, color=dict_color["slategray"])
#plt.plot(vector_Iext_eCompleteSecondLeft,vector_epsCompleteSecondLeft, linewidth=2, color=dict_color["slategray"])

plt.axhline(0,color="black", ls="-")


# A. Total Diagram
plt.xlim([-10,20])
plt.ylim([-1,200])

# B. Oscillatory Region
plt.xlim([0,18])
plt.ylim([0,45])

# C. Period Doubling
plt.xlim([7.5,12.5])
plt.ylim([8,16])

# D. Bistability
plt.xlim([7.5,12.5])
plt.ylim([8,25])

# E. UpperPart
plt.xlim([0,10])
plt.ylim([17,36])

# F. MostUpperPart
plt.xlim([-4,5])
plt.ylim([25,140])

# Originals
plt.xlim([-5,15]) #Originals
plt.ylim([0,140])

plt.xlim([8,12]) #Originals
plt.ylim([9,15])

#plt.xlim([-10,20])
#plt.ylim([-10,200])
#plt.xlim([7.75,12.25])
#plt.ylim([7,16])
#plt.xlim([0.5,10])
#plt.ylim([20,36])
#plt.xlim([-4,6])
#plt.ylim([20,122])
plt.savefig('Definitive/2D_Diagram/Pdf_PeriodDoubling.png', dpi=600,bbox_inches=Bbox([[0,-1],fig2.get_size_inches()]))
plt.show()

