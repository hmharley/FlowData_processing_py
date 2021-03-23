import sys
# sys.path.append("D:/Ming_DOCUMENTS/Working/FlowData_processing_py/source/")
# sys.path.append("F:/HM/SURFDRIVE/Working/FlowData_processing_py/source/")
sys.path.append("G:/SURFdrive/Working/FlowData_processing_py/source/")
import math
import plt2pandas as p2p
import numpy as np
from timer_post import timer
from time import time
from MyPlot import MyPlot
from scipy.interpolate import griddata
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import variable_analysis as VA

# %%
# --------------------------- Control parameters --------------------------#
R = 0.15
H = 0.3
U_inf = 5
rho = 1.225
A = 2*R*H
planevalue = [3,8,15]
# planevalue = [1,3,5,10]
# planevalue = np.linspace(1,10,10,dtype=int)
# planevalue =np.r_[[-5,-3,-1],np.linspace(1,10,10,dtype=int)] 
# planevalue = [-5,-3,-1,1,2,3,4,5,6,20,25,30,35,40]
# limit for ploting
xlim = [-1.5,1.5]
ylim = [-1.5,1.5]
# limit for integral
ylim_int = [-3,3]
zlim_int = [-2.99,2.99]
size = [8,2]
grid = False
# cflevel = np.linspace(0.4, 1.0, 8, endpoint=True)   # contourf level
clevel = np.linspace(0.3,1.1,9)
# clevel =np.linspace(0.5,1.0,6)
clblevel = np.linspace(0.5,1.1,7)
cmap = 'coolwarm'
contourline=True
ContourlineWidth = 0.6
ContourlineColors = 'k'
# cllevel = np.linspace(0.4, 1.0, 8, endpoint=True)
# -------------------------------------------------------------------------#
# %% Class MyPlot: Data Loading

# FoldPath = "D:\\Ming_DOCUMENTS\\Working\\Numerical_Work\\Openfoam\\Cases\\ActuatorCylinder\\Ct067NewMesh\\postProcessing\\sampleDict\\39\\"

# -----------------------VAWT Plus 5-------------------------
# FoldPath = "F:\\HM\\SURFDRIVE\\Working\\Numerical_Work\\Openfoam\\Cases\\ActuatorLine\\TecioneVAWT\\Pitch_P5\\postProcessing\\sampleDict\\1.999829516\\"

# -----------------------VAWT Minus 5-------------------------
# FoldPath = "Z:\\scratch\\T_VAWT\\PitchMinus5\\postProcessing\\sampleDict\\1.9040081\\"

# -----------------------VAWT no pitch-------------------------
# FoldPath = "Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\TecioneVAWT\\postProcessing\\sampleDict\\3.404131987905258\\"

# FoldPath = "Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\SmallVAWT\\postProcessing\\sampleDict\\2.004\\"

# ----------------------AC NewMesh-----------------------------
FoldPath = "Z:\\scratch\\S_VAWT\\DoubleP10\\postProcessing\\sampleDict\\1.20044579502\\"
# FoldPath = "Z:\\OpenFOAM\\minghuang-4.1\\run\\ActuatorCylinder\\3DAC\\postProcessing\\sampleDict\\2.2\\"


# -----------------------Output path --------------------------
# OutPath  = "D:\\Ming_DOCUMENTS\\Working\\Paper_writing\\SimulationReport\\FIgures\\"
casename = "AC"
plane = "U_planeX"
saveplot = 0

# VarList  = ['x [mm]', 'y [mm]', 'z [mm]', 'Vx [m/s]', 'Vy [m/s]', 'Vz [m/s]', 'isValid']
VarList  = ['X','Y','Z','u','v','w']
Mp = [None]*(len(planevalue))
with timer("Load Meanflow data"):
   for i in range(len(planevalue)):
      Mp[i] = MyPlot()
      Mp[i].LoadRaw(FoldPath+plane+str(planevalue[i])+'.raw',VarList)
      Mp[i].CoordNorm(2*R)

# %% Ct U_plane
with timer("Integration"):
    uw = [None]*(len(planevalue))
    Drag = [None]*(len(planevalue))
    KEnergy = [None]*(len(planevalue))
    for i in range(len(planevalue)):
        df = Mp[i].flowdata
        # round(df['u'],2)
        # VarInt = (U_inf - df['u'])*df['u']
        n = 100
        y = np.linspace(ylim_int[0], ylim_int[1], n)
        z = np.linspace(zlim_int[0], zlim_int[1], n)
        # dfinZone = df.loc[(df['Y']>=np.min(y)) & (df['Y']<=np.max(y)) &
                # (df['Z']>=np.min(z)) & (df['Z']<=np.max(z))]
        DragInt = (U_inf - round(df['u'],2))*round(df['u'],2)
        KEnergyInt = 0.5*rho*(round(df['u'],2)**3)                      

        uw[i] = VA.integral_db(df['Y'], df['Z'], df['u'],
                                     range1=y, range2=z, opt=2)
        Drag[i] = VA.integral_db(df['Y'], df['Z'], DragInt,
                                     range1=y, range2=z, opt=2)
        KEnergy[i] = VA.integral_db(df['Y'], df['Z'], KEnergyInt,
                                     range1=y, range2=z, opt=2)
# Drag = VA.integral_db(df['Y'], df['Z'], VarInt,
                                    #  range1=df['Y'], range2=df['Z'], opt=2)
# Drag = VA.integral_db(df['Y'], df['Z'], VarInt,
                                    #  range1=y, range2=z, opt=3)                                     
Ct = np.true_divide(Drag,0.5*rho*U_inf**2*A)  
Ct                 
# %%
plt.plot(planevalue,Ct,label='Ct')
legend = plt.legend(shadow=True, fontsize='x-large')
plt.plot(planevalue,uw,label='$U_w$')
legend = plt.legend(shadow=True, fontsize='x-large')
plt.plot(planevalue,KEnergy,label='Kinetic Energy')
legend = plt.legend(shadow=True, fontsize='x-large')

legend = plt.legend(shadow=True, fontsize='x-large')
# yg, zg = np.meshgrid(y, z)
# ux = griddata((df['Y'], df['Z']), df['u'], (yg, zg), method='cubic')
# plt.contourf( yg, zg, ux)

# %%

plane = "p_planeX"

# VarList  = ['x [mm]', 'y [mm]', 'z [mm]', 'Vx [m/s]', 'Vy [m/s]', 'Vz [m/s]', 'isValid']
VarList  = ['X','Y','Z','p']
Mpp = [None]*(len(planevalue))
with timer("Load Meanflow data"):
   for i in range(len(planevalue)):
      Mpp[i] = MyPlot()
      Mpp[i].LoadRaw(FoldPath+plane+str(planevalue[i])+'.raw',VarList)
      Mpp[i].CoordNorm(2*R)
with timer("Integration"):
    p = [None]*(len(planevalue))
    for i in range(len(planevalue)):
        df = Mpp[i].flowdata
        # round(df['u'],2)
        # VarInt = (U_inf - df['u'])*df['u']
        n = 100
        y = np.linspace(ylim_int[0], ylim_int[1], n)
        z = np.linspace(zlim_int[0], zlim_int[1], n)
        # dfinZone = df.loc[(df['Y']>=np.min(y)) & (df['Y']<=np.max(y)) &
                # (df['Z']>=np.min(z)) & (df['Z']<=np.max(z))]                 
        p[i] = VA.integral_db(df['Y'], df['Z'], df['p'],
                                     range1=y, range2=z, opt=2)    
plt.plot(planevalue,p,label='p')


# %% Plot
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r'\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother'
plt.rc('grid', linestyle="dashed", color='black',linewidth=0.8,alpha=.5)
plt.set_cmap('coolwarm')
# plt.tight_layout()
# plt.xlabel('$\it{X/R}$')
fig, ax = plt.subplots(1,len(planevalue),constrained_layout=True,sharey=True)
fig.set_size_inches(size[0],size[1])
cp = [None]*(len(ax))
cl = [None]*(len(ax))
for i in range(len(ax)):
    df = Mp[i].flowdata
    y = np.arange(np.min(df['Y']), np.max(df['Y']), 0.01)
    z = np.arange(np.min(df['Z']), np.max(df['Z']), 0.01)
    # y = np.linspace(-0.5, 0.5, n)
    # z = np.linspace(-0.5, 0.5, n)    
    yg, zg = np.meshgrid(y, z)
    ux = griddata((df['Y'], df['Z']), df['u'], (yg, zg), method='cubic')
    ux = ux/U_inf
    cp[i] = ax[i].contourf( yg, zg, ux,clevel)
    cl[i] = ax[i].contour( yg, zg, ux,clevel,colors='k', linestyles = 'dashed', linewidths = 0.5)
    
    ax[i].set_xlim(xlim[0], xlim[1])
    ax[i].set_ylim(ylim[0], ylim[1])
    rect = patches.Rectangle((-0.5,-0.5),1,1,linewidth=1,edgecolor='k',facecolor='none',linestyle='dashed')
    ax[i].add_patch(rect)
ax[0].set_ylabel('$\it{Z/R}$') 

fig.text(0.5, -0.05, '$\it{X/R}$', ha='center')
fig.text(0, 0.9, r'C\rom{3}')
clb = fig.colorbar(cp[len(planevalue)-1])

if saveplot:
    plt.savefig(OutPath + casename + '_' + plane + '.png', dpi=300, bbox_inches = "tight")
    plt.savefig(OutPath + casename + '_' + plane + '.png', dpi=300)

#  %%
