import sys
sys.path.append("D:/Ming_DOCUMENTS/Working/FlowData_processing_py/source/")
# sys.path.append("F:/HM/SURFDRIVE/Working/FlowData_processing_py/source/")
# sys.path.append("G:/SURFdrive/Working/FlowData_processing_py/source/")
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

# %%
# --------------------------- Control parameters --------------------------#
R = 0.5
U_inf = 1
# planevalue = [3,5,15,25]
planevalue = [1,3,5]
xlim = [-1,1]
ylim = [-1,1]
size = [7,2]
grid = False
# cflevel = np.linspace(0.4, 1.0, 8, endpoint=True)   # contourf level
clevel = np.linspace(0.3,1.1,9)
clblevel = np.linspace(0.5,1.1,7)
cmap = 'coolwarm'
contourline=True
ContourlineWidth = 0.6
ContourlineColors = 'k'
# cllevel = np.linspace(0.4, 1.0, 8, endpoint=True)
# -------------------------------------------------------------------------#
# %% Class MyPlot: Data Loading
# FoldPath = "DATA\\"
FoldPath = "D:\\Ming_DOCUMENTS\\Working\\Numerical_Work\\Openfoam\\Cases\\ActuatorCylinder\\Ct067coarse\\postProcessing\\sampleDict\\34\\"
OutPath  = "FIGURE\\"
casename = "CI"
plane = "U_planeX"
saveplot = 0

# VarList  = ['x [mm]', 'y [mm]', 'z [mm]', 'Vx [m/s]', 'Vy [m/s]', 'Vz [m/s]', 'isValid']
VarList  = ['X','Y','Z','u','v','w']
Mp = [None]*(len(planevalue))
with timer("Load Meanflow data"):
   for i in range(len(planevalue)):
      Mp[i] = MyPlot()
      Mp[i].LoadRaw(FoldPath+plane+str(planevalue[i])+'.raw',VarList)

# %% Plot
plt.rc('text', usetex=True)
plt.rc('grid', linestyle="dashed", color='black',linewidth=0.8,alpha=.5)
plt.set_cmap('coolwarm')
# plt.xlabel('$\it{X/R}$')
fig, ax = plt.subplots(1,len(planevalue),constrained_layout=True,sharey=True)
fig.set_size_inches(size[0],size[1])
cp = [None]*(len(ax))
cl = [None]*(len(ax))
for i in range(len(ax)):
    df = Mp[i].flowdata
    y = np.arange(np.min(df['Y']), np.max(df['Y']), 0.01)
    z = np.arange(np.min(df['Z']), np.max(df['Z']), 0.01)
    yg, zg = np.meshgrid(y, z)
    ux = griddata((df['Y'], df['Z']), df['u'], (yg, zg), method='cubic')
    cp[i] = ax[i].contourf( yg, zg, ux,clevel)
    cl[i] = ax[i].contour( yg, zg, ux,clevel,colors='k', linestyles = 'dashed', linewidths = 0.5)
    
    ax[i].set_xlim(xlim[0], xlim[1])
    ax[i].set_ylim(ylim[0], ylim[1])
    rect = patches.Rectangle((-0.5,-0.5),1,1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
    ax[i].add_patch(rect)
ax[0].set_ylabel('$\it{Z/R}$') 

fig.text(0.5, -0.05, '$\it{X/R}$', ha='center')
clb = fig.colorbar(cp[len(planevalue)-1])

if saveplot:
    plt.savefig(OutPath + casename + '_' + plane + '.png', dpi=300)








# %%
