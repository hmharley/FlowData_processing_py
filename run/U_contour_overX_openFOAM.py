import sys
# sys.path.append("D:/Ming_DOCUMENTS/Working/FlowData_processing_py/source/")
sys.path.append("F:/HM/SURFDRIVE/Working/FlowData_processing_py/source/")
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
U_inf = 9.3
H = 1
# planevalue = [3,5,15,25]
AR = H/2/R # aspect ratio
planevalue = [1,2,3,4,5,6]
xlim = [-1,1]
ylim = [-H/2-0.5,H/2+0.5]
size = [2*len(planevalue)+1,2+0.5*(H-1)]
grid = False
# cflevel = np.linspace(0.4, 1.0, 8, endpoint=True)   # contourf level
clevel = np.linspace(0.1,1.1,11)
clblevel = np.linspace(0.5,1.1,7)
cmap = 'coolwarm'
contourline=True
ContourlineWidth = 0.6
ContourlineColors = 'k'
# cllevel = np.linspace(0.4, 1.0, 8, endpoint=True)
# -------------------------------------------------------------------------#
# %% Class MyPlot: Data Loading

FoldPath = "Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\TecioneVAWT\\AR1_k2em3_eplson_1em2\\postProcessing\\sampleDict\\5.004323150815\\"

# FoldPath = "Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\TecioneVAWT\\AR2\\postProcessing\\sampleDict\\3.604596474409887\\"

# FoldPath ="Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\TecioneVAWT\\AR3\\postProcessing\\sampleDict\\2.6041686237315\\"

OutPath  = "G:\\SURFdrive\\Working\\Paper_writing\\FirstPaper\\3rd_edition\\Figures\\"
casename = "AR"+str(AR)
plane = "UMean_planeXDstar" 
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
plt.rcParams['text.latex.preamble']=r'\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother'
plt.rc('grid', linestyle="dashed", color='black',linewidth=0.8,alpha=.5)
plt.set_cmap('coolwarm')
# plt.tight_layout()
# plt.xlabel('$\it{X/R}$')
fig, ax = plt.subplots(1,len(planevalue),constrained_layout=True,sharey=True)
fig.set_size_inches(size[0],size[1])
cp = [None]*(len(ax))
cl = [None]*(len(ax))
Ua = [None]*(len(ax))
for i in range(len(ax)):
    df = Mp[i].flowdata
    y = np.arange(np.min(df['Y']), np.max(df['Y']), 0.01)
    z = np.arange(np.min(df['Z']), np.max(df['Z']), 0.01)
    yg, zg = np.meshgrid(y, z)
    ux = griddata((df['Y'], df['Z']), df['u'], (yg, zg), method='cubic')
    ux = ux/U_inf
    cp[i] = ax[i].contourf( yg, zg, ux,clevel)
    cl[i] = ax[i].contour( yg, zg, ux,clevel,colors='k', linestyles = 'dashed', linewidths = 0.5)
    
    ax[i].set_xlim(xlim[0], xlim[1])
    ax[i].set_ylim(ylim[0], ylim[1])
    rect = patches.Rectangle((-R,-H/2),2*R,H,linewidth=1,edgecolor='k',facecolor='none',linestyle='dashed')
    ax[i].add_patch(rect)

    # averaged U in the frontal area
    ya = y[(y<=R)&(y>=-R)]
    za = z[(z<=(H/2))&(z>=(-H/2))]
    yga, zga = np.meshgrid(ya, za)
    uxa = griddata((df['Y'], df['Z']), df['u'], (yga, zga), method='cubic')
    uxa = uxa/U_inf
    Ua[i] = np.mean(uxa)

ax[0].set_ylabel('$\it{Z/R}$') 



fig.text(0.5, -0.05, '$\it{X/R}$', ha='center')
# fig.text(0, 0.9, r'C\rom{3}')
clb = fig.colorbar(cp[len(planevalue)-1])

if saveplot:
    plt.savefig(OutPath + casename + '_' + plane + '.png', dpi=300, bbox_inches = "tight")
    # plt.savefig(OutPath + casename + '_' + plane + '.png', dpi=300)

#  %%
