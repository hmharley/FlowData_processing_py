import sys
# sys.path.append("D:/Ming_DOCUMENTS/Working/FlowData_processing_py/source/")
sys.path.append("F:/HM/SURFDRIVE/Working/FlowData_processing_py/source/")
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
R = 150
U_inf = 5.0
planevalue = [1,2,3,4]
xlim = [-2,2]
ylim = [-2,2]
size = [9,2]
grid = False
cflevel = np.linspace(0.4, 1.2, 8, endpoint=True)   # contourf level
cmap = 'coolwarm'
contourline=True
ContourlineWidth = 0.6
ContourlineColors = 'k'
# cllevel = np.linspace(0.4, 1.0, 8, endpoint=True)
# -------------------------------------------------------------------------#
# %% Class MyPlot: Data Loading
FoldPath = "plotTest/"
OutFile  = "plotTest/"
# VarList  = ['x [mm]', 'y [mm]', 'z [mm]', 'Vx [m/s]', 'Vy [m/s]', 'Vz [m/s]', 'isValid']
VarList  = ['X','Y','Z','u','v','w']

with timer("Load Meanflow data"):
   Mp = MyPlot()
   Mp.LoadPlt(FoldPath, VarList)

with timer('Normalization'):
   Mp.CoordNorm(R)
   Mp.VarNorm('u',U_inf)
   Mp.VarNorm('v',U_inf)

# %% Plot
plt.rc('text', usetex=True)
plt.rc('grid', linestyle="dashed", color='black',linewidth=0.8,alpha=.5)
fig, ax = plt.subplots(1,4,constrained_layout=True)

cp = [None]*(len(ax))
cl = [None]*(len(ax))
for i in range(len(ax)):
    Mp.PlotStyle(fig, ax[i], xlim = xlim, ylim = ylim, grid = grid, size = size)

    cp[i],cl[i] = Mp.PlaneContour(fig, ax[i], 'X', planevalue[i], 'u', cflevels=cflevel, cllevels=None, cmap = cmap, axislabel = None,contourline=contourline, linewidths=ContourlineWidth, colors=ContourlineColors)

    rect = patches.Rectangle((-1,-1),2,2,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
    ax[i].add_patch(rect)
ax[0].set_ylabel('$\it{Z/R}$')    
clb = fig.colorbar(cp[3])









# %%
