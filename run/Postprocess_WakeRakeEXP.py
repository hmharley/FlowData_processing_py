import sys
sys.path.append("D:/Ming_DOCUMENTS/Working/FlowData_processing_py/source/")
# sys.path.append("F:/HM/SURFDRIVE/Working/FlowData_processing_py/source/")
import math
import pandas as pd
import plt2pandas as p2p
import numpy as np
from timer_post import timer
from time import time
from MyPlot import MyPlot
from scipy.interpolate import griddata
from matplotlib.ticker import AutoMinorLocator

import matplotlib.pyplot as plt
import matplotlib.patches as patches
# plt.rc('text', usetex=True)
# plt.rc('grid', linestyle="dashed", color='black',linewidth=0.5)
FoldPath = "plotTest/TSR2p5/"
OutFile  = "plotTest/TSR2p5/"
VarList  = ['X','Y','Z','u']
U_inf = 5
Xplane = 2  # X/D

# %% load data
with timer("Load Meanflow data"):
   Mp = MyPlot()
   Mp.LoadPlt(FoldPath, VarList)


# Normalization
with timer('Normalization'):
    Mp.CoordNorm(300.0)
    Mp.VarNorm('u',U_inf)
    # Mp.VarNorm('u',U_inf)

# %% Save data
# with timer('Save flow data'):
#     # Mp.flowdata.to_csv('FlowData.csv')
#     store = pd.HDFStore('FlowData.h5')
#     store['df'] = Mp.flowdata  # save it
#     store['df']  # load it

# %%
# Plot
with timer('U contour at X/D = '+ str(Xplane)):
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-1,1], ylim = [-1,1], grid = True, size = [3.6,3])
#    v = np.linspace(0.4, 1.1, 8, endpoint=True) # normalized
   # v = np.linspace(2, 5.5, 8, endpoint=True)
   v = np.linspace(0.4,1.1,8,endpoint=True)

   cp,cl = Mp.PlaneContour(fig, ax, 'X', Xplane, 'u', cflevels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')   
   # cp = Mp.PlaneContour(fig, ax, 'X', 0, 'u', cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')   
   clb = fig.colorbar(cp)
   # ax.set_ylabel('$\it{Y/D}$')
   clb.ax.set_title('$\mathrm{U} / \mathrm{U}_{\infty} $')
   rect = patches.Rectangle((-0.5,-0.5),1,1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(rect)


# %% Integration
yg, zg, df = Mp.ExtractPlane('X',Xplane)
df['u'][np.isnan(df['u'])] = 0
VarInt = (U_inf - round(df['u'],2))*round(df['u'],2)
# Ug = griddata((df['Y'], df['Z']),df['u'],(yg,zg),method='linear')
varg = griddata((df['Y'], df['Z']), VarInt, (yg, zg), method='linear')
# varg[np.isnan(varg)] = 0
# grid size for uniform dicretezation
dy = yg[0,1]-yg[0,0]
dz = zg[1,0]-zg[0,0]

InteJ = np.zeros(varg.shape[0])
for j in range(varg.shape[0]):
    InteJ[j] = np.trapz(varg[j,:],dx=dy)
Integral = np.trapz(InteJ,dx=dz)
Ct = Integral/(0.5*1.2*5**2)
Ct
# meanVal = varg.mean()
# meanVal




# %%
