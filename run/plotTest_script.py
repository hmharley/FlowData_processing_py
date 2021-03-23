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
plt.rc('text', usetex=True)
# plt.rc('grid', linestyle="dashed", color='black',linewidth=0.5)


# %% Class MyPlot: Data Loading
FoldPath = "plotTest/"
OutFile  = "plotTest/"
# VarList  = ['x [mm]', 'y [mm]', 'z [mm]', 'Vx [m/s]', 'Vy [m/s]', 'Vz [m/s]', 'isValid']
VarList  = ['X','Y','Z','u','v','w','u_std','Om_y','Om_z']
U_inf = 5.0
with timer("Load Meanflow data"):
   Mp = MyPlot()
   Mp.LoadPlt(FoldPath, VarList)

# Processing 
with timer('Normalization'):
   # Mp.CoordNorm(150.0)
   Mp.VarNorm('X',150)
   Mp.VarNorm('Y',150)
   Mp.VarNorm('Z',300)
   Mp.VarNorm('u',U_inf)
   Mp.VarNorm('v',U_inf)
   Mp.VarNorm('w',U_inf)
   Mp.VarNorm('u_std',U_inf)
   Mp.VarNorm('Om_y',U_inf/0.03)
   Mp.VarNorm('Om_z',U_inf/0.03)


plt.rc('text', usetex=True)
plt.rc('grid', linestyle="dashed", color='black',linewidth=0.8,alpha=.5)

# %% Plotting
# plot1
with timer('U contour at Z/D = 0'):
   
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-1.5,6], ylim = [-2,1], grid = True, size = [9,3])
   # v = np.linspace(0.0, 1.2, 12, endpoint=True)
   v = np.arange(0, 1.28, 0.08 )
   # v = np.arange(0,1.4,0.2)
   cp,cl = Mp.PlaneContour(fig, ax, 'Z', 0, 'u', cflevels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')
   t = np.arange(0,1.4,0.2)   
   clb = fig.colorbar(cp,ticks=t)
   clb.ax.set_title('$\mathrm{U} / \mathrm{U}_{\infty} $')
   cir = patches.Circle((0,0),1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(cir)   
   # fig.savefig('U_Zr0', dpi=330)

# %% plot2
with timer('V contour at Z/D = 0'):
      
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-1.5,6], ylim = [-2,1], grid = True, size = [9,3])
   v = np.arange(-0.3, 0.4, 0.1, )
   # cp = Mp.PlaneContour(fig, ax, 'Z', 0, 'v', levels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')   
   cp,cl = Mp.PlaneContour(fig, ax, 'Z', 0, 'v', cflevels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k') 
   t = np.arange(-0.3,0.4,0.1)    
   clb = fig.colorbar(cp,ticks=t)
   clb.ax.set_title('$\mathrm{V} / \mathrm{V}_{\infty} $')
   cir = patches.Circle((0,0),1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(cir)   
   fig.savefig('V_Zr0', dpi=330)   

# %% plot3
with timer('U\' contour at Z/D = 0'):
      
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-1.5,6], ylim = [-2,1], grid = True, size = [9,3])
   v = np.arange(0, 0.22, 0.02 )
   # cp = Mp.PlaneContour(fig, ax, 'Z', 0, 'v', levels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')   
   cp,cl = Mp.PlaneContour(fig, ax, 'Z', 0, 'u_std', cflevels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k') 
   # t = np.arange(-0.3,0.4,0.1) 
   t = v   
   clb = fig.colorbar(cp,ticks=t)
   clb.ax.set_title('$\mathrm{u}\' / \mathrm{U}_{\infty} $')
   cir = patches.Circle((0,0),1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(cir)   
   fig.savefig('u\'_Zr0', dpi=330)   

# %% plot4
with timer('Om_Z contour at Z/D = 0'):
      
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-1.5,6], ylim = [-2,1], grid = True, size = [9,3])
   v = np.arange(-0.4, 0.5, 0.1 )
   # cp = Mp.PlaneContour(fig, ax, 'Z', 0, 'v', levels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')   
   cp,cl = Mp.PlaneContour(fig, ax, 'Z', 0, 'Om_z', cflevels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k') 
   # t = np.arange(-0.3,0.4,0.1) 
   t = v   
   clb = fig.colorbar(cp,ticks=t)
   clb.ax.set_title('$\omega_{z} \mathrm{c} / \mathrm{U}_{\infty}$')
   cir = patches.Circle((0,0),1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(cir)   
   fig.savefig('Omz_Zr0', dpi=330)   


# %% plot5
with timer('U contour at Y/D = 0'):
   
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-1.5,6], ylim = [-0,1], grid = True, size = [9,2])
   # v = np.linspace(0.0, 1.2, 12, endpoint=True)
   v = np.arange(0.0, 1.28, 0.05 )
   # v = np.arange(0,1.4,0.2)
   cp,cl = Mp.PlaneContour(fig, ax, 'Y', 0.5, 'u', cflevels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')
   t = np.arange(0,1.4,0.2)   
   clb = fig.colorbar(cp,ticks=t)
   clb.ax.set_title('$\mathrm{U} / \mathrm{U}_{\infty} $')
   rect = patches.Rectangle((-1,-0.5),2,1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(rect)
   ax.set_ylabel('$\it{Z/H}$')
   ax.set_yticks([0,0.5,1.0])
   # fig.savefig('U_Yr0', dpi=330)

# %% plot6
with timer('V contour at Y/D = 0'):
   
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-1.5,6], ylim = [-0,1], grid = True, size = [9,2])
   # v = np.linspace(0.0, 1.2, 12, endpoint=True)
   v = np.arange(-0.25, 0.3, 0.05 )
   # v = np.arange(0,1.4,0.2)
   cp,cl = Mp.PlaneContour(fig, ax, 'Y', 0.5, 'v', cflevels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')
   t = np.arange(-0.3, 0.4, 0.1)   
   clb = fig.colorbar(cp,ticks=t)
   clb.ax.set_title('$\mathrm{V} / \mathrm{V}_{\infty} $')
   rect = patches.Rectangle((-1,-0.5),2,1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(rect)
   ax.set_ylabel('$\it{Z/H}$')
   ax.set_yticks([0,0.5,1.0])
   fig.savefig('V_Yr0', dpi=330)

# %% plot7
with timer('W contour at Y/D = 0'):
   
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-1.5,6], ylim = [-0,1], grid = True, size = [9,2])
   # v = np.linspace(0.0, 1.2, 12, endpoint=True)
   v = np.arange(-0.25, 0.3, 0.05 )
   # v = np.arange(0,1.4,0.2)
   cp,cl = Mp.PlaneContour(fig, ax, 'Y', 0.5, 'w', cflevels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')
   t = np.arange(-0.3, 0.4, 0.1)   
   clb = fig.colorbar(cp,ticks=t)
   clb.ax.set_title('$\mathrm{W} / \mathrm{U}_{\infty} $')
   rect = patches.Rectangle((-1,-0.5),2,1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(rect)
   ax.set_ylabel('$\it{Z/H}$')
   ax.set_yticks([0,0.5,1.0])
   # fig.savefig('W_Yr0', dpi=330)

# %% plot8
with timer('Omy contour at Y/D = 0'):
   
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-1.5,6], ylim = [-0,1], grid = True, size = [9,2])
   # v = np.linspace(0.0, 1.2, 12, endpoint=True)
   v = np.arange(-0.25, 0.3, 0.1 )
   # v = np.arange(0,1.4,0.2)
   cp,cl = Mp.PlaneContour(fig, ax, 'Y', 0.5, 'Om_y', cflevels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')
   t = np.arange(-0.3, 0.4, 0.1)   
   clb = fig.colorbar(cp,ticks=t)
   clb.ax.set_title('$\omega_{y} \mathrm{c} / \mathrm{U}_{\infty} $')
   rect = patches.Rectangle((-1,-0.5),2,1,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(rect)
   ax.set_ylabel('$\it{Z/H}$')
   ax.set_yticks([0,0.5,1.0])
   fig.savefig('Omy_Yr0', dpi=330)











# %% plot6
with timer('U contour at X/D = 0.5'):
      
   fig, ax = plt.subplots()
   Mp.PlotStyle(fig, ax, xlim = [-2,2], ylim = [-2,2], grid = True, size = [3.6,3])
   v = np.linspace(0.4, 1.1, 8, endpoint=True)
   cp,cl = Mp.PlaneContour(fig, ax, 'X', 1, 'u', levels=v, cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')   
   # cp = Mp.PlaneContour(fig, ax, 'X', 0, 'u', cmap = 'coolwarm', contourline=True, linewidths=0.6, colors='k')   
   clb = fig.colorbar(cp)
   # ax.set_ylabel('$\it{Y/D}$')
   clb.ax.set_title('$\mathrm{U} / \mathrm{U}_{\infty} $')
   rect = patches.Rectangle((-1,-1),2,2,linewidth=2,edgecolor='k',facecolor='none',linestyle='dashed')
   ax.add_patch(rect)



# %% Func Test

xg, yg, ExtDf = Mp.ExtractPlane('Z',0, 'intp', 100, 250)
varg = griddata((ExtDf['X'], ExtDf['Y']), ExtDf['u'], (xg, yg), method='linear')
U = ExtDf['u']
V = ExtDf['v']

