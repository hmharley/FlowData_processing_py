import sys
# sys.path.append("D:/Ming_DOCUMENTS/Working/FlowData_processing_py/source/")
sys.path.append("F:/HM/SURFDRIVE/Working/FlowData_processing_py/source/")
# sys.path.append("G:/SURFDRIVE/Working/FlowData_processing_py/source/")
import math
import pandas as pd
import plt2pandas as p2p
import numpy as np
# from findiff import FinDiff, coefficients, Coefficient
from timer_post import timer
from time import time
from MyPlot import MyPlot
from scipy.interpolate import griddata
from matplotlib.ticker import AutoMinorLocator
import variable_analysis as VA
import matplotlib.pyplot as plt
import matplotlib.patches as patches


# fpath = 'Z:\\scratch\\S_VAWT\\DoubleP10\\postProcessing\\sampleDict\\1.20044579502\\'
# fpath = 'Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\SmallVAWT\\InlineCases\\DoubleBaseline\\postProcessing\\sampleDict\\1.20080847147\\'
# fpath = 'Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\SmallVAWT\\PitchedVAWT\\p2\\postProcessing\\sampleDict\\1.800226404308\\'
# fpath = 'Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\SmallVAWT\\PitchedVAWT\\p5\\postProcessing\\sampleDict\\1.601114575433\\'
# fpath = 'Z:\\scratch\\S_VAWT\\DoubleM10\\postProcessing\\sampleDict\\0.8009353457\\'

fpath = 'F:\\HM\\SURFDRIVE\\Working\\'
# casename = 'Baseline'
# fpath = 'Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\SmallVAWT\\PitchedVAWT\\Polar_NACA0012_2e4\\M10_coarse\\postProcessing\\sampleDict\\1.4\\'
# casename = 'M10'
# fpath = 'Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\SmallVAWT\\PitchedVAWT\\Polar_NACA0012_2e4\\P10_coarse\\postProcessing\\sampleDict\\0.8\\'
# fpath = 'Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\SmallVAWT\\PitchedVAWT\\Polar_NACA0012_2e4\\P10Test\\postProcessing\\sampleDict\\1.0006199404\\'
# casename = 'P10'
# fpath = 'Z:\\OpenFOAM\\minghuang-4.1\\run\\AL_VAWT\\SmallVAWT\\PitchedVAWT\\Polar_NACA0018\\M10\\postProcessing\\sampleDict\\1.6016609707864\\'
casename = 'M10NACA0018TKE1e3'
saveplot = 0
VectorPlot = 1
StreamLinePlot = 0
U_inf = 5
D = 0.3

# %% Script Xplane
plane = 'UMean_planeX15'
fname = fpath + plane + '.raw'
varlist = ['x','y','z','Ux','Uy','Uz']

size = [5,6]
clevel = np.linspace(0.3,1.2,10)
clblevel = np.linspace(0.5,1.1,7)
cmap = 'jet'
contourline=True
xlim = [-2,2]
ylim = [-2,2]


plot1 = MyPlot()
plot1.LoadRaw(fname,varlist)
df = plot1.flowdata
y = np.arange(np.min(df['y']), np.max(df['y']), 0.01)/D
z = np.arange(np.min(df['z']), np.max(df['z']), 0.01)/D
yg, zg = np.meshgrid(y, z)
ux = griddata((df['y']/D, df['z']/D), df['Ux'], (yg, zg), method='cubic')
ux = ux/U_inf
fig, ax = plt.subplots()

cp = ax.contourf( yg, zg, ux,clevel, extend='both')
# cp = ax.contourf( yg, zg, ux)
clevel2 = np.delete(clevel,5)
cpl = ax.contour( yg, zg, ux,clevel2,colors='k', linestyles = 'solid', linewidths = 1)
cl1 = ax.contour(yg, zg, ux, [0.8], linestyles='dashed', colors=('k'),linewidths=(2)   ) 
fig.set_size_inches(size[0], size[1])

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r'\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother'
plt.rc('grid', color='black',linewidth=0.8,alpha=.5)
# plt.set_cmap('jet')
plt.set_cmap(cmap)
# ax.set(xlim=xlim, ylim=ylim, aspect=1)
ax.set(xlim=xlim, ylim=ylim, aspect=1)
# ax.set_xlim(-1, 1)
# ax.set_ylim(-1, 1)
ax.set_xlabel('$Y/D$', fontsize=12)
ax.set_ylabel('$Z/D$', fontsize=12)
VAWT = patches.Rectangle((-0.5,-0.5),1,1,
                        linewidth=2,edgecolor='darkgray',facecolor='none')
Pole = patches.Rectangle((-0.015,-1),0.03,2,
                        linewidth=0.1,edgecolor='darkgray',facecolor='darkgray')                       
ax.add_patch(VAWT)
ax.add_patch(Pole)


clb = fig.colorbar(cp,extendrect = True)
cax = clb.ax
cax.hlines(0.556, 0, 1, colors = 'k', linewidth = 2, linestyles ='dashed')
clb.ax.set_title('$U_{x}/U_{\infty}$')
# ax.set_title('$X/D=1$')
# ax.set_xticks(np.linspace(-1,1,5))
# ax.set_yticks(np.linspace(-1,1,5))
ax.set_xticks(ax.get_xticks()[1:])
ax.set_title('$X/D = 15 $', fontsize=12)
# ax[i].set_xlabel('$Y/D$', fontsize=12)
# ax[i].set_ylabel('$Z/D$', fontsize=12)
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)


if VectorPlot:  
    uy = griddata((df['y']/D, df['z']/D), df['Uy'],(yg, zg),method='cubic')
    uz = griddata((df['y']/D, df['z']/D), df['Uz'],(yg, zg),method='cubic') 
    uy = uy/U_inf 
    uz = uz/U_inf     
    skip = 3
    Q = ax.quiver(yg[::skip, ::skip], zg[::skip,::skip], uy[::skip, ::skip], uz[::skip,::skip] ,
               pivot='mid', units='inches',linewidths = 0.5, edgecolors='k')

# if StreamLinePlot:
#     ux = griddata((df['x'], df['y']), df['Vx'],(xg, yg),method='cubic')
#     uy = griddata((df['x'], df['y']), df['Vy'],(xg, yg),method='cubic')
#     isValid = griddata((Mp.flowdata['x'], Mpflowdata['y']), Mp.flowdata['isValid'], (xg, yg,  method='cubic')
#     ux[isValid<0.01] = None
#     uy[isValid<0.01] = None
#     ux = ux/U_inf   
#     uy = uy/U_inf 
#     skip = 8    
#     # stream_points = np.array(zip(np.arange(-9,9,5), -np.arange(-9,9,.5)))
#     SL = ax[i].streamplot(xg, yg, ux, uy, color ='k',density=2)



# %% Script Yplane
plane = 'UMean_planeY'
fname = fpath + plane + '.raw'
varlist = ['x','y','z','Ux','Uy','Uz']

size = [15,2.5]
clevel = np.linspace(0.3,1.2,10)
clblevel = np.linspace(0.5,1.1,7)
cmap = 'jet'
contourline=True
xlim = [-1,15]
ylim = [-1.5,1.5]


plot1 = MyPlot()
plot1.LoadRaw(fname,varlist)
df = plot1.flowdata
z = np.arange(np.min(df['z']), np.max(df['z']), 0.01)/D
x = np.arange(np.min(df['x']), np.max(df['x']), 0.01)/D
xg, zg = np.meshgrid(x, z)
ux = griddata((df['x']/D, df['z']/D), df['Ux'], (xg, zg), method='cubic')
ux = ux/U_inf
fig, ax = plt.subplots()

cp = ax.contourf( xg, zg, ux,clevel, extend='both')
# cp = ax.contourf( yg, zg, ux)
clevel2 = np.delete(clevel,5)
cpl = ax.contour( xg, zg, ux,clevel2,colors='k', linestyles = 'solid', linewidths = 1)
cl1 = ax.contour(xg, zg, ux, [0.8], linestyles='dashed', colors=('k'),linewidths=(2)   ) 
fig.set_size_inches(size[0], size[1])

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r'\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother'
plt.rc('grid', color='black',linewidth=0.8,alpha=.5)
# plt.set_cmap('jet')
plt.set_cmap(cmap)
ax.set(xlim=xlim, ylim=ylim, aspect=1)
# ax.set_xlim(-1, 1)
# ax.set_ylim(-1, 1)
ax.set_xlabel('$X/D$', fontsize=12)
ax.set_ylabel('$Z/D$', fontsize=12)
VAWT1 = patches.Rectangle((-0.5,-0.5),1,1,
                        linewidth=2,edgecolor='darkgray',facecolor='none')
Pole1 = patches.Rectangle((-0.015,-1),0.03,2,
                        linewidth=0.1,edgecolor='darkgray',facecolor='darkgray')     

VAWT2 = patches.Rectangle((4.5,-0.5),1,1,
                        linewidth=2,edgecolor='darkgray',facecolor='none')
Pole2 = patches.Rectangle((4.985,-1),0.03,2,
                        linewidth=0.1,edgecolor='darkgray',facecolor='darkgray')                                              
ax.add_patch(VAWT1)
ax.add_patch(VAWT2)
ax.add_patch(Pole1)
ax.add_patch(Pole2)


clb = fig.colorbar(cp,extendrect = True)
cax = clb.ax
cax.hlines(0.556, 0, 1, colors = 'k', linewidth = 2, linestyles ='dashed')
clb.ax.set_title('$U_{x}/U_{\infty}$')
# ax.set_title('$X/D=1$')
# ax.set_xticks(np.linspace(-1,1,5))
# ax.set_yticks(np.linspace(-1,1,5))
# ax.set_xticks(ax.get_xticks()[1:])
# ax.set_title('$X/D = 1 $', fontsize=12)
# ax[i].set_xlabel('$Y/D$', fontsize=12)
# ax[i].set_ylabel('$Z/D$', fontsize=12)
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)


# if VectorPlot:  
#     uy = griddata((df['x']/D, df['z']/D), df['Uy'],(xg, zg),method='cubic')
#     uz = griddata((df['x']/D, df['z']/D), df['Uz'],(xg, zg),method='cubic') 
#     uy = uy/U_inf 
#     uz = uz/U_inf     
#     skip = 3
#     Q = ax.quiver(xg[::skip, ::skip], zg[::skip,::skip], uy[::skip, ::skip], uz[::skip,::skip] ,
#                pivot='mid', units='inches',linewidths = 0.5, edgecolors='k')

# %% Plane Z
plane = 'UMean_planeZ'
fname = fpath + plane + '.raw'
varlist = ['x','y','z','Ux','Uy','Uz']

size = [15,2.5]
clevel = np.linspace(0.3,1.2,10)
clblevel = np.linspace(0.5,1.1,7)
cmap = 'jet'
contourline=True
xlim = [-1,15]
ylim = [-3.5,2]


plot1 = MyPlot()
plot1.LoadRaw(fname,varlist)
df = plot1.flowdata
y = np.arange(np.min(df['y']), np.max(df['y']), 0.01)/D
x = np.arange(np.min(df['x']), np.max(df['x']), 0.01)/D
xg, yg = np.meshgrid(x, y)
ux = griddata((df['x']/D, df['y']/D), df['Ux'], (xg, yg), method='cubic')
ux = ux/U_inf
fig, ax = plt.subplots()

cp = ax.contourf( xg, yg, ux,clevel, extend='both')
# cp = ax.contourf( yg, zg, ux)
clevel2 = np.delete(clevel,5)
cpl = ax.contour( xg, yg, ux,clevel2,colors='k', linestyles = 'solid', linewidths = 1)
cl1 = ax.contour(xg, yg, ux, [0.8], linestyles='dashed', colors=('k'),linewidths=(2)   ) 
fig.set_size_inches(size[0], size[1])

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r'\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother'
plt.rc('grid', color='black',linewidth=0.8,alpha=.5)
# plt.set_cmap('jet')
plt.set_cmap(cmap)
ax.set(xlim=xlim, ylim=ylim, aspect=1)

# ax.set_xlim(-1, 1)
# ax.set_ylim(-1, 1)
ax.set_xlabel('$X/D$', fontsize=12)
ax.set_ylabel('$Y/D$', fontsize=12)
VAWT1 = patches.Circle((0,0),0.5,
                        linewidth=2,edgecolor='darkgray',facecolor='none')
VAWT2 = patches.Circle((5,0),0.5,
                        linewidth=2,edgecolor='darkgray',facecolor='none')                        
# Pole = patches.Rectangle((-0.015,-1),0.03,2,
                        # linewidth=0.1,edgecolor='darkgray',facecolor='darkgray')                       
ax.add_patch(VAWT1)
ax.add_patch(VAWT2)


clb = fig.colorbar(cp,extendrect = True)
cax = clb.ax
cax.hlines(0.556, 0, 1, colors = 'k', linewidth = 2, linestyles ='dashed')
clb.ax.set_title('$U_{x}/U_{\infty}$')
# ax.set_title('$X/D=1$')
# ax.set_xticks(np.linspace(-1,1,5))
# ax.set_yticks(np.linspace(-1,1,5))
# ax.set_xticks(ax.get_xticks()[1:])
# ax.set_title('$X/D = 1 $', fontsize=12)
# ax[i].set_xlabel('$Y/D$', fontsize=12)
# ax[i].set_ylabel('$Z/D$', fontsize=12)
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)


# if VectorPlot:  
#     uy = griddata((df['x']/D, df['z']/D), df['Uy'],(xg, zg),method='cubic')
#     uz = griddata((df['x']/D, df['z']/D), df['Uz'],(xg, zg),method='cubic') 
#     uy = uy/U_inf 
#     uz = uz/U_inf     
#     skip = 3
#     Q = ax.quiver(xg[::skip, ::skip], zg[::skip,::skip], uy[::skip, ::skip], uz[::skip,::skip] ,
#                pivot='mid', units='inches',linewidths = 0.5, edgecolors='k')
if saveplot:
    plt.savefig(casename + '_' + plane + '.png', dpi=150)
# %%
