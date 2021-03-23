import sys
# -------- Append the source file directory -------------
sys.path.append("D:/Ming_DOCUMENTS/Working/FlowData_processing_py/source/")
# sys.path.append("F:/HM/SURFDRIVE/Working/FlowData_processing_py/source/")
#---------------------------------------------
import numpy as np
import pandas as pd
import os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from MyPlot import MyPlot
plt.rc('text', usetex=True)
plt.set_cmap('coolwarm')
# from scipy.interpolate import griddata bv


# %% Script Xplane
# Directory of the .raw file
fname = 'D:\\Ming_DOCUMENTS\\Working\\FlowData_processing_py\\run\\AL_DS_SmallVAWT_TSR1p5\\DATA\\U_planeX5.raw'
varlist = ['x','y','z','Ux','Uy','Uz']

plot1 = MyPlot()
plot1.LoadRaw(fname,varlist)
df = plot1.flowdata
y = np.arange(np.min(df['y']), np.max(df['y']), 0.01)
z = np.arange(np.min(df['z']), np.max(df['z']), 0.01)
yg, zg = np.meshgrid(y, z)
ux = griddata((df['y'], df['z']), df['Ux'], (yg, zg), method='cubic')
plt.contourf(yg, zg, ux,11)


# plt.xlim(-1, 1)
# plt.ylim(0, 2)
plt.axis('equal')
plt.axis([-1, 1, -1, 3])
plt.colorbar()

plt.show


# %% residual
fname = 'UxFinalRes_0'
pt = 'CT0p87_corse/logs/'
varlist = ['t','Res']
plot = MyPlot()
plot1.LoadLog(pt+fname,varlist)

# %%
