import sys
# sys.path.append("D:/Ming_DOCUMENTS/Working/Data_processing_py/source/")
# sys.path.append("F:/HM/SURFDRIVE/Working/Data_processing_py/source/")
import plt2pandas as p2p
from timer_post import timer
from time import time
from MyPlot import MyPlot
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.set_cmap('coolwarm')

# %% Class method
FoldPath = "plotTest/"
OutFile  = "plotTest/"
# VarList  = ['x [mm]', 'y [mm]', 'z [mm]', 'Vx [m/s]', 'Vy [m/s]', 'Vz [m/s]', 'isValid']
VarList  = ['X','Y','Z','u']
with timer("Load Meanflow data"):
   Mp = MyPlot()
   Mp.LoadPlt(FoldPath, VarList)
   
with timer('Normalization'):
   Mp.CoordNorm(300.0)


# %%
with timer("Plot Ux on Z/D = 0"):
   # xg, yg = Mp.extract_plane('Z', 0.0,intp='intp', m=200, n=200) 
   xg, yg, ExtDf = Mp.ExtractPlane('Z', 0)    # extract Z plane 
   # Mp.flowdata = Mp.flowdata.fillna(0)
   Ug = griddata((ExtDf['X'], ExtDf['Y']), ExtDf['u'], (xg, yg), method='linear')   
   fig1,ax1 = plt.subplots()
   plt.contourf(xg, yg, Ug, 11)
   plt.colorbar()
   plt.show

# %%
with timer("Plot Ux on Y/D = 0"):
   # xg, yg = Mp.extract_plane('Y', 0.0,intp='intp', m=200, n=200) 
   xg, yg, ExtDf = Mp.ExtractPlane('Y', 0)    # extract Z plane 
   # Mp.flowdata = Mp.flowdata.fillna(0)

   Ug = griddata((ExtDf['X'], ExtDf['Z']), ExtDf['u'], (xg, yg), method='linear')   
   plt.contourf(xg, yg, Ug, 11)
   plt.colorbar()
   plt.show







# %% Read plt data

# FoldPath = "plotTest/"
# OutFile  = "plotTest/"
# # VarList  = ['x [mm]', 'y [mm]', 'z [mm]', 'Vx [m/s]', 'Vy [m/s]', 'Vz [m/s]', 'isValid']
# VarList  = ['X','Y','Z','u']
# with timer("Read Meanflow data"):
#    stat = p2p.ReadAllINCAResults(FoldPath, VarList)
#    # stat = p2p.ReadAllINCAResults(FoldPath)

# %% Plot
# with timer("Plot"):
#    Ucontour_Z0 = MyPlot()
#    xg, yg = Ucontour_Z0.extract_plane(stat, 'Z', 0.0)  
#    Ug = griddata((stat['X'], stat['Y']), stat['u'], (xg, yg), method='linear')   
#    plt.contourf(xg, yg, Ug, 11)
#    plt.colorbar()
#    plt.show


# %% write data file
# stat.to_hdf(OutFile+"Meanflow" + ".h5", 'w', format='fixed')
# stat.to_csv(OutFile+"Meanflow.dat", index=False, sep = '\t')