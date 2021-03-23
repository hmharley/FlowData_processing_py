class MyPlot:
    """ This is a class for visualizing the EXP/numerical data.

    """
    # Class Attribute

    # Initializer
    def __init__(self):
        # self.data = []
        self.flowdata = []
        self.postdata = []
        # self.data = []
    # Instance method: LoadRaw
    def LoadRaw(self, filename, varlist):
        """Load *.raw* file generated from OpenFoam 
        
        :param filename: [File name to load.]
        :type filename: [file, str, pathlib.Path, list of str]
        :param varlist: [Variable list in the laoded data.]
        :type varlist: [Index or array-like]
        """        
        import pandas as pd
        import numpy as np
        data = np.genfromtxt(filename, delimiter=' ', comments='#')
        self.flowdata = pd.DataFrame(data,columns = varlist)

    def LoadPlt(self, FoldPath, VarList, SavePath=None, Equ=None,
                           FileName=None, SpanAve=None, OutFile=None):
        """Load *.plt* file in a folder.

        :param FoldPath: [description]
        :type FoldPath: [type]
        :param VarList: [description]
        :type VarList: [type]
        :param SavePath: [description], defaults to None
        :type SavePath: [type], optional
        :param Equ: [description], defaults to None
        :type Equ: [type], optional
        :param FileName: [description], defaults to None
        :type FileName: [type], optional
        :param SpanAve: [description], defaults to None
        :type SpanAve: [type], optional
        :param OutFile: [description], defaults to None
        :type OutFile: [type], optional
        """        
        import tecplot as tp
        import pandas as pd
        import os
        import sys
        import warnings
        import numpy as np
        from scipy.interpolate import griddata
        from timer_post import timer
        from time import time
        import logging as log
        from glob import glob
        
        if FileName is None:
            # files = os.listdir(FoldPath)
            # FileName = [os.path.join(FoldPath, name) for name in files]
            FileName = glob(FoldPath + '*plt')
        if(isinstance(FileName, list)):
            szplt = FileName[0].find('szplt')
        else:
            szplt = FileName.find('szplt')
        if(szplt != -1):
            dataset = tp.data.load_tecplot_szl(FileName, read_data_option=2)
        else:
            dataset = tp.data.load_tecplot(FileName, read_data_option=2)
        if Equ is not None:
            for i in range(np.size(Equ)):
                tp.data.operate.execute_equation(Equ[i])
        if VarList is None:
            VarList = [v.name for v in dataset.variables()]
        # else:        
        # VarList = [v.name for v in dataset.variables()]
        df = pd.DataFrame(columns=VarList)
        if (np.size(dataset.solution_times) == 0):
            SolTime = 0.0
        else:
            SolTime = dataset.solution_times[0]
        for zone in dataset.zones('*'):
            for i in range(np.size(VarList)):
                if i == 0:
                    VarCol = zone.values(VarList[i]).as_numpy_array()
                else:
                    Var_index = zone.values(VarList[i]).as_numpy_array()
                    VarCol = np.column_stack((VarCol, Var_index))
            df1 = pd.DataFrame(data=VarCol, columns=VarList)
            df = df.append(df1, ignore_index=True)
        self.flowdata=df

        del FileName, dataset, zone
        # df = df.drop_duplicates(keep='last')
        if SpanAve is not None:
            grouped = df.groupby(['x', 'y'])
            df = grouped.mean().reset_index()
            # df = df.loc[df['z'] == 0.0].reset_index(drop=True)
        if SavePath is not None and OutFile is not None:
            st = "%08.2f" % SolTime
            df.to_hdf(SavePath+OutFile+'_'+st+".h5", 'w', format='fixed')
        self.flowdata=df

    def ExtractPlane(self, plane, planeValue, intp='intp', m=100, n=100):
        """This fuction returns the mesh grid of desired plane, which is used for interpolation of the .plt file.

        :param plane: 'X','Y','Z'
        :type plane: string
        :param planeValue: value of the plane
        :type planeValue: float64
        :param m: number of grid on horizontal axis
        :type m: int
        :param n: number of grid on vertical axis
        :type n: int        
        :return: grid data of extracted plane
        :rtype: meshgrid
        """        
        import numpy as np
        import pandas as pd
        # XD: a downstream position
        planeValue = float(planeValue)
        if plane == 'X':
            df = self.flowdata.loc[np.round(self.flowdata['X'], 1)==round(planeValue,1)]
            if intp == 'intp':
                x = np.linspace(np.min(df['Y']), np.max(df['Y']), m)
                y = np.linspace(np.min(df['Z']), np.max(df['Z']), n)
            else:
                x = df['Y']
                y = df['Z']
        elif plane == 'Y':
            df = self.flowdata.loc[np.round(self.flowdata['Y'], 1)==round(planeValue,1)]
            if intp == 'intp':
                x = np.linspace(np.min(df['X']), np.max(df['X']), m)
                y = np.linspace(np.min(df['Z']), np.max(df['Z']), n)
            else:
                x = df['X']
                y = df['Z']            
        else:
            df = self.flowdata.loc[np.round(self.flowdata['Z'], 1)==round(planeValue,1)]
            if intp == 'intp':
                x = np.linspace(np.min(df['X']), np.max(df['X']), m)
                y = np.linspace(np.min(df['Y']), np.max(df['Y']), n)
            else:
                x = df['X']
                y = df['Y']    
        xg, yg = np.meshgrid(x, y)
        return xg, yg, df                
    
    def CoordNorm(self, D=1.0):
        """Coordinate Normalization

        :param D: characterization length
        :type D: float 64
        """        
        self.flowdata['X'] = self.flowdata['X']/D
        self.flowdata['Y'] = self.flowdata['Y']/D
        self.flowdata['Z'] = self.flowdata['Z']/D

    def VarNorm(self, var, CharacterizationLenth):
        """Variable Normalization

        :param var: Variable to be normalized
        :type var: str
        :param CharacterizingLenth: characterization length
        :type CharacterizingLenth: float
        """        
        self.flowdata[var] = self.flowdata[var]/CharacterizationLenth

    def ExtractLine(self, xg, yg, v, l, tp = 'h',rd = 4):
        """Extract line in a plane.

        :param xg: gridded data horizontal
        :type xg: array
        :param yg: gridded data vertical
        :type yg: array
        :param v: [description]
        :type v: extracted variable
        :param l: locator
        :type l: float64
        :param tp: extraction type: 'v' vertical, 'h' horizontal
        :type tp: string        
        :return: xl,vl
        :rtype: float64
        """        
        import numpy as np
        if tp == 'h':
            xl = xg[np.round(yg,rd)==round(l,rd)]
            vl = v[np.round(yg,rd)==round(l,rd)]
        elif tp == 'v':
            xl = yg[np.round(xg,rd)==round(l,rd)]
            vl = v[np.round(xg,rd)==round(l,rd)]
        return xl,vl
        # Instance method: Need modification

    def PlaneContour(self, fig, ax, plane, planeValue, var, intp='intp', m=100, n=100, *, axislabel=1, cflevels=10, cllevels=None, contourline=False, linewidths=1.5, linestyles='solid', colors='k', vector=False, **param_dict):
        """Plane Contour plotter.

        :param fig: Figure for plot.
        :type fig: Figure
        :param ax: Axes for plot.
        :type ax: AxeSubplot
        :param plane: plane for plot: 'X','Y','Z'
        :type plane: str
        :param planeValue: plane value
        :type planeValue: float
        :param var: plotted variable name
        :type var: str
        :param intp: interpolation flag, defaults to 'intp', interpolated linearly. If assigned to other values, no interpolation will be implemented.
        :type intp: str, optional
        :param m: number of interpolation points on horizontal axis, defaults to 100
        :type m: int, optional
        :param n: number of interpolation points on vertical axis, defaults to 100
        :type n: int, optional
        :param levels: Determines the number and positions of the contour lines / regions, defaults to 10
        :type levels: int or array-like, optional
        :return: cp, contourf set
        :rtype: QuadContourSet
        :return: cl, contour set
        :rtype: QuadContourSet        
        """        
        from scipy.interpolate import griddata
        import matplotlib.pyplot as plt
        from matplotlib.ticker import AutoMinorLocator
        import numpy as np
        
        plt.rc('text', usetex=True)
        xg, yg, ExtDf = self.ExtractPlane(plane,planeValue, intp, m, n)
        
        # Plot 
           
        if plane == 'X':
            x = ExtDf['Y']
            y = ExtDf['Z']
        elif plane == 'Y':
            x = ExtDf['X']
            y = ExtDf['Z']                       
        else:
            x = ExtDf['X']
            y = ExtDf['Y']
        varg = griddata((x, y), ExtDf[var], (xg, yg), method='linear')
        cp = ax.contourf(xg, yg, varg, cflevels, **param_dict) 
        # Plot axis label 
        if axislabel !=None:
            if plane == 'X':
                ax.set_xlabel('$\it{Y/R}$')
                ax.set_ylabel('$\it{Z/R}$')
            elif plane == 'Y':
                ax.set_xlabel('$\it{X/R}$')
                ax.set_ylabel('$\it{Z/R}$')
            else:
                ax.set_xlabel('$\it{X/R}$')
                ax.set_ylabel('$\it{Y/R}$')                
        if contourline != False:
            if cllevels is None:
                cllevels = cflevels
            cl = ax.contour(xg, yg, varg, cllevels, linewidths=linewidths, linestyles=linestyles, colors=colors)
        if vector != False: # For now only support in plane velocity vectors
            if plane == 'X':
                U = ExtDf['v']
                V = ExtDf['w']
            elif plane == 'Y':
                U = ExtDf['u']
                V = ExtDf['w']
            else:
                U = ExtDf['u']
                V = ExtDf['v']
            Q = ax.quiver(xg, yg, U, V, units='width')
        if contourline != False:               
            return cp, cl
        return cp

    def PlotStyle(self, fig, ax, *, xlim=None, ylim=None, size=None, grid=False):
        import matplotlib.pyplot as plt
        from matplotlib.ticker import AutoMinorLocator
        # fig.colorbar(cp)
        if xlim != None:
            ax.set_xlim(xlim)
        if ylim != None:
            ax.set_ylim(ylim)
        if size != None:
            fig.set_size_inches(size[0], size[1])
        if grid != False:
            ax.grid(True)
        ax.tick_params(direction='in')


        






