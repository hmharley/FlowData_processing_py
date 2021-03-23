import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from scipy.interpolate import griddata
import copy
# import tecplot as tp

# with open('Rectangle_EXP.dat') as Rectangle_EXP:
#     all_data = 
# D rectangle = 100

def load_data(fname):
    # To load tecplot dat to dataframe
    # fname: file name
    data = np.genfromtxt(fname, delimiter=' ', skip_header=10) 
    dataset = pd.DataFrame(data,columns=['X_D','Y_D','Z_D','U_Uinf'])
    return dataset

def D_area(D,H):
    # Calculate D_area for rectangle
    D_area = ((4/math.pi)*D*H)**0.5
    return D_area

def D_star(D,H):
    # Calculate D_star for rectangle
    D_star = 2*D*H/(D+H)
    return D_star

def wake_scaling(dataset, D_gometry, D_scaling):
    # To scale the wake. 
    # D_geometry: geometry diameter 
    # D_scaling: scale length
    df = copy.copy(dataset)
    df['X_D'] = df['X_D'] * D_gometry / D_scaling
    df['Y_D'] = df['Y_D'] * D_gometry / D_scaling
    df['Z_D'] = df['Z_D'] * D_gometry / D_scaling
    return df

def extract_plane(dataset, xd):
    # To interplate results in desired plane.
    # dataset: normalized dataframe 
    # XD: a downstream position
    df = dataset.loc[np.round(dataset['X_D'], 1)==round(xd,1)]
    y = np.arange(np.min(df['Y_D']), np.max(df['Y_D']), 0.01)
    z = np.arange(np.min(df['Z_D']), np.max(df['Z_D']), 0.01)
    yg, zg = np.meshgrid(y, z)
    u = griddata((df['Y_D'], df['Z_D']), df['U_Uinf'], (yg, zg), method='linear')
    return yg, zg, u

def extract_line(yg,zg,u, zd):
    yl = yg[np.round(zg,2)==round(zd,2)]
    ul = u[np.round(zg,2)==round(zd,2)]
    return yl,ul


# def grid_interplation(Rectangle):
#     XX, YY, ZZ = np.meshgrid(Rectangle[['X_D']],Rectangle[['Y_D']],Rectangle[['Z_D']],sparse=True)


# %% Scaling length calculation
D_cir = 200; D_squ = 200; D_rec = 100; H_rec = 300
# Area-based Shammensodin and Port-Agel
D_area_cir = D_cir
D_area_squ = D_area(D_squ,D_squ)
D_area_rec = D_area(D_rec,H_rec)
# Area and perimeter based scaling length
D_star_cir = D_cir
D_star_squ = D_star(D_squ,D_squ)
D_star_rec = D_star(D_rec,H_rec) 


# %% Read tecplot .dat file
f1 = 'Rectangle_EXP.dat'
f2 = 'Circle_EXP.dat'
f3 = 'Square_EXP.dat'

Rectangle = load_data(f1)
Circle = load_data(f2)
Square = load_data(f3)

#  wake scaling
# D: loaded already 
# %%
# D_area
Cir_area = wake_scaling(Circle,D_cir,D_area_cir)
Squ_area = wake_scaling(Square,D_squ,D_area_squ)
Rec_area = wake_scaling(Rectangle,D_rec,D_area_rec)
# %%
# D_star
Cir_star = wake_scaling(Circle,D_cir,D_star_cir)
Squ_star = wake_scaling(Square,D_squ,D_star_squ)
Rec_star = wake_scaling(Rectangle,D_rec,D_star_rec)

# %% Interpolation at postions wanted:
xd = np.linspace(1.,5.,5)
zd = np.linspace(0.5,0.,3)


# for i,x in enumerate(xd):
#     yg, zg, u = extract_plane(Cir_area,x)
#     for j,z in enumerate(zd):
#         plt.subplot(len(zd),len(xd),len(zd)*i + j+1)
#         yl,ul = extract_line(yg,zg,u, z)
#         ul = ul[yl<0.6]; yl = yl[yl<0.6]
#         ul = ul[yl>-1]; yl = yl[yl>-1]
#         plt.plot(yl,ul)

# for i,x in enumerate(xd):
#     yg, zg, u = extract_plane(Squ_area,x)
#     for j,z in enumerate(zd):
#         plt.subplot(len(zd),len(xd),len(zd)*i + j+1)
#         yl,ul = extract_line(yg,zg,u, z)
#         ul = ul[yl<1]; yl = yl[yl<1]
#         ul = ul[yl>-1]; yl = yl[yl>-1]
#         plt.plot(yl,ul)

# for i,x in enumerate(xd):
#     yg, zg, u = extract_plane(Rec_area,x)
#     for j,z in enumerate(zd):
#         plt.subplot(len(zd),len(xd),len(zd)*i + j+1)
#         yl,ul = extract_line(yg,zg,u, z)
#         ul = ul[yl<1]; yl = yl[yl<1]
#         ul = ul[yl>-1]; yl = yl[yl>-1]
#         plt.plot(yl,ul)

# plt.show()

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
# %% plot D_area
fig1, axs1 = plt.subplots(len(zd),len(xd), sharex=True, sharey=True, gridspec_kw={'hspace': 0, 'wspace': 0.1}, figsize=(30,10))
fig1.suptitle('U/U_inf')
axs1[0,0].set_ylabel('Z/Darea=0.5')
axs1[1,0].set_ylabel('Z/Darea=0.25')
axs1[2,0].set_ylabel('Z/Darea=0')
# xlabel
for i,x in enumerate(xd):
    axs1[2,i].set_xlabel('Y/Darea')
axs1[0,0].set_title('X/Darea=1')
axs1[0,1].set_title('X/Darea=2')
axs1[0,2].set_title('X/Darea=3')
axs1[0,3].set_title('X/Darea=4')
axs1[0,4].set_title('X/Darea=5')
plt.setp(axs1, xlim=(-0.99,0.99))

for i,x in enumerate(xd):
    yg, zg, u = extract_plane(Cir_area,x)
    for j,z in enumerate(zd):
        yl,ul = extract_line(yg,zg,u, z)
        ul = ul[yl<0.6]; yl = yl[yl<0.6]
        ul = ul[yl>-1]; yl = yl[yl>-1]
        axs1[j, i].plot(yl,ul,color='grey',marker='o', fillstyle='full', markevery=8,markersize=12,linestyle='None')

for i,x in enumerate(xd):
    yg, zg, u = extract_plane(Squ_area,x)
    for j,z in enumerate(zd):
        yl,ul = extract_line(yg,zg,u, z)
        ul = ul[yl<0.9]; yl = yl[yl<0.9]
        ul = ul[yl>-1]; yl = yl[yl>-1]
        axs1[j, i].plot(yl,ul,color='red',marker='s', fillstyle='none', markevery=8,markersize=12,linestyle='None',markeredgewidth = 2)

for i,x in enumerate(xd):
    yg, zg, u = extract_plane(Rec_area,x)
    for j,z in enumerate(zd):
        yl,ul = extract_line(yg,zg,u, z)
        ul = ul[yl<0.9]; yl = yl[yl<0.9]
        ul = ul[yl>-1]; yl = yl[yl>-1]
        axs1[j, i].plot(yl,ul,color='blue',marker='d', fillstyle='none', markevery=8,markersize=12,linestyle='None',markeredgewidth = 2)

# %% plot D_star
fig2, axs2 = plt.subplots(len(zd),len(xd), sharex=True, sharey=True, gridspec_kw={'hspace': 0, 'wspace': 0.1}, figsize=(30,10))
fig1.suptitle('U/U_inf')
axs2[0,0].set_ylabel('Z/Dstar=0.5')
axs2[1,0].set_ylabel('Z/Dstar=0.25')
axs2[2,0].set_ylabel('Z/Dstar=0')
# xlabel
for i,x in enumerate(xd):
    axs2[2,i].set_xlabel('Y/Dstar')
axs2[0,0].set_title('X/Dstar=1')
axs2[0,1].set_title('X/Dstar=2')
axs2[0,2].set_title('X/Dstar=3')
axs2[0,3].set_title('X/Dstar=4')
axs2[0,4].set_title('X/Dstar=5')
plt.setp(axs2, xlim=(-0.99,0.99))

for i,x in enumerate(xd):
    yg, zg, u = extract_plane(Cir_star,x)
    for j,z in enumerate(zd):
        yl,ul = extract_line(yg,zg,u, z)
        ul = ul[yl<0.6]; yl = yl[yl<0.6]
        ul = ul[yl>-1]; yl = yl[yl>-1]
        axs2[j, i].plot(yl,ul,color='grey',marker='o', fillstyle='full', markevery=8,markersize=12,linestyle='None')

for i,x in enumerate(xd):
    yg, zg, u = extract_plane(Squ_star,x)
    for j,z in enumerate(zd):
        yl,ul = extract_line(yg,zg,u, z)
        ul = ul[yl<0.9]; yl = yl[yl<0.9]
        ul = ul[yl>-1]; yl = yl[yl>-1]
        axs2[j, i].plot(yl,ul,color='red',marker='s', fillstyle='none', markevery=8,markersize=12,linestyle='None',markeredgewidth = 2)

for i,x in enumerate(xd):
    yg, zg, u = extract_plane(Rec_star,x)
    for j,z in enumerate(zd):
        yl,ul = extract_line(yg,zg,u, z)
        ul = ul[yl<0.9]; yl = yl[yl<0.9]
        ul = ul[yl>-1]; yl = yl[yl>-1]
        axs2[j, i].plot(yl,ul,color='blue',marker='d', fillstyle='none', markevery=8,markersize=12,linestyle='None',markeredgewidth = 2)



# %% save fig
fig1.savefig('U_Darea.svg', dip = 300)
fig2.savefig('U_Dstar.svg', dip = 300)


# %%
