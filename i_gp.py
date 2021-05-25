# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Autorki: Julia Urbanek, Julia Baranowska

import geopandas
import numpy as np
import matplotlib.pyplot as plt

gdf = geopandas.read_file('pd.shp')
gdf = gdf.to_crs('EPSG:4326')
gdf.plot('TOT', legend = True)
gdf['centroid'] = gdf.centroid

# %%

import shapely
#
xmin, ymin, xmax, ymax = [13, 48, 25, 56]
#
n_cells = 30
cell_size = (xmax - xmin)/n_cells
#
grid_cells = []

# %%

for x0 in np.arange(xmin, xmax + cell_size, cell_size):
    for y0 in np.arange(ymin, ymax + cell_size, cell_size):
        # bounds
        x1 = x0 - cell_size
        y1 = y0 + cell_size
        grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))
cell = geopandas.GeoDataFrame(grid_cells, columns = ['geometry'])

# %%

ax = gdf.plot(markersize = .1, figsize = (12, 8), column = 'TOT', cmap = 'jet')

plt.autoscale(False)
cell.plot(ax = ax, facecolor = 'none', edgecolor = 'grey')
ax.axis('off')

# %%

merged = geopandas.sjoin(gdf, cell, how = 'left', op = 'within')
dissolve = merged.dissolve(by = 'index_right', aggfunc = 'sum')
cell.loc[dissolve.index, 'TOT'] = dissolve.TOT.values

ax = cell.plot(column = 'TOT', figsize = (12, 8), cmap = 'viridis', vmax = 700000, edgecolor = 'grey', legend = True)
plt.autoscale(False)
ax.set_axis_off()
plt.axis('equal')
plt.title('liczba ludności w siatce')

# %%

a = sum(gdf.TOT_0_14)
b = sum(gdf.TOT_15_64)
c = sum(gdf.TOT_65__)
d = sum(gdf.MALE_0_14 + gdf.MALE_15_64 + gdf.MALE_65__)
e = sum(gdf.FEM_0_14 + gdf.FEM_15_64 + gdf.FEM_65__)

print('Przedział wiekowy 0-14:', a, '\n'
      'Przedział wiekowy 15-64:', b, '\n'
      'Przedział wiekowy >65:', c, '\n'
      'Ludność męska w przedziałach wiekowych a-c:', d, '\n'
      'Ludność żeńska w przedziałach wiekowych a-c:', e)

# nowy plik shp z powierzchniami województw
gdf2 = geopandas.read_file('wn')

# powierzchnie województw obliczone w odwzorowaniu WGS 84 / World Mercator
slaskie = (gdf2.area[0])/1e+6
opolskie = (gdf2.area[1])/1e+6
swietokrzyskie = (gdf2.area[2])/1e+6
lodzkie = (gdf2.area[11])/1e+6
mazowieckie = (gdf2.area[12])/1e+6
kujawsko_pomorskie = (gdf2.area[13])/1e+6
lubelskie = (gdf2.area[14])/1e+6
wielkopolskie = (gdf2.area[7])/1e+6
podkarpackie = (gdf2.area[8])/1e+6
malopolskie = (gdf2.area[9])/1e+6
warminsko_mazurskie = (gdf2.area[10])/1e+6
pomorskie = (gdf2.area[3])/1e+6
podlaskie = (gdf2.area[4])/1e+6
zachodniopomorskie = (gdf2.area[5])/1e+6
lubuskie = (gdf2.area[15])/1e+6
dolnoslaskie = (gdf2.area[6])/1e+6

l_l = sum(gdf.TOT)   # liczba ludnosci

f = [l_l/slaskie,
     l_l/opolskie,
     l_l/swietokrzyskie,
     l_l/lodzkie,
     l_l/mazowieckie,
     l_l/kujawsko_pomorskie,
     l_l/lubelskie,
     l_l/wielkopolskie,
     l_l/podkarpackie,
     l_l/malopolskie,
     l_l/warminsko_mazurskie,
     l_l/pomorskie,
     l_l/podlaskie,
     l_l/zachodniopomorskie,
     l_l/lubuskie,
     l_l/dolnoslaskie]

print('Ratio liczby ludności do powierzchni dla województwa: \n'
      'śląskiego:', int(f[0]), '\n'
      'opolskiego:', int(f[1]), '\n'
      'świętokrzyskiego:', int(f[2]), '\n'
      'łódzkiego:', int(f[3]), '\n'
      'mazowieckiego:', int(f[4]), '\n'
      'kujawsko-pomorskiego:', int(f[5]), '\n'
      'lubelskiego:', int(f[6]), '\n'
      'wielkopolskiego:', int(f[7]), '\n'
      'podkarpackiego:', int(f[8]), '\n'
      'małopolskiego:', int(f[9]), '\n'
      'warmińsko-mazurskiego:', int(f[10]), '\n'
      'pomorskiego:', int(f[11]), '\n'
      'podlaskiego:', int(f[12]), '\n'
      'zachodniopomorsiego:', int(f[13]), '\n'
      'lubuskiego:', int(f[14]), '\n'
      'dolnośląskiego:', int(f[15]))


