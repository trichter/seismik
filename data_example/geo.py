#!/usr/bin/env python
# (C) 2018, Tom Eulenfeld, MIT license

import seisgeo as sg
import seisutil as su
import matplotlib.pyplot as plt

names, nums, x, y, z = sg.read_coords('GEO/geo_eubabrunn_tag3.txt')
line = sg.odr_regr(x, y)
xp, yp = sg.project_coordinates(x, y, *line)
sg.plot_projection(x, y, xp, yp, *line, fname='GEO/map0.png')
sps = su.read_shotpoints('GEO/shotpoints.txt', verbose=True)
sg.write_geo(xp, yp, z, xp[0], yp[0], sps, 'GEO/rec_xz.dat', 'GEO/sou_xz.dat')

plt.show()
