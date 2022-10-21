# (C) 2018-2020, Tom Eulenfeld, MIT license

import numpy as np
import matplotlib.pyplot as plt
import re


def _f(B, x):
    """Linear function y=m*x+b"""
    return B[0]*x + B[1]


def odr_regr(x, y):
    """Orthogonal distance regression, return m and b of line y=m*x+b"""
    from scipy.odr import Data, Model, ODR
    linear = Model(_f)
    mydata = Data(x, y, wd=1, we=1)
    myodr = ODR(mydata, linear, beta0=[1., 2.])
    res = myodr.run()
    return res.beta


def project_coordinates(x, y, m, b):
    """Project coordinates x,y onto line y=mx+b"""
    # direction of line (1, m)
    c = (1 * x + m * (y - b)) / (1 + m ** 2)
    return 1 * c, m * c + b


def plot_projection(x, y, xp, yp, m, b, fname=None):
    """Plot projection

    x,y: coordinates
    xp,yp: projected coordinates
    m,b: line y=m*x+b
    """
    plt.figure()
    plt.plot(x, y, '.')
    plt.plot(x, m*np.array(x)+b)
    plt.plot(xp, yp, '.')
    plt.axis('scaled')
    if fname:
        plt.savefig(fname)


def plot_html(fname_in, fname_out, zone=None, delimiter=',', encoding=None,
              color=None, color_func=None, reverse=False, **kw):
    """Plot html of coordinates in fname_in"""
    import folium
    import utm
    data = np.genfromtxt(fname_in, dtype=None, delimiter=delimiter, encoding=encoding, **kw)
    names, x, y, z, *_ = list(zip(*reversed(data))) if reverse else zip(*data)
    x0 = np.mean(x)
    y0 = np.mean(y)
    latlon0 = utm.to_latlon(x0, y0, *zone) if zone else (y0, x0)
    m = folium.Map(location=latlon0, zoom_start=15, max_zoom=20, control_scale=True)
    m.add_child(folium.LatLngPopup())
    for name, xp, yp, zp in zip(names, x, y, z):
        if color is None and color_func is None:
            c = '#3186cc'
        elif color is None:
            c = color_func(name)
        else:
            c = color
        latlon = utm.to_latlon(xp, yp, *zone) if zone else (yp, xp)
        folium.CircleMarker(latlon, radius=1,
                            popup='%s %.5f %.5f %.2f' % (name, latlon[0], latlon[1], zp),
                            color=c, fill_color=c
                            ).add_to(m)
    m.save(fname_out)


def read_coords(fname, delimiter=',', sort=True):
    """Read coordinates in fname

    name, x, y, z
    """
    data = np.genfromtxt(fname, dtype=None, encoding=None, delimiter=delimiter)
    names, x, y, z, *_ = zip(*data)
    nums = [int(re.sub('[a-zA-Z]+', '', name)) for name in names]
    if sort:
        ind = sorted(range(len(names)), key=lambda i: nums[i])
        nums = [nums[i] for i in ind]
        names = [names[i] for i in ind]
        x = [x[i] for i in ind]
        y = [y[i] for i in ind]
        z = [z[i] for i in ind]
    return names, np.array(nums), np.array(x), np.array(y), np.array(z)


def test_project_coords():
    x =  np.arange(10)
    y = 0.5 * x + 10 + 2*np.random.random(10)
    m, b = odr_regr(x, y)
    xp, yp = project_coordinates(x, y, m, b)
    plot_projection(x, y, xp, yp, m, b)


def write_geo(x, y, z, x0, y0, sps, fname1, fname2, rec0=1, sort=False,
              xyzs=None, shot0=1):
    """Write profil files (receivers and shots)

    x,y,z: profile coordinates
    x0,y0: profile in m relative to this coordinates
    fname1,fname2: file name for receivers, shots
    sps: shotpoints - geofon dictionary from read_shotpoints
    rec0,shot0: Start number of receivers/shots
    xyzs: sps can be set to None. Then, use xyzs.
          tuple of profile coordinates of shotpoints,
          (x, ys, zs)
    """
    fmt = '%3d  %8.3f  %8.3f'
    pos = ((x - x0) ** 2 + (y - y0) ** 2) ** 0.5
    num = np.arange(len(x)) + rec0

    if xyzs is not None:
        x2, y2, z2 = xyzs
        pos2 = ((x2 - x0) ** 2 + (y2 - y0) ** 2) ** 0.5
        num2 = np.arange(len(x2)) + shot0
    elif sps is not None:
        shot_ind=np.array(list(sps.values())) - 1
        num2 = sps.keys()
        pos2 = pos[shot_ind]
        z2 = z[shot_ind]
    else:
        raise

    if sort:
        z = [a for _, a in sorted(zip(pos, z))]
        z2 = [a for _, a in sorted(zip(pos2, z2))]
        pos = sorted(pos)
        pos2 = sorted(pos2)
    np.savetxt(fname2, list(zip(num2, pos2, z2)), fmt=fmt)
    np.savetxt(fname1, list(zip(num, pos, z)), fmt=fmt)


def profilemeter2xy(meters, x, y, print_lines=False):
    dy = y[-1] - y[0]
    dx = x[-1] - x[0]
    alpha = np.arctan2(dy, dx)
    res = [(meter * np.cos(alpha) + x[0],
            meter * np.sin(alpha) + y[0]) for meter in meters]
    if print_lines:
        for i, (x2, y2), in enumerate(res):
            print(f'x{i}={x2:.3f} y{i}={y2:.3f} z{i}=',)
    return res


def interpolate_geophones(nums, x, y, z, nnums=None):
    if nnums is None:
        nnums = np.arange(nums[0], nums[-1]+1)
    res = [np.interp(nnums, nums, coordpart) for coordpart in (x, y, z)]
    return res


if __name__ == '__main__':
    test_project_coords()
