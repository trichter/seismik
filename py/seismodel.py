# Copyright 2022 Tom Eulenfeld, MIT license
"""
Possible invocations:
    seismodel  (needs seisconf.json)


    -n    plot negative offsets
    -v    verbose mode
    -h    show help

version 2022.08 first version

Upon start set terminal to "always on top" and enlarge matplotlib figure.
Press h for help on hot keys.
"""

import json
import sys
import matplotlib.pyplot as plt
import numpy as np

import seisutil as su
from seispick import _try_read_file, load_picks, calc_backshots

_key_doc = """
h: print help for keys
c: clear shot selction
1-6: fit n-D model
t: enter mode to set travel time curve manually,
   second time: fit model
r: enter mode to set travel time legs manually,
   second time: fit model
x: replot, clear axis from tt picks, clear shot selection
I: start IPython session
   *Danger zone* - program must be killed afterwards

Click marker for info on picks and to select specific shots
"""


def convert_picks(picks, stuff, sc):
    ch2rec = stuff.ch2rec

    def spreadch2ch(sp, ch):
        return ch - 1 + stuff.spread_from_sp[sp]['channels'][0]
    sht, channel, pick = zip(
        *[(sp, spreadch2ch(sp, ch), t) for sp in sorted(picks)
          for ch, t in sorted(picks[sp].items())
          if spreadch2ch(sp, ch) not in stuff.missing_channels
          ])
    recs = np.genfromtxt(sc['info'] + 'rec_xz.dat', dtype=None,
                         encoding=None)
    shots = np.genfromtxt(sc['info'] + 'sou_xz.dat', dtype=None,
                          encoding=None)
    shot_pos = {sp: (sx, sz) for sp, sx, sz in shots}
    rec_pos = {sp: (sx, sz) for sp, sx, sz in recs}
    picksnew = list(zip(
        *[(sp,
           rec_pos[ch2rec[spreadch2ch(sp, ch)]][0] - shot_pos[sp][0],
           t #+ sc.get('delay', 0.0)
           )
          for sp in sorted(picks)
          for ch, t in sorted(picks[sp].items())
          if spreadch2ch(sp, ch) not in stuff.missing_channels
          ]))
    return picksnew


def load_stuff(verbose=False):
    with open('seisconf.json') as f:
        sc = json.load(f)
    shotpoints = _try_read_file(
        su.read_shotpoints, sc['info'] + 'shotpoints.txt', verbose=verbose)
    spreads = _try_read_file(
        su.read_spreads, sc['info'] + 'spreads.txt', verbose=verbose)
    stuff = calc_backshots(shotpoints, spreads, verbose=verbose)
    picks = load_picks('picks1.txt')
    return convert_picks(picks, stuff, sc)


def linear_fit_through(point, x, y):
    x0, y0 = point
    A = np.transpose([x-x0])
    result = np.linalg.lstsq(A, y-y0, rcond=None)
    m = result[0][0]
    c = y0 - m * x0
    return (m, c)


def fit_legs(x2, x, t):
    mcs = []
    x2 = np.hstack((0, x2, np.max(x)+1))
    point = (0, 0)
    for i in range(len(x2)-1):
        ind = np.logical_and(x >= x2[i], x < x2[i+1])
        mc = linear_fit_through(point, x[ind], t[ind])
        point = (x2[i+1], x2[i+1] * mc[0] + mc[1])
        mcs.append(mc)
    return mcs, x2


def linear_fit(x, y):
    A = np.transpose(np.vstack([x, np.ones(len(x))]))
    result = np.linalg.lstsq(A, y, rcond=None)
    return result[0]


def fit_points(x2, t2, x, t):
    mcs = []
    for i in range(len(x2)-1):
        mc = linear_fit(x2[i:i+2], t2[i:i+2])
        mcs.append(mc)
    return mcs


def error_func(x2, x, t):
    if any(x2[i] >= x2[i+1] for i in range(len(x2) - 1)):
        return 100000
    mcs, x2 = fit_legs(x2, x, t)
    rms = calc_rms(mcs, x2, x, t)
    return rms


def fit_ndmodel(n, x, t):
    from scipy.optimize import brute
    assert n >= 2
    ranges = [(np.min(x)+1, np.max(x)-1) for _ in range(n-1)]
    x2 = brute(error_func, ranges, args=(x, t), workers=-1)
    return fit_legs(x2, x, t)


def tt2model(mcs):
    v = [1/m for m, _ in mcs]
    tau = [c for _, c in mcs]
    h = []
    N = len(mcs)

    def term(v1, v2):
        return (1/v1**2 - 1/v2**2)**0.5

    for n in range(N - 1):
        tau0 = sum(2*h[j]*term(v[j], v[n]) for j in range(n))
        hn = (tau[n+1]-tau0) / 2 / term(v[n], v[n+1])
        h.append(hn)
    for i in range(N-1):
        print(f'  <{h[i]:.1f}m: v={v[i]/1000:.2f}km/s')
    h2 = h[-1] if len(h) > 0 else 0
    print(f'  >{h2:.1f}m: v={v[-1]/1000:.2f}km/s\n')
    return h, v


def calc_rms(mcs, x2, x, t):
    residuals = []
    for i in range(len(x2)-1):
        ind = np.logical_and(x >= x2[i], x < x2[i+1])
        m, c = mcs[i]
        residuals.append(t[ind] - m * x[ind] - c)
    # print(residuals)
    # print(np.shape(residuals))
    residuals = np.hstack(residuals)

    rms = (np.mean(residuals ** 2)) ** 0.5
    return rms


def fit_model(x, t, n=None, ttpicks=None, mode=None):
    x = np.abs(x)
    if ttpicks and mode == 't':
        x2, t2 = zip(*ttpicks)
        mcs = fit_points(x2, t2, x, t)
    elif ttpicks and mode == 'r':
        x2, _ = zip(*ttpicks)
        mcs, x2 = fit_legs(x2, x, t)
    else:
        assert n is not None
        if n == 1:
            mcs, x2 = fit_legs([], x, t)
        else:
            mcs, x2 = fit_ndmodel(n, x, t)
    # plot travel time curve
    for i in range(len(x2)-1):
        plt.plot(x2[i:i+2], mcs[i][0]*np.array(x2[i:i+2])+mcs[i][1], 'C0')
    print(f'Fitted {len(mcs)}d model')
    rms = calc_rms(mcs, x2, x, t) * 1000
    print(f'RMS: {rms:.2f}ms')
    tt2model(mcs)


def plot_and_show(n, x, t, twosides=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('distance (m)')
    ax.set_ylabel('travel time (s)')
    n = np.array(n)
    x = np.array(x)
    t = np.array(t)
    if not twosides:
        x = np.abs(x)
    markers = ax.scatter(x, t, color='k', marker='x', picker=True, alpha=0.5)
    fit_data = [x, t]
    mode = None
    ttpicks = []

    def onpick(event):
        if mode is not None:
            return
        ind = event.ind
        print(f'Picked {len(ind)} points\n'
              f'select shots {np.unique(n[ind])}\n'
              f'position {np.median(x[ind]):.1f}m\n'
              f'ttime    {np.median(t[ind]):.4f}s\n')
        # for n0 in np.unique(n[ind]):
        #     colors[n == n0] = 'm'
        ind2 = np.logical_or.reduce([n == n0 for n0 in n[ind]])
        colors = np.array(['k'] * len(n))
        colors[ind2] = 'm'
        fit_data[:] = [x[ind2], t[ind2]]
        markers.set_color(colors)
        fig.canvas.draw()

    def keypress(event):
        nonlocal mode, markers
        if event.inaxes != ax:
            return
        print('You pressed ' + event.key)
        if event.key == 'h':
            print(_key_doc)
        elif event.key == 'c':
            markers.set_color('k')
            fit_data[:] = [x, t]
        elif event.key == 't':
            if mode == 't':
                print('Exit pick mode for travel time curves')
                fit_model(*fit_data, ttpicks=ttpicks, mode=mode)
                mode = None
            else:
                print('Enter pick mode for travel time curves')
                mode = 't'
                ttpicks[:] = [(0, 0)]
                ax.scatter(0, 0, marker='o', color='C0')
        elif event.key == 'r':
            if mode == 'r':
                print('Exit pick mode for travel time legs')
                fit_model(*fit_data, ttpicks=ttpicks, mode=mode)
                mode = None
            else:
                print('Enter pick mode for travel time legs')
                mode = 'r'
                ttpicks[:] = []
        elif event.key in '123456':
            if mode is not None:
                print(f'Please exit {mode}-mode first.')
                return
            fit_model(*fit_data, n=int(event.key))
        elif event.key == 'x':
            print('Clear axis from tt picks, deselect shots')
            mode = None
            ax.clear()
            markers = ax.scatter(
                x, t, color='k', marker='x', picker=True, alpha=0.5)
            fit_data[:] = [x, t]
        elif event.key == 'I':
            # Danger zone - program must be killed afterwards
            from IPython import embed
            embed()
            sys.exit()
        fig.canvas.draw()

    def buttonpress(event):
        if event.inaxes != ax or mode is None:
            return
        xy = (event.xdata, event.ydata)
        ttpicks.append(xy)
        if mode == 't':
            ax.scatter(*xy, marker='o', color='C0')
        elif mode == 'r':
            ax.axvline(xy[0], color='C0')
        fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)
    fig.canvas.mpl_connect('key_press_event', keypress)
    fig.canvas.mpl_connect('button_press_event', buttonpress)
    plt.show()


def main():
    args = sys.argv
    if '-h' in args or '--help' in args:
        print(__doc__)
        sys.exit()
    picks = load_stuff('-v' in args)
    plot_and_show(*picks, twosides='-n' in args)


if __name__ == '__main__':
    main()
