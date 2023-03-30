# (C) 2018-2020, Tom Eulenfeld, MIT license
"""
Possible invocations:
    seispick  (needs seisconf.json)
    seispick shot_number (needs seisconf.json)
    seispick RAW/{}.dat
    seispick RAW/{}.dat shot_number
    seispick RAW/1001.dat  (single file)

    -v    verbose mode
    -h    show help

version 2023.03 add export of picks to pygimli format
version 2022.09 fix bug: order config option in spreads.txt usable again
version 2021.06 intercept of line fit now given at shot point
                add least-squares line fit
version 2020.10 options to reload config, plot individual traces
                robust parsing of seisconf file and command line options
                display all possible backshot picks
version 2020.09 create textboxes in Tkinter for faster plotting
                supress ObsPy warning when reading data
version 2019.09 add delay config option for seisconf.json
version 2019.05 first version

Upon start set terminal to "always on top" and enlarge matplotlib figure.
"""

_key_doc = """
k: print help for keys
+/-: next/previous shot
t: set picks continuously
d: delete picks continuously
r: set/delete a single pick
s/w: write picks to file
F1,F2,..: switch pick type
o: goto shot number of backshot
h: hide trace
i: show individual traces of stack
left/right: go left/right
up/down: go up/down

g/b: increase/decrease amplitude (gain)
c: clip on/off
n: normalize on/off
p: switch polarity
f: apply/remove filter
a/A: apply autopicker for traces on the right/left
m: fit line to two points
M: least-squares line fit

x: replot
l: reload configuration
I: start IPython session

Set shot number in textbox 1
Set filter in textbox 2: HP/LP/BP freq1 (freq2) (order) (zerophase 0/1)
Set autopicker in textbox 3: threshold min_time_template max_time_template min_time_data max_time_data
Click wiggle for info on trace
Click backshot pick for info
"""

from collections import defaultdict
from copy import copy
import json
import os.path
import sys
from types import SimpleNamespace
import warnings
from warnings import warn

import matplotlib as mpl
MS = mpl.rcParams['lines.markersize']
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import numpy as np
from obspy import read
from obspy.signal.cross_correlation import correlate_template
from scipy.linalg import lstsq
from scipy.stats import linregress
import tkinter as tk


_DEFAULT_PICK_COLORS = ['C0', 'C1', 'C2', 'C3', 'C4']


def load_picks(fname):
    picks = defaultdict(dict)
    if os.path.exists(fname):
        with open(fname) as f:
            data = f.read()
        for line in data.splitlines():
            sp, i, t = line.split()
            picks[int(sp)][int(i)] = float(t)
    return picks


def write_picks(picks, fname1, fname2=None, fname3=None,
                stuff=None, conf=None):
    with open(fname1, 'w') as f:
        for sp in sorted(picks):
            for i, t in sorted(picks[sp].items()):
                f.write(f'{sp}   {i}  {t:.7f}\n')
    if fname2 is None or fname3 is None:
        return
    ch2rec = stuff.ch2rec


    def spreadch2ch(sp, ch):
        return ch - 1 + stuff.spread_from_sp[sp]['channels'][0]
    sht, channel, pick = zip(
        *[(sp, spreadch2ch(sp, ch), t) for sp in sorted(picks)
          for ch, t in sorted(picks[sp].items())
          if spreadch2ch(sp, ch) not in stuff.missing_channels
          ])
    recs = np.genfromtxt(conf.info + 'rec_xz.dat', dtype=None,
                         encoding=None)
    shots = np.genfromtxt(conf.info + 'sou_xz.dat', dtype=None,
                          encoding=None)
    shot_pos = {sp: (sx, sz) for sp, sx, sz in shots}
    rec_pos = {rec: (sx, sz) for rec, sx, sz in recs}
    profil = [rec_pos[ch2rec[ch]][0] for ch in channel]
    profil_shots = [shot_pos[sp][0] for sp in sht]
    offset = [p - s for p, s in zip(profil, profil_shots)]
    depth = -np.array([rec_pos[ch2rec[ch]][1] for ch in channel])
    pick1 = np.maximum(0, pick)
    if fname2:
        error = 1 + np.abs(offset) * conf.error
        out = [sht, profil, conf.ref_depth + depth, 1000 * pick1, error]
        fmt = '%d %.3f %.3f %.3f %.2f'
        np.savetxt(fname2, list(zip(*out)), fmt=fmt)
    if fname3:
        depth_shots = [-shot_pos[sp][1] for sp in sht]
        out = [pick1, profil_shots, depth_shots, profil, depth]
        fmt = '%.8f 0 %.3f %.3f %.3f %.3f'
        np.savetxt(fname3, list(zip(*out)), fmt=fmt)


def calc_backshots(shotpoints, spreads, verbose=True):
    all_channels = sorted(set().union(
        *[set(spread['channels']) for spread in spreads]))
    missing_channels = [ch for ch in range(
        all_channels[0], all_channels[-1]+1) if ch not in all_channels]
    ch2rec = {ch: ch - np.sum(np.array(missing_channels) < ch)
              for ch in all_channels}
    rec2ch = {ch - np.sum(np.array(missing_channels) < ch):
              ch for ch in all_channels}
    spread_from_sp = {sp: spread for spread in spreads
                      for sp in spread['shotpoints']}
    _spread_id = lambda sp: spread_from_sp[sp]['shotpoints'][0]
    shotpoints_rec = {(rec, _spread_id(sp)): sp
                      for sp, rec in shotpoints.items()}
    shotpoints_rec2 = {}
    for sp, rec in shotpoints.items():
        shotpoints_rec2[rec] = shotpoints_rec2.get(rec, ()) + (sp,)
    backshots = {}
    backshots2 = {}
    for sp, recshot in shotpoints.items():
        for ch in spread_from_sp[sp]['channels']:
            spreadch = ch - spread_from_sp[sp]['channels'][0] + 1
            rec = ch2rec[ch]
            if recshot != rec:
                sp2 = shotpoints_rec.get((rec, _spread_id(sp)))
                if sp2 is not None:
                    if rec2ch[recshot] in spread_from_sp[sp2]['channels']:
                        spreadch2 = rec2ch[recshot] - \
                                spread_from_sp[sp2]['channels'][0] + 1
                        backshots[(sp, spreadch)] = (sp2, spreadch2)
                sps2 = shotpoints_rec2.get(rec, ())
                for sp2 in sps2:
                    if rec2ch[recshot] in spread_from_sp[sp2]['channels']:
                        spreadch2 = rec2ch[recshot] - \
                                spread_from_sp[sp2]['channels'][0] + 1
                        ind = (sp, spreadch)
                        ind2 = (sp2, spreadch2)
                        backshots2[ind] = backshots2.get(ind, ()) + (ind2,)
        for sp2, recshot2 in shotpoints.items():
            if sp != sp2 and recshot == recshot2:  # second shot point at the same position
                for ch in spread_from_sp[sp]['channels']:
                    if ch in spread_from_sp[sp2]['channels']:
                        spreadch = ch - spread_from_sp[sp]['channels'][0] + 1
                        spreadch2 = ch - spread_from_sp[sp2]['channels'][0] + 1
                        ind = (sp, spreadch)
                        ind2 = (sp2, spreadch2)
                        backshots2[ind] = backshots2.get(ind, ()) + (ind2,)
    if verbose:
        print('\nbackshots')
        print(backshots)
        print('\nall backshots')
        print(backshots2)
    return SimpleNamespace(backshots=backshots, backshots2=backshots2,
                           spread_from_sp=spread_from_sp,
                           missing_channels=missing_channels, ch2rec=ch2rec)


def interp(x1, y1, x2, y2, n):
    if x1 is None or x2 is None or y1 is None or y2 is None:
        return (), ()
    if x1 > x2:
        x1, x2, y1, y2 = x2, x1, y2, y1
    xp = np.arange(int(np.ceil(x1)), int(np.floor(x2))+1)
    yp = np.interp(xp, (x1, x2), (y1, y2))
    ind = np.logical_and(xp >= 1, xp <= n)
    return xp[ind], yp[ind]


def fit_line(x1, y1, x2, y2, x0, fname=None):
    if x1 > x2:
        x1, x2, y1, y2 = x2, x1, y2, y1
    if fname is not None:
        recs = np.genfromtxt(fname, dtype=None, encoding=None)
        rec_pos = {sp: (sx, sz) for sp, sx, sz in recs}
        x0 = rec_pos[x0][0]
        x1 = rec_pos[int(round(x1))][0] - x0
        x2 = rec_pos[int(round(x2))][0] - x0
    a = (y1 - y2) / (x1 - x2)
    b = y1 - a * x1
    print('line: y = a * x + b')
    print(f'a = {a:.4g}')
    print(f'1 / a = {1/a:.1f}')
    print(f'b = {b:.3g}')


def fit_line_ls(pr, x0, fname=None):
    x, t = zip(*pr)
    if fname is not None:
        recs = np.genfromtxt(fname, dtype=None, encoding=None)
        rec_pos = {sp: (sx, sz) for sp, sx, sz in recs}
        x0 = rec_pos[x0][0]
        xf = [rec_pos[xx][0] - x0 for xx in x]
    else:
        xf = x
    a, b, *_ = linregress(xf, t)
    print('least squares fit: y = a * x + b')
    print(f'a = {a:.4g}')
    print(f'1 / a = {1/a:.1f}')
    print(f'b = {b:.3g}')
    # just for plotting a line
    y0 = a * xf[0] + b
    y1 = a * xf[-1] + b
    a2 = (y1 - y0) / (x[-1] - x[0])
    b2 = y0 - a2 * x[0]
    return lambda x: a2 * x + b2


def _try_read_file(func, fname, verbose=True):
    try:
        return func(fname, verbose=verbose)
    except FileNotFoundError as ex:
        print("Don't worry, but", ex)
        return None


def _oo(b):
    return 'on' if b else 'off'


def _get_points(picks, nshot):
    return (zip(*picks[nshot].items()) if len(picks[nshot]) > 0 else (0, 0))


class MPLSeisPicker(object):
    def __init__(self, *args):
        if '-h' in args or '--help' in args:
            print(__doc__)
            sys.exit()
        expr = None
        nshot = None
        verbose = False
        for arg in args:
            try:
                nshot = int(arg)
            except:
                if arg == '-v':
                    verbose = True
                else:
                    expr = arg
        self.nshot = 1 if nshot is None else nshot
        self.opt = SimpleNamespace(
            norm=1, normalize=True, cut=False,
            apply_filter=False, filter=None,
            polarity=1,
            downsample=1,
            individual=False,
            hide_trace=defaultdict(bool)
        )
        self.state = SimpleNamespace(
            last_key=None, last_xy=None,
            pickmarkers=[None, None], wiggles=None)
        self.stream = None
        self.stream_all = None
        self.ostream = None
        self.ostream_all = None
        self.load_conf(expr=expr, verbose=verbose)
        self.all_picks = [load_picks('picks{}.txt'.format(i+1))
                          for i in range(self.conf.num_pt)]
        self.picktype = 0
        self.picks = self.all_picks[self.picktype]
        filter0 = 'LP 100 2 1'
        picker0 = '0.8 -1e-3 1e-3 -2e-3 3e-3'

        self.root = root = tk.Tk()
        self.root.wm_title('seispick')
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = canvas = FigureCanvasTkAgg(self.fig, master=root)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
        toolbar.update()
        l1 = tk.Label(root, text='shot')
        self.e1 = e1 = tk.Entry(root)
        e1.insert(0, str(self.nshot))
        l2 = tk.Label(root, text='filter')
        self.e2 = e2 = tk.Entry(root)
        e2.insert(0, filter0)
        l3 = tk.Label(root, text='autopicker')
        self.e3 = e3 = tk.Entry(root)
        e3.insert(0, picker0)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        l1.pack(side=tk.LEFT)
        e1.pack(side=tk.LEFT)
        l2.pack(side=tk.LEFT)
        e2.pack(side=tk.LEFT)
        l3.pack(side=tk.LEFT)
        e3.pack(side=tk.LEFT)
        toolbar.pack(side=tk.LEFT)
        self.plot('first_call')

        e1.bind('<Return>', self.new_shot)
        e1.bind('<FocusOut>', self.new_shot)
        e2.bind('<Return>', self.set_filter)
        e2.bind('<FocusOut>', self.set_filter)
        e3.bind('<Return>', self.set_picker)
        e3.bind('<FocusOut>', self.set_picker)
        canvas.mpl_connect('key_press_event', self.keypress)
        canvas.mpl_connect('button_press_event',
                           lambda event:canvas._tkcanvas.focus_set())
        canvas.mpl_connect('pick_event', self.on_pick)
        self.set_filter()
        self.set_picker()

        print(__doc__)
        print(_key_doc)

    def load_conf(self, expr=None, verbose=False):
        import seisutil as su
        try:
            with open('seisconf.json') as f:
                sc = json.load(f)
        except FileNotFoundError as ex:
            print("Dont't worry, but", ex)
            sc = dict(expr=expr)
            self.synthetic = None
            kw = dict(nshot_max = 9999)
            self.stuff = SimpleNamespace(spread_from_sp=None)
        else:
            if expr is not None:
                raise ValueError('expr given as arg and in seisconf file')
            self.synthetic = sc.get('synthetic', None)
            if self.synthetic is not None:
                self.synthetic = load_picks(sc['synthetic'])
            shotpoints = _try_read_file(su.read_shotpoints, sc['info'] + 'shotpoints.txt', verbose=verbose)
            spreads = _try_read_file(su.read_spreads, sc['info'] + 'spreads.txt', verbose=verbose)
            filenumbers = _try_read_file(su.read_filenumbers, sc['info'] + 'filenumbers.txt', verbose=verbose)
            # calculation for backshots and for writing fb_all.dat
            try:
                self.stuff = calc_backshots(shotpoints, spreads, verbose=verbose)
            except Exception:
                import traceback
                traceback.print_exc()
                print("Don't worry, not showing picks of backshots")
                if spreads is not None:
                    spread_from_sp = {sp: spread for spread in spreads
                          for sp in spread['shotpoints']}
                else:
                    spread_from_sp = None
                self.stuff = SimpleNamespace(spread_from_sp=spread_from_sp)
            nshot_max = 9999 if shotpoints is None else max(shotpoints)
            kw = dict(filenumbers=filenumbers, shotpoints=shotpoints,
                      spreads=spreads, nshot_max=nshot_max)
        self.conf = SimpleNamespace(
            expr=sc['expr'], stack=sc.get('stack', False), info=sc.get('info'),
            error=sc.get('error', 0.008), delay=sc.get('delay', 0.0),
            ref_depth=sc.get('ref_depth', 0),
            num_pt=sc.get('number_pick_types', 1),
            pick_colors=sc.get('pick_colors', _DEFAULT_PICK_COLORS),
            **kw)

    def on_pick(self, event):
        line = event.artist
        for i, l in self.state.wiggles.items():
            if l == line:
                stream = self.stream_all if self.opt.individual else self.stream
                print('You clicked on trace')
                tr = stream[i]
                time = tr.stats.endtime - tr.stats.starttime
                sr = tr.stats.sampling_rate
                print(tr.id[:-1], f'| {time:.2f}s  {sr:.1f}Hz')
                break
        else:
            print('You clicked on backshot')
            x = int(round(event.artist.get_xdata()[event.ind][0]))
            bss = self.stuff.backshots2.get((self.nshot, x), ())
            for ind in range(self.conf.num_pt):
                times = []
                for bs in bss:
                    picks = self.all_picks[ind]
                    t = picks[bs[0]].get(bs[1])
                    if t is not None:
                        times.append((t, bs[0], bs[1]))
                for t, sh, rec in sorted(times):
                    print(f'Pick{ind} for shot {sh} receiver {rec} at {t:.4f}s')

    def keypress(self, event):
        if event.inaxes != self.ax:
            return
        print('You pressed ' + event.key)
        if event.key == 'k':
            print(_key_doc)
        elif event.key in 'gb':
            factor = 1.2 if self.opt.normalize else 2
            if event.key == 'g':
                self.opt.norm *= factor
            else:
                self.opt.norm /= factor
            print(f'Normalize by {self.opt.norm:.2f}')
            self.plot('update_wiggles')
        elif event.key == 'c':
            self.opt.cut = not self.opt.cut
            print(f'Clip is {_oo(self.opt.cut)}')
            self.plot('update_wiggles')
        elif event.key == 'p':
            self.opt.polarity = (self.opt.polarity + 2) % 3 - 1
            print(f'Polarity set to {self.opt.polarity}')
            self.plot('update_wiggles')
        elif event.key == 'h' and event.xdata is not None:
            i = int(round(event.xdata))
            self.opt.hide_trace[i] = not self.opt.hide_trace[i]
            self.plot('update_wiggles')
        elif event.key == 'n':
            self.opt.normalize = not self.opt.normalize
            self.opt.norm = 1
            print(f'Normalization turned {_oo(self.opt.normalize)}')
            print(f'Normalize by {self.opt.norm:.2f}')
            self.plot('update_stream')
        elif event.key == 'f':
            self.opt.apply_filter = not self.opt.apply_filter
            print(f'Filter turned {_oo(self.opt.apply_filter)}')
            self.plot('update_stream')
        elif event.key == 'x':
            self.plot('new_shot')
        elif event.key == 'i':
            if self.conf.stack:
                self.opt.individual = not self.opt.individual
                self.plot('new_shot')
            else:
                print('Cannot plot individual shots for stack=False')
        elif event.key == 'l':
            self.load_conf()
            self.plot('new_shot')
        elif event.key in ('-', '+'):
            if event.key == '+':
                self.nshot += 1
            else:
                self.nshot -= 1
            if self.nshot > self.conf.nshot_max:
                self.nshot = 1
            if self.nshot < 1:
                self.nshot = self.conf.nshot_max
            self.new_shot()
        elif event.key == self.state.last_key and event.key in 'td':
            xp, yp = interp(event.xdata, event.ydata, *self.state.last_xy,
                            len(self.stream))
            for x, y in zip(xp, yp):
                if event.key == 't':
                    self.picks[self.nshot][x] = y
                else:
                    self.picks[self.nshot].pop(x, None)
            self.plot('update_picks')
        elif event.key == 'r':
            x = int(round(event.xdata))
            try:
                self.picks[self.nshot].pop(x)
            except KeyError:
                self.picks[self.nshot][x] = event.ydata
            self.plot('update_picks')
        elif event.key == 'o':
            x = int(round(event.xdata))
            y = event.ydata
            try:
                bss = self.stuff.backshots2[(self.nshot, x)]
            except KeyError:
                print(f'No backshot for shot {self.nshot} and receiver {x}')
            else:
                nshot, xn = self.stuff.backshots.get((self.nshot, x),
                                                     (None, None))
                for bs in bss:
                    if (nshot is None or
                            abs(self.picks[bs[0]].get(bs[1], -1000) - y) <
                            abs(self.picks[nshot].get(xn, -1001) - y)):
                        nshot, xn = bs
                self.nshot = nshot
                print(f'Go to shot {nshot} receiver {xn}')
                # center receiver of backshot to cursor
                x1, x2 = self.ax.get_xlim()
                self.ax.set_xlim(x1 + xn - x, x2 + xn - x)
                self.new_shot()
        elif event.key in 'aA':
            x = int(round(event.xdata))
            d = -1 if event.key == 'A' else 1
            self.apply_autopicker(x, d)
            self.plot('update_picks')
        elif event.key in ('left', 'right'):
            x1, x2 = self.ax.get_xlim()
            w = int(0.4 * (x2 - x1))
            d = -1 if event.key == 'left' else 1
            self.ax.set_xlim(x1 + d * w, x2 + d * w)
            self.fig.canvas.draw()
        elif event.key in ('up', 'down'):
            y1, y2 = self.ax.get_ylim()
            h = 0.4 * (y2 - y1)
            d = -1 if event.key == 'down' else 1
            self.ax.set_ylim(y1 + d * h, y2 + d * h)
            self.fig.canvas.draw()
        elif event.key in 'sw':
            print('Save picks')
            for i in range(self.conf.num_pt):
                try:
                    write_picks(self.all_picks[i], f'picks{i+1}.txt',
                                    f'fb_all{i+1}.dat',
                                    f'fb_all_gimli{i+1}.dat',
                                    stuff=self.stuff, conf=self.conf)
                except Exception:
                    import traceback
                    traceback.print_exc()
                    print("Don't worry, picks are anyway written to picks?.txt file")
        elif event.key == self.state.last_key == 'm':
            x1, y1, x2, y2 = event.xdata, event.ydata, *self.state.last_xy
            x0 = self.conf.shotpoints[self.nshot]
            if not (x1 is None or x2 is None or y1 is None or y2 is None):
                fit_line(x1, y1, x2, y2, x0, fname=self.conf.info + 'rec_xz.dat')
                self.ax.plot((x1, x2), (y1, y2), color='C1')
                self.fig.canvas.draw()
            event.key = None
        elif event.key == self.state.last_key == 'M':
            x1, x2 = event.xdata, self.state.last_xy[0]
            x0 = self.conf.shotpoints[self.nshot]
            if not (x1 is None or x2 is None):
                if x1 > x2:
                    x1, x2 = x2, x1
                x1 = int(round(x1))
                x2 = int(round(x2))
                pickrange = [(x, p) for x, p in self.picks[self.nshot].items()
                             if x in range(x1, x2+1)]
                if len(pickrange) < 3:
                    print('Not enough picks for least squares fit.')
                else:
                    func = fit_line_ls(pickrange, x0, fname=self.conf.info + 'rec_xz.dat')
                    xr = np.array([x1, x2])
                    self.ax.plot(xr, func(xr), color='C1')
                    self.fig.canvas.draw()
            event.key = None
        elif len(event.key) == 2 and event.key[0] == 'f':
            pt = int(event.key[1])
            if pt <= self.conf.num_pt:
                self.picktype = pt - 1
                self.picks = self.all_picks[self.picktype]
                print(f'Switch to pick type {pt}')
        elif event.key == 'I':
            # Danger zone - program must be killed afterwards
            from IPython import embed
            embed()
        self.state.last_key = event.key
        self.state.last_xy = (event.xdata, event.ydata)

    def set_filter(self, event=None):
        """Called if filter box changes"""
        filt = [_.lower() for _ in self.e2.get().split()]
        map_ = {'bp': 'bandpass', 'lp': 'lowpass', 'hp': 'highpass'}
        try:
            ftype = map_[filt[0]]
            if ftype == 'bandpass':
                f1, f2, *filt = filt[1:]
                kw = dict(freqmin=float(f1), freqmax=float(f2))
            else:
                f1, *filt = filt[1:]
                kw = dict(freq=float(f1))
            if len(filt) > 0:
                kw['corners'] = int(filt[0])
            if len(filt) > 1:
                kw['zerophase'] = bool(int(filt[1]))
        except Exception:
            print('Error in filter textbox')
            self.e2.focus_set()
        else:
            self.opt.filter = dict(type=ftype, **kw)
            print(f'You chose filter {self.opt.filter}')
            self.canvas.get_tk_widget().focus_set()

    def set_picker(self, event=None):
        """Called if autopicker box changes"""
        try:
            thres, t1, t2, t3, t4 = [float(_) for _ in self.e3.get().split()]
        except Exception:
            print('Error in picker textbox')
            self.e3.focus_set()
        else:
            self.opt.autopick = SimpleNamespace(thres=thres, t1=t1, t2=t2,
                                                t3=t3, t4=t4)
            print(f'You chose autopicker {self.opt.autopick}')
            self.canvas.get_tk_widget().focus_set()

    def new_shot(self, event=None):
        """Called if shot box changes or if +/- pressed"""
        if event is None:
            self.e1.delete(0, tk.END)
            self.e1.insert(0, str(self.nshot))
        else:
            try:
                nshot = int(self.e1.get())
                if nshot > self.conf.nshot_max or nshot < 1:
                    raise Exception
            except Exception:
                print('Error in nshot textbox')
                self.e1.focus_set()
            else:
                self.nshot = nshot
                self.canvas.get_tk_widget().focus_set()
        self.opt.hide_trace = defaultdict(bool)
        self.plot('new_shot')

    def plot(self, mode):
        modes = ['first_call', 'new_shot', 'update_stream', 'update_wiggles',
                 'update_picks']
        assert mode in modes
        mind = modes.index(mode)
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        if mind < 2:
            self.fig.suptitle('shot %s' % self.nshot)
            try:
                self.read_stream()
            except FileNotFoundError as ex:
                print(ex)
                self.ax.plot(self.ax.get_xlim(), self.ax.get_ylim(), 'r', lw=5)
                self.fig.canvas.draw()
                return
            self.ostream = self.stream
            self.ostream_all = self.stream_all
        if mind < 3:
            stream = self.ostream.copy()
            if self.opt.apply_filter:
                stream.filter(**self.opt.filter)
            if self.opt.normalize:
                stream.normalize()
            else:
                stream.normalize(global_max=True)
            self.stream = stream
            if self.conf.stack:
                stream_all = self.ostream_all.copy()
                if self.opt.apply_filter:
                    stream_all.filter(**self.opt.filter)
                if self.opt.normalize:
                    stream_all.normalize()
                else:
                    stream_all.normalize(global_max=True)
                self.stream_all = stream_all
        if self.opt.individual:
            stream = self.stream_all
        else:
            stream = self.stream
        if mind < 2:
            self.ax.clear()
            self.state.wiggles = {}
            self.ax.axhline(color='0.8')
        elif mind < 4:
            self.ax.collections.clear()
        # plot wiggles
        if mind < 4:
            for i, tr in enumerate(stream):
                gf = i % len(self.stream) + 1 #tr.stats.geofon
                ds = self.opt.downsample
                out = stream[i].data[::ds] * self.opt.norm
                if self.opt.cut:
                    out = np.clip(out, -0.45, 0.45)
                times = tr.times('relative')[::ds] + self.conf.delay
                if mind < 2:
                    if self.opt.individual:
                        lw = 1
                        color = 'gray'
                    else:
                        lw = 0.5
                        color = 'k'
                    self.state.wiggles[i], = self.ax.plot(gf + out, times,
                        color=color, lw=lw,
                        picker=True, pickradius=0.5)
                else:
                    self.state.wiggles[i].set_data(gf + out, times)
                self.state.wiggles[i].set_visible(not self.opt.hide_trace[i+1])
                if (self.opt.polarity != 0 and not self.opt.hide_trace[i+1]
                        and not self.opt.individual):
                    self.ax.fill_betweenx(
                        times, gf, gf+out, where=out*self.opt.polarity > 0,
                        facecolor='k')
        # plot picks and picks of backshots
        if mind < 2:
            if self.synthetic:
                spick_points = list(
                    zip(*self.synthetic[self.nshot].items()))
                if len(spick_points) > 0:
                    self.ax.plot(*spick_points, 'o',
                                 color=self.conf.pick_colors[-1])
            for ind in range(self.conf.num_pt):
                picks = self.all_picks[ind]
                pick_points = _get_points(picks, self.nshot)
                self.state.pickmarkers[ind], = self.ax.plot(
                    *pick_points, 'o', color=self.conf.pick_colors[2*ind])
                self.state.pickmarkers[ind].set_visible(
                    pick_points != (0, 0))
                if hasattr(self.stuff, 'backshots'):
                    bpicks = []
                    for i in range(len(stream)+1):
                        bs = self.stuff.backshots.get((self.nshot, i))
                        if bs is not None and bs[1] in picks[bs[0]]:
                            bpicks.append((i, picks[bs[0]][bs[1]]))
                    bpicks2 = []
                    for i in range(len(stream)+1):
                        bss = self.stuff.backshots2.get((self.nshot, i), ())
                        for bs in bss:
                            if bs[1] in picks[bs[0]]:
                                bpicks2.append((i, picks[bs[0]][bs[1]]))
                    if len(bpicks2) > 0:
                        self.ax.plot(*zip(*bpicks2), 'o', ms = MS//2,
                                     color=self.conf.pick_colors[2*ind+1],
                                     alpha=0.5, picker=True)
                    if len(bpicks) > 0:
                        self.ax.plot(*zip(*bpicks), 'o', ms = MS//2,
                                     color=self.conf.pick_colors[2*ind+1])
        elif mind == 4:
            pick_points = _get_points(self.picks, self.nshot)
            self.state.pickmarkers[self.picktype].set_data(*pick_points)
            self.state.pickmarkers[self.picktype].set_visible(
                pick_points != (0, 0))
        if mind > 0:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        self.fig.canvas.draw()

    def apply_autopicker(self, x, d):
        pick0 = self.picks[self.nshot][x]
        tr0 = self.stream[x-1]
        autopick = self.opt.autopick
        trel = tr0.stats.starttime - self.conf.delay
        t1 = pick0 + autopick.t1
        t2 = pick0 + autopick.t2
        template = tr0.slice(trel + t1, trel + t2)
        rect_kw = dict(edgecolor=None, alpha=0.5)
        rect = Rectangle((x-0.5, t1), 1, t2 - t1, facecolor='C1', **rect_kw)
        self.ax.add_patch(rect)
        for tr in self.stream[x-1+d::d]:
            x += d
            t1 = pick0 + autopick.t1 + autopick.t3
            t2 = pick0 + autopick.t2 + autopick.t4
            data = tr.slice(trel + t1, trel + t2)
            rect = Rectangle((x-0.5, t1), 1, t2 - t1,
                             facecolor='C0', **rect_kw)
            self.ax.add_patch(rect)
            cc = correlate_template(data, template, demean=False)
            index = np.argmax(cc)
            print(f'receiver {x:3d}: cc {cc[index]:.3f}')
            if cc[index] < autopick.thres:
                break
            #pick0 = pick0 + p.t3 + index * tr0.stats.delta
            # sub-sample adjustment by fitting a*i**2 + b*i + c = cc[i]
            # best index then given by -b / 2 / a
            i = np.arange(index-3, index+4)
            i = i[i >= 0]
            i = i[i < len(cc)]
            M = i[:, np.newaxis]**[2, 1, 0]
            (a, b, _), _, _, _ = lstsq(M, cc[i])
            imax = -b / 2 / a
            if a > 0:
                imax = index
            elif imax < i[0]:
                imax = i[0]
            elif imax > i[-1]:
                imax = i[-1]
            pick = pick0 + autopick.t3 + imax * tr0.stats.delta
            print(f'receiver {x:3d}: set pick at {pick:.3f}s')
            self.picks[self.nshot][x] = pick
            pick0 = pick0 + autopick.t3 + index * tr0.stats.delta

    def read_stream(self):
        conf = self.conf
        if conf.stack and conf.filenumbers is None:
            raise ValueError('Need filenumbers.txt file when stacking')
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            if conf.stack:
                streams = [read(conf.expr.format(num))
                           for num in conf.filenumbers[self.nshot]]
            else:
                streams = [read(conf.expr.format(self.nshot))]
            for j, stream in enumerate(streams):
                for i, tr in enumerate(stream):
                    tr.stats.network = f'SP{self.nshot}'
                    tr.stats.station = f'GF{i+1}'
                    if conf.stack:
                        tr.stats.location = f'F{conf.filenumbers[self.nshot][j]}'
                    tr.stats.geofon = i + 1
                    tr.stats.shot = self.nshot
                    try:
                        assert int(tr.stats.seg2.CHANNEL_NUMBER) == i+1
                    except:
                        warn('SEG2 channel number does not match, check data')
        if self.stuff.spread_from_sp is not None:
            spread = self.stuff.spread_from_sp[self.nshot]
            for stream in streams:
                if 'mute' in spread:
                    for r in spread['mute']:
                        stream[r-1].data[:] = 0
                if 'polarity' in spread:
                    for r in spread['polarity']:
                        stream[r-1].data = -stream[r-1].data
                if 'order' in spread:
                    for r1, r2 in spread['order']:
                        if r1 == 1:
                            stream[:r2] = stream[r2-1::-1]
                        else:
                            stream[r1-1:r2] = stream[r2-1:r1-2:-1]
        self.stream_all = sum(streams[1:], streams[0])
        self.stream = copy(self.stream_all).stack('{network}.{station}') if conf.stack else self.stream_all
        return self.stream


def run_picker():
    p = MPLSeisPicker(*sys.argv[1:])
    tk.mainloop()


if __name__ == '__main__':
    run_picker()
