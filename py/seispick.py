# (C) 2018, Tom Eulenfeld, MIT license
"""
usage picking tool: seispick  (needs seisconf.json)
usage plotting tool: seispick RAW/{}.dat
usage plotting tool: seispick RAW/{}.dat shot_number

version 2020.09.03 add downsample option for faster? plotting
                   supress ObsPy warning when reading data
version 2019.09 add delay config option for seisconf.json
version 2019.05
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
left/right: go left/right
up/down: go up/down

g/b: increase/decrease amplitude (gain)
c: clip on/off
n: normalize on/off
p: switch polarity

f: apply/remove filter
a/A: apply autopicker for traces on the right/left
m: fit line to two points
x: replot

Specify shot number in textbox 1
Specify downsample factor in textbox 2
Specify filter in textbox 3: HP/LP/BP freq1 (freq2) (order) (zerophase 0/1)
Specify autopicker in textbox 4: threshold min_time_template max_time_template min_time_data max_time_data
"""


from collections import defaultdict
import json
import os.path
import sys
from types import SimpleNamespace
import warnings

import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
from matplotlib.patches import Rectangle
import numpy as np
from obspy import read
from obspy.signal.cross_correlation import correlate_template
from scipy.linalg import lstsq

import matplotlib as mpl
MS = mpl.rcParams['lines.markersize']


def load_picks(fname):
    picks = defaultdict(dict)
    if os.path.exists(fname):
        with open(fname) as f:
            data = f.read()
        for line in data.splitlines():
            sp, i, t = line.split()
            picks[int(sp)][int(i)] = float(t)
    return picks


def write_picks(picks, fname1, fname2=None, stuff=None, conf=None):
    with open(fname1, 'w') as f:
        for sp in sorted(picks):
            for i, t in sorted(picks[sp].items()):
                f.write(f'{sp}   {i}  {t:.7f}\n')
    if fname2 is None:
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
    rec_pos = {sp: (sx, sz) for sp, sx, sz in recs}
    profil = [rec_pos[ch2rec[ch]][0] for ch in channel]
    profil_shots = [shot_pos[sp][0] for sp in sht]
    offset = [p - s for p, s in zip(profil, profil_shots)]
    depth = conf.ref_depth - \
        np.array([rec_pos[ch2rec[ch]][1] for ch in channel])
    error = 1 + np.abs(offset) * conf.error
    pick = [1000 * max(0, p) for p in pick]
    out = [sht, profil, depth, pick, error]
    fmt = '%d %.3f %.3f %.3f %.2f'
    np.savetxt(fname2, list(zip(*out)), fmt=fmt)


def calc_backshots(shotpoints, spreads, verbose=True):
    shotpoints_r = {r: sp for sp, r in shotpoints.items()}
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
    backshots = {}
    for sp, rshot in shotpoints.items():
        for ch in spread_from_sp[sp]['channels']:
            spreadch = ch - spread_from_sp[sp]['channels'][0] + 1
            if (shotpoints[sp] != ch2rec[ch] and ch in ch2rec and
                    ch2rec[ch] in shotpoints_r):
                sp2 = shotpoints_r[ch2rec[ch]]
                if rec2ch[rshot] in spread_from_sp[sp2]['channels']:
                    spreadch2 = rec2ch[rshot] - \
                        spread_from_sp[sp2]['channels'][0] + 1
                    backshots[(sp, spreadch)] = (sp2, spreadch2)
    if verbose:
        print('\nbackshots')
        print(backshots)
    return SimpleNamespace(backshots=backshots, spread_from_sp=spread_from_sp,
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


def fit_line(x1, y1, x2, y2, fname=None):
    if x1 > x2:
        x1, x2, y1, y2 = x2, x1, y2, y1
    if fname is not None:
        recs = np.genfromtxt(fname, dtype=None, encoding=None)
        rec_pos = {sp: (sx, sz) for sp, sx, sz in recs}
        x1 = rec_pos[int(round(x1))][0]
        x2 = rec_pos[int(round(x2))][0]
        # xmax = sorted(rec_pos)[-1]
        # xmax = rec_pos[xmax][0]
    a = (y1 - y2) / (x1 - x2)
    b = y1 - a * x1
    print('line: y = a * x + b')
    print(f'a = {a}')
    print(f'1 / a = {1/a}')
    print(f'b = {b}')


def _oo(b):
    return 'on' if b else 'off'


def _get_points(picks, nshot):
    return (zip(*picks[nshot].items()) if len(picks[nshot]) > 0 else (0, 0))


class MPLSeisPicker(object):
    def __init__(self, expr=None, nshot=None):
        if expr in ('-h', '--help'):
            print(__doc__)
            sys.exit()
        self.pickmode = expr is None
        self.nshot = 1 if nshot is None else nshot
        self.opt = SimpleNamespace(
            norm=1, normalize=True, cut=False,
            apply_filter=False, filter=None,
            polarity=1,
            downsample=1,
            hide_trace=defaultdict(bool)
        )
        self.state = SimpleNamespace(
            last_key=None, last_xy=None,
            pickmarkers=[None, None], wiggles=None)
        self.stream = None
        self.ostream = None
        if not self.pickmode:
            self.conf = SimpleNamespace(expr=expr, stack=False, delay=0,
                                        nshot_max=9999)
        else:
            import seisutil as su
            with open('seisconf.json') as f:
                sc = json.load(f)
            self.synthetic = sc.get('synthetic', None)
            if self.synthetic is not None:
                self.synthetic = load_picks(sc['synthetic'])
            shotpoints = su.read_shotpoints(sc['info'] + 'shotpoints.txt',
                                            verbose=True)
            spreads = su.read_spreads(sc['info'] + 'spreads.txt', verbose=True)
            filenumbers = su.read_filenumbers(sc['info'] + 'filenumbers.txt',
                                              verbose=True)
            self.conf = SimpleNamespace(
                expr=sc['expr'], stack=sc['stack'], info=sc['info'],
                error=sc['error'], delay=sc.get('delay', 0.0),
                filenumbers=filenumbers, shotpoints=shotpoints,
                #                    nshot_max=sorted(shotpoints)[-1],
                nshot_max=max(shotpoints),
                spreads=spreads,
                ref_depth=sc['ref_depth'],
                num_pt=sc['number_pick_types'],
                pick_colors=sc['pick_colors'])
            self.all_picks = [load_picks('picks{}.txt'.format(i+1))
                              for i in range(self.conf.num_pt)]
            self.picktype = 0
            self.picks = self.all_picks[self.picktype]
            # calculation for backshots and for writing fb_all.dat
            self.stuff = calc_backshots(shotpoints, spreads)
        self.fig, self.ax = fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.1)
        self.plot('first_call')
        fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
        fig.canvas.mpl_connect('key_press_event', self.keypress)

        axbox0 = fig.add_axes([0.15, 0.02, 0.03, 0.04])
        textbox0 = TextBox(axbox0, 'Shot', initial=str(self.nshot))
        textbox0.on_submit(self.new_shot)

        axbox1 = fig.add_axes([0.3, 0.02, 0.25, 0.04])
        filter0 = 'LP 100 2 1'
        textbox1 = TextBox(axbox1, 'Filter', initial=filter0)
        textbox1.on_submit(self.set_filter)
        if self.pickmode:
            axbox2 = plt.axes([0.65, 0.02, 0.25, 0.04])
            picker0 = '0.8 -1e-3 1e-3 -2e-3 3e-3'
            textbox2 = TextBox(axbox2, 'Autopicker', initial=picker0)
            textbox2.on_submit(self.set_picker)
        else:
            textbox2 = None
        axbox3 = fig.add_axes([0.20, 0.02, 0.03, 0.04])
        textbox3 = TextBox(axbox3, 'DS', initial=1)
        textbox3.on_submit(self.set_downsample)
        self.boxes = SimpleNamespace(box0=textbox0, box1=textbox1,
                                     box2=textbox2, box3=textbox3)
        self.set_filter(filter0)
        if self.pickmode:
            self.set_picker(picker0)
        print(__doc__)
        print(_key_doc)

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
            try:
                self.nshot, xn = self.stuff.backshots[(self.nshot, x)]
            except KeyError:
                print(f'No backshot for shot {self.nshot} and receiver {x}')
            else:
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
            w = int(0.75 * (x2 - x1))
            d = -1 if event.key == 'left' else 1
            self.ax.set_xlim(x1 + d * w, x2 + d * w)
            self.fig.canvas.draw()
        elif event.key in ('up', 'down'):
            y1, y2 = self.ax.get_ylim()
            h = 0.75 * (y2 - y1)
            d = -1 if event.key == 'down' else 1
            self.ax.set_ylim(y1 + d * h, y2 + d * h)
            self.fig.canvas.draw()
        elif event.key in 'sw':
            print('Save picks')
            for i in range(self.conf.num_pt):
                try:
                    write_picks(self.all_picks[i], f'picks{i+1}.txt',
                                f'fb_all{i+1}.dat', stuff=self.stuff, conf=self.conf)
                except Exception as ex:
                    print(ex)
        elif event.key == self.state.last_key == 'm':
            x1, y1, x2, y2 = event.xdata, event.ydata, *self.state.last_xy
            if not (x1 is None or x2 is None or y1 is None or y2 is None):
                fit_line(x1, y1, x2, y2, fname=self.conf.info + 'rec_xz.dat')
                self.ax.plot((x1, x2), (y1, y2), color='C1')
                self.fig.canvas.draw()
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

    def set_filter(self, text):
        """Called if filter box changes"""
        filt = [_.lower() for _ in text.split()]
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
        else:
            self.opt.filter = dict(type=ftype, **kw)
            print(f'You chose filter {self.opt.filter}')
        self.boxes.box2.stop_typing()

    def set_picker(self, text):
        """Called if autopicker box changes"""
        try:
            thres, t1, t2, t3, t4 = [float(_) for _ in text.split()]
        except Exception:
            print('Error in picker textbox')
        else:
            self.opt.autopick = SimpleNamespace(thres=thres, t1=t1, t2=t2,
                                                t3=t3, t4=t4)
            print(f'You chose autopicker {self.opt.autopick}')
        self.boxes.box2.stop_typing()

    def set_downsample(self, text):
        """Called if downsample box changes"""
        try:
            ds = int(text)
        except Exception:
            print('Error in downsample textbox')
        else:
            self.opt.downsample = ds
            print(f'You chose downsample factor {self.opt.downsample}')
        self.boxes.box3.stop_typing()

    def new_shot(self, text=None):
        """Called if shot box changes or if +/- pressed"""
        if text is None:
            self.boxes.box0.set_val(self.nshot)
#            self.boxes.box0.text = str(self.nshot)
#            self.boxes.box0.text_disp.remove()
#            self.boxes.box0.text_disp = self.boxes.box0._make_text_disp(
#                self.boxes.box0.text)
        else:
            try:
                self.nshot = int(text)
            except Exception:
                print('Error in nshot textbox')
        self.boxes.box0.stop_typing()
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
                self.ostream = self.read_stream()
            except FileNotFoundError as ex:
                print(ex)
                self.ax.plot(self.ax.get_xlim(), self.ax.get_ylim(), 'r', lw=5)
                self.fig.canvas.draw()
                return
        if mind < 3:
            stream = self.ostream.copy()
            if self.opt.apply_filter:
                stream.filter(**self.opt.filter)
            if self.opt.normalize:
                stream.normalize()
            else:
                stream.normalize(global_max=True)
            self.stream = stream
        else:
            stream = self.stream
        if mind < 2:
            self.ax.clear()
            self.state.wiggles = {}
        elif mind < 4:
            self.ax.collections.clear()
        if mind < 4:
            for i, tr in enumerate(stream):
                ds = self.opt.downsample
                out = stream[i].data[::ds] * self.opt.norm
                if self.opt.cut:
                    out = np.clip(out, -0.45, 0.45)
                times = tr.times('relative')[::ds] + self.conf.delay
                if mind < 2:
                    self.state.wiggles[i], = self.ax.plot(
                        1 + i + out, times, 'k', lw=0.5)
                else:
                    self.state.wiggles[i].set_data(1 + i + out, times)
                self.state.wiggles[i].set_visible(not self.opt.hide_trace[i+1])
                if self.opt.polarity != 0 and not self.opt.hide_trace[i+1]:
                    self.ax.fill_betweenx(
                        times, 1+i, 1+i+out, where=out*self.opt.polarity > 0,
                        facecolor='k')
        if self.pickmode:
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
                    bpicks = []
                    for i in range(len(stream)):
                        if (self.nshot, i) in self.stuff.backshots:
                            bs = self.stuff.backshots[(self.nshot, i)]
                            if bs[1] in picks[bs[0]]:
                                bpicks.append((i, picks[bs[0]][bs[1]]))
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
        trel = tr0.stats.starttime
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
        if conf.stack:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                streams_in = [read(conf.expr.format(num))
                              for num in conf.filenumbers[self.nshot]]
            for traces in zip(*streams_in):
                set1 = {tr.stats.seg2.CHANNEL_NUMBER for tr in traces}
                set2 = {tr.stats.seg2.RECEIVER_LOCATION for tr in traces}
                assert len(set1) == 1
                assert len(set2) == 1
                traces[0].data = np.mean([tr.data for tr in traces], axis=0)
            stream = streams_in[0]
        else:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                stream = read(conf.expr.format(self.nshot))
        if self.pickmode:
            spread = self.stuff.spread_from_sp[self.nshot]
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
        return stream


def run_picker():
    MPLSeisPicker(*sys.argv[1:])
    plt.show()
