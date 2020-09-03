#!/usr/bin/env python
# (C) 2018, Tom Eulenfeld, MIT license
"""
usage picking tool: seispick  (needs seisconf.json)
usage plotting tool: seispick RAW/{}.dat
usage plotting tool: seispick RAW/{}.dat shot_number

version 2019_02

Keys:
k: print this help
+/-: next/previous shot
up/down: set amplitude
c: clip on/off
n: normalize on/off
p: switch polarity
t: set picks
d: delete picks
h: hide trace
l/r: go left/right
s/w: write picks
f: apply/remove filter
a/y: apply autopicker for traces on the right/left
x: replot
m: fit line to two points
42enter: go to shot 42

Specify filter in textbox: HP/LP/BP freq1 (freq2) (order) (zerophase 0/1)
Specify autopicker in textbox: threshold min_template max_template min_data max_data
"""

from collections import defaultdict
import json
import os.path
import sys
from types import SimpleNamespace

import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
from matplotlib.patches import Rectangle
import numpy as np
from obspy import read
from obspy.signal.cross_correlation import correlate_template
from scipy.linalg import lstsq


def write_picks():
    with open('picks.txt', 'w') as f:
        for sp in sorted(PICKS):
            for i, t in sorted(PICKS[sp].items()):
                f.write(f'{sp}   {i}  {t:.4f}\n')
    def spreadch2ch(sp, ch):
        return ch - 1 + spread_from_sp[sp]['channels'][0]
    sht, channel, pick = zip(*[(sp, spreadch2ch(sp, ch), t)
                               for sp in sorted(PICKS)
                               for ch, t in sorted(PICKS[sp].items())
                               if spreadch2ch(sp, ch) not in missing_channels
                               ])
    recs = np.genfromtxt(sc['info'] + 'rec_xz.dat', dtype=None, encoding=None)
    shots = np.genfromtxt(sc['info'] + 'sou_xz.dat', dtype=None, encoding=None)
    shot_pos = {sp: (sx, sz) for sp, sx, sz in shots}
    rec_pos = {sp: (sx, sz) for sp, sx, sz in recs}
    profil = [rec_pos[ch2rec[ch]][0] for ch in channel]
    profil_shots = [shot_pos[sp][0] for sp in sht]
    offset = [p - s for p, s in zip(profil, profil_shots)]
    depth = conf.ref_depth - np.array([rec_pos[ch2rec[ch]][1] for ch in channel])
    error = 1 + np.abs(offset) * sc['error']
    out = [sht, profil, depth, pick, error]
    fmt = '%d %.3f %.3f %.3f %.2f'
    np.savetxt('fb_all.dat', list(zip(*out)), fmt=fmt)


def load_picks():
    picks = defaultdict(dict)
    if os.path.exists('picks.txt'):
        with open('picks.txt') as f:
            data = f.read()
        for line in data.splitlines():
            sp, i, t = line.split()
            picks[int(sp)][int(i)] = float(t)
    return picks


def autopick(x, d):
    pick0 = PICKS[opt.nshot][x]
    tr0 = state.stream[x-1]
    trel = tr0.stats.starttime
    p = opt.pick
    t1 = pick0 + p.t1
    t2 = pick0 + p.t2
    template = tr0.slice(trel + t1, trel + t2)
    rect_kw = dict(edgecolor=None, alpha=0.5)
    rect = Rectangle((x-0.5, t1), 1, t2 - t1, facecolor='C1', **rect_kw)
    AX.add_patch(rect)
    for tr in state.stream[x-1+d::d]:
        x += d
        t1 = pick0 + p.t1 + p.t3
        t2 = pick0 + p.t2 + p.t4
        data = tr.slice(trel + t1, trel + t2)
        rect = Rectangle((x-0.5, t1), 1, t2 -t1, facecolor='C0', **rect_kw)
        AX.add_patch(rect)
        cc = correlate_template(data, template, demean=False)
        index = np.argmax(cc)
        print(f'cc {cc[index]}')
        if cc[index] < opt.pick.thres:
            break
        #pick0 = pick0 + p.t3 + index * tr0.stats.delta
        # sub-sample adjustment by fitting a*i**2 + b*i + c = cc[i]
        # best index then given by -b / 2 / a
        i = np.arange(index-3, index+4)
        i = i[i>=0]
        i = i[i<len(cc)]
        M = i[:, np.newaxis]**[2, 1, 0]
        (a, b, _), _, _, _ = lstsq(M, cc[i])
        imax = -b / 2 / a
        if a > 0:
            imax = index
        pick = pick0 + p.t3 + imax * tr0.stats.delta
        PICKS[opt.nshot][x] = pick
        pick0 = pick0 + p.t3 + index * tr0.stats.delta


def read_stream():
    if conf.stack:
        streams_in = [read(conf.expr.format(num)) for num in conf.filenumbers[opt.nshot]]
        for traces in zip(*streams_in):
            set1 = {tr.stats.seg2.CHANNEL_NUMBER for tr in traces}
            set2 = {tr.stats.seg2.RECEIVER_LOCATION for tr in traces}
            assert len(set1) == 1
            assert len(set2) == 1
            traces[0].data = np.mean([tr.data for tr in traces], axis=0)
        stream = streams_in[0]
    else:
        stream = read(conf.expr.format(opt.nshot))
    if PICKMODE:
        spread = spread_from_sp[opt.nshot]
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


def plot(mode):
    modes = ['first_call', 'new_shot', 'update_stream', 'update_wiggles',
             'update_picks']
    assert mode in modes
    mind = modes.index(mode)
    xlim = AX.get_xlim()
    ylim = AX.get_ylim()
    if mind < 2:
        FIG.suptitle('shot %s' % opt.nshot)
        state.ostream = read_stream()
    if mind < 3:
        stream = state.ostream.copy()
        if opt.apply_filter:
            stream.filter(**opt.filter)
        if opt.normalize:
            stream.normalize()
        else:
            stream.normalize(global_max=True)
        state.stream = stream
    else:
        stream = state.stream
    if mind < 2:
        AX.clear()
        state.wiggles = {}
    elif mind < 4:
        AX.collections.clear()
    if mind < 4:
        for i, tr in enumerate(stream):
            out = stream[i].data * opt.norm
            if opt.cut:
                out = np.clip(out, -0.45, 0.45)
            times = tr.times('relative')
            if mind < 2:
                state.wiggles[i], = AX.plot(1 + i + out, times, 'k', lw=0.5)
            else:
                state.wiggles[i].set_data(1 + i + out, times)
            state.wiggles[i].set_visible(not opt.hide_trace[i+1])
            if opt.polarity != 0 and not opt.hide_trace[i+1]:
                AX.fill_betweenx(times, 1+i, 1+i+out,
                                 where=out*opt.polarity>0, facecolor='k')
    if PICKMODE:
        if mind in (0, 1, 4):
            picks = zip(*PICKS[opt.nshot].items()) if len(PICKS[opt.nshot]) > 0 else (0, 0)
        if mind < 2:
            state.pickmarkers, = AX.plot(*picks, 'o', color='C0')
            bpicks = []
            for i in range(len(stream)):
                if (opt.nshot, i) in BS:
                    bs = BS[(opt.nshot, i)]
                    if bs[1] in PICKS[bs[0]]:
                        bpicks.append((i, PICKS[bs[0]][bs[1]]))
            if len(bpicks) > 0:
                AX.plot(*zip(*bpicks), 'o', color='C1')
        if mind == 4:
            state.pickmarkers.set_data(*picks)
    if mind > 0:
        AX.set_xlim(xlim)
        AX.set_ylim(ylim)
    FIG.canvas.draw()


def _interp(x1, y1, x2, y2):
    if x1 is None or x2 is None or y1 is None or y2 is None:
        return (), ()
    if x1 > x2:
        x1, x2, y1, y2 = x2, x1, y2, y1
    xp = np.arange(int(np.ceil(x1)), int(np.floor(x2))+1)
    yp = np.interp(xp, (x1, x2), (y1, y2))
    ind = np.logical_and(xp >= 1, xp <= len(state.stream))
    return xp[ind], yp[ind]


def fit_line(x1, y1, x2, y2):
    if x1 is None or x2 is None or y1 is None or y2 is None:
        return
    if x1 > x2:
        x1, x2, y1, y2 = x2, x1, y2, y1
    AX.plot((x1, x2), (y1, y2), 'o', color='C4')
    FIG.canvas.draw()
    recs = np.genfromtxt(sc['info'] + 'rec_xz.dat', dtype=None, encoding=None)
    rec_pos = {sp: (sx, sz) for sp, sx, sz in recs}
    x1 = rec_pos[int(round(x1))][0]
    x2 = rec_pos[int(round(x2))][0]
    xmax = sorted(rec_pos)[-1]
    xmax = rec_pos[xmax][0]
    a = (y1 - y2) / (x1 - x2)
    b = y1 - a * x1
    print('line: y = a * x + b')
    print(f'a = {a}')
    print(f'1 / a = {1/a}')
    print(f'b = {b}')


def press(event):
    if event.inaxes != AX:
        return
    print('pressed ' + event.key)

    if event.key == 'k':
        print(__doc__)
    elif event.key == 'up':
        opt.norm *= 1.2
        plot('update_wiggles')
    elif event.key == 'down':
        opt.norm /= 1.2
        plot('update_wiggles')
    elif event.key == 'c':
        opt.cut = not opt.cut
        plot('update_wiggles')
    elif event.key == 'p':
        opt.polarity = (opt.polarity + 2) % 3 - 1
        print(f'Polarity set to {opt.polarity}')
        plot('update_wiggles')
    elif event.key == 'h' and event.xdata is not None:
        i = int(round(event.xdata))
        opt.hide_trace[i] = not opt.hide_trace[i]
        plot('update_wiggles')
    elif event.key == 'n':
        opt.normalize = not opt.normalize
        print('Normalization turned ' + ('on' if opt.normalize else 'off'))
        plot('update_stream')
    elif event.key == 'f':
        opt.apply_filter = not opt.apply_filter
        print('Filter turned ' + ('on' if opt.apply_filter else 'off'))
        plot('update_stream')
    elif event.key == 'x':
        plot('new_shot')
    elif event.key in ('-', '+'):
        if event.key == '+':
            opt.nshot += 1
        else:
            opt.nshot -= 1
        if opt.nshot > conf.nshot_max:
            opt.nshot = 1
        if opt.nshot < 1:
            opt.nshot = conf.nshot_max
        new_shot()
    elif event.key == state.last_key and event.key in 'td':
        xp, yp = _interp(event.xdata, event.ydata, *state.last_xy)
        for x, y in zip(xp, yp):
            if event.key == 't':
                PICKS[opt.nshot][x] = y
            else:
                PICKS[opt.nshot].pop(x, None)
        plot('update_picks')
    elif event.key in 'ay':
        x = int(round(event.xdata))
        d = -1 if event.key == 'y' else 1
        autopick(x, d)
        plot('update_picks')
    elif event.key in 'lr':
        x1, x2 = AX.get_xlim()
        w = int(0.75 * (x2 - x1))
        d = -1 if event.key == 'l' else 1
        AX.set_xlim(x1 + d * w, x2 + d * w)
        FIG.canvas.draw()
    elif event.key in 'sw':
        print('Save picks')
        write_picks()
    elif event.key == state.last_key == 'm':
        fit_line(event.xdata, event.ydata, *state.last_xy)
    state.last_key = event.key
    state.last_xy = (event.xdata, event.ydata)


def set_filter(text):
    filt = [_.lower() for _ in text.split()]
    ftype = {'bp': 'bandpass', 'lp': 'lowpass', 'hp': 'highpass'}[filt[0]]
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
    opt.filter = dict(type=ftype, **kw)
    print(f'You chose filter {opt.filter}')


def set_picker(text):
    thres, t1, t2, t3, t4 = [float(_) for _ in text.split()]
    opt.pick = SimpleNamespace(thres=thres, t1=t1, t2=t2, t3=t3, t4=t4)
    print(f'You chose picker {opt.pick}')


def new_shot(text=None):
    if text is None:
#        text_box0.set_val(opt.nshot)
        text_box0.text = str(opt.nshot)
        text_box0.text_disp.remove()
        text_box0.text_disp = text_box0._make_text_disp(text_box0.text)
    else:
        opt.nshot = int(text)
    opt.hide_trace=defaultdict(bool)
    plot('new_shot')


PICKMODE = len(sys.argv) == 1

opt = SimpleNamespace(nshot=1, norm=1, normalize=True, cut=False,
                      apply_filter=False, filter=None,
                      polarity=1,
                      hide_trace=defaultdict(bool)
                      )
state = SimpleNamespace(last_key=None, last_xy=None,
                        stream=None, ostream=None,
                        pickmarkers=None, wiggles=None)

if not PICKMODE:
    if sys.argv[1] in ('-h', '--help'):
        print(__doc__)
        sys.exit()
    conf = SimpleNamespace(expr=sys.argv[1], stack=False, nshot_max=9999)
    if len(sys.argv) > 2:
        opt.nshot = int(sys.argv[2])
if PICKMODE:
    import seisutil as su

    with open('seisconf.json') as f:
        sc = json.load(f)
    PICKS = load_picks()
    shotpoints = su.read_shotpoints(sc['info'] + 'shotpoints.txt', verbose=True)
    spreads = su.read_spreads(sc['info'] + 'spreads.txt', verbose=True)
    conf = SimpleNamespace(
            expr=sc['expr'], stack=sc['stack'],
            filenumbers=su.read_filenumbers(sc['info'] + 'filenumbers.txt', verbose=True),
            shotpoints=shotpoints,
            nshot_max=sorted(shotpoints)[-1],
            spreads=spreads,
            ref_depth=sc['ref_depth'])


    # calculation for backshots and for writing fb_all.dat
    shotpoints_r = {r: sp for sp, r in shotpoints.items()}
    all_channels = sorted(set().union(*[set(spread['channels']) for spread in spreads]))
    missing_channels = [ch for ch in range(all_channels[0], all_channels[-1]+1) if ch not in all_channels]
    ch2rec = {ch: ch - np.sum(np.array(missing_channels) < ch) for ch in all_channels}
    rec2ch = {ch - np.sum(np.array(missing_channels) < ch): ch for ch in all_channels}
    spread_from_sp = {sp: spread for spread in spreads
                            for sp in spread['shotpoints']}
    backshots = {}
    for sp, rshot in shotpoints.items():
        for ch in spread_from_sp[sp]['channels']:
            spreadch = ch - spread_from_sp[sp]['channels'][0] + 1
            if shotpoints[sp] != ch2rec[ch] and ch in ch2rec and ch2rec[ch] in shotpoints_r:
                sp2 = shotpoints_r[ch2rec[ch]]
                if rec2ch[rshot] in spread_from_sp[sp2]['channels']:
                    spreadch2 = rec2ch[rshot] - spread_from_sp[sp2]['channels'][0] + 1
                    backshots[(sp, spreadch)] = (sp2, spreadch2)
    print()
    print('backshots')
    print(backshots)
    BS = backshots

FIG, AX = plt.subplots()
plt.subplots_adjust(bottom=0.1)
plot('first_call')
FIG.canvas.mpl_disconnect(FIG.canvas.manager.key_press_handler_id)
FIG.canvas.mpl_connect('key_press_event', press)

axbox0 = plt.axes([0.15, 0.02, 0.05, 0.04])
text_box0 = TextBox(axbox0, 'Shot', initial=str(opt.nshot))
text_box0.on_submit(new_shot)

axbox1 = plt.axes([0.3, 0.02, 0.25, 0.04])
filter0 = 'LP 100 2 1'
text_box1 = TextBox(axbox1, 'Filter', initial=filter0)
text_box1.on_submit(set_filter)
set_filter(filter0)
if PICKMODE:
    axbox2 = plt.axes([0.65, 0.02, 0.25, 0.04])
    picker0 = '0.8 -1e-3 1e-3 -2e-3 3e-3'
    text_box2 = TextBox(axbox2, 'Autopicker', initial=picker0)
    text_box2.on_submit(set_picker)
    set_picker(picker0)

print(__doc__)
plt.show()



