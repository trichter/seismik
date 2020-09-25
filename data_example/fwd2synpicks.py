from glob import glob
import numpy as np
import seisutil as su

globexpr = 'INV/FWD2/fd*.synth.asc'
soufname = 'GEO/sou_xz.dat'
recfname = 'GEO/rec_xz.dat'
spredfname = 'GEO/spreads.txt'
fname_out = 'synthetic_picks.txt'

sou = np.loadtxt(soufname)[:,:2]
rec = np.loadtxt(recfname)[:,:2]
spreads = su.read_spreads(spredfname, verbose=True)
spread_from_sp = {sp: spread for spread in spreads
                  for sp in spread['shotpoints']}

out = []

for fname in sorted(glob(globexpr)):
    picks = np.loadtxt(fname)
    assert picks[0, -1] == -1
    assert all(picks[1:, -1] == 1)
    picks = picks[:, [0, 3]]
    for srcind in np.nonzero(np.abs(picks[0, 0] - sou[:, 1]) < 0.1)[0]:
        srcnum = sou[srcind, 0]
        rec0 = spread_from_sp[srcnum]['channels'][0]
        recind = np.nonzero(np.abs(picks[1:, 0] - rec[:, 1, np.newaxis]) < 0.1)[0]
        recnums = rec[recind, 0] - rec0 + 1
        o = np.vstack([np.ones(len(recnums))*srcnum, recnums, picks[1:,1] / 1000])
        out.append(o)

np.savetxt(fname_out, np.transpose(np.hstack(out)), fmt=('%d', '%d', '%.5f'), delimiter='  ')
