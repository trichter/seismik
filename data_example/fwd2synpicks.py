from glob import glob
import numpy as np


globexpr = 'INV/FWD/fd*.synth.asc'
soufname = 'GEO/sou_xz.dat'
recfname = 'GEO/rec_xz.dat'
fname_out = 'synthetic_picks.txt'

sou = np.loadtxt(soufname)[:,:2]
rec = np.loadtxt(recfname)[:,:2]

out = []

for fname in sorted(glob(globexpr)):
    picks = np.loadtxt(fname)
    assert picks[0, -1] == -1
    assert all(picks[1:, -1] == 1)
    picks = picks[:, [0, 3]]
    srcind = np.nonzero(np.abs(picks[0, 0] - sou[:, 1]) < 0.1)[0][0]
    srcnum = sou[srcind, 0]
    recind = np.nonzero(np.abs(picks[1:, 0] - rec[:, 1, np.newaxis]) < 0.1)[0]
    recnums = rec[recind, 0]
    o = np.vstack([np.ones(len(recnums))*srcnum, recnums, picks[1:,1] / 1000])
    out.append(o)

np.savetxt(fname_out, np.transpose(np.hstack(out)), fmt=('%d', '%d', '%.5f'), delimiter='  ')