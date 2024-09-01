# Update (02/21/2023):
# Should only need to run once, then just read the binary data
# python make_binary.py

import numpy as np
import pickle

fname = 'avgong1995-2009rls'
f = open(fname + '.txt', 'r')
lines = f.readlines()
colat = []
radius = []
om = []
om_err = []

# Read in the (flattened) arrays
for line in lines:
    li = line.split()
    colat.append(float(li[0]))
    radius.append(float(li[1]))
    om.append(float(li[2]))
    om_err.append(float(li[3]))

# Convert lists to arrays
tt_2d = np.array(colat)
rr_2d = np.array(radius)
om = np.array(om)
om_err = np.array(om_err)

# This is the shape of the 2D array
nt = 49
nr = 51

# Convert flattened arrays to 2D arrays
tt_2d = tt_2d.reshape((nt, nr))
rr_2d = rr_2d.reshape((nt, nr))

om = om.reshape((nt, nr))
om_err = om_err.reshape((nt, nr))

# Invert the axes of the arrays to order them like Rayleigh AZ_Avgs data
tt_2d = tt_2d[::-1, ::-1]
rr_2d = rr_2d[::-1, ::-1]
om = om[::-1, ::-1]
om_err = om_err[::-1, ::-1]

# other grid info
tt = tt_2d[:, 0]
rr = rr_2d[0,:]

# cut the data off at a bottom radius (avoid over-interpretation)
# and at a high latitude

#  don't cut the high latitudes
itcut = nt
rbot = 0.5 # radius
irbot = np.argmin(np.abs(rr - rbot))

# cut the data
tt = tt[:itcut+1]
rr = rr[:irbot+1]

nt = len(tt)
nr = len(rr)

tt_2d = tt_2d[:itcut+1, :irbot+1]
rr_2d = rr_2d[:itcut+1, :irbot+1]
om = om[:itcut+1, :irbot+1]
om_err = om_err[:itcut+1, :irbot+1]

# more grid info (now that everything is cut)
cost = np.cos(tt)
sint = np.sin(tt)

cost_2d = np.cos(tt_2d)
sint_2d = np.sin(tt_2d)

xx = rr_2d*sint_2d
zz = rr_2d*cost_2d

fsave = open(fname + '_full.pkl', 'wb')
pickle.dump({
    'om': om, 'om_err': om_err,
    'nt': nt, 'nr': nr, 
    'tt': tt, 'rr': rr, 'sint': sint, 'cost': cost,
    'rr_2d': rr_2d, 'tt_2d': tt_2d, 'sint_2d': sint_2d, 'cost_2d': cost_2d,
    'xx': xx, 'zz': zz
    },
        fsave, protocol=4)
fsave.close()
f.close()
