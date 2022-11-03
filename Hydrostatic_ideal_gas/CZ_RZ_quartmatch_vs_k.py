import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['idref'])
sys.path.append(os.environ['co'])
sys.path.append(os.environ['pl'])
from plotref import plotref
from arbitrary_atmosphere import arbitrary_atmosphere
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import basic_constants as bc

# Define basic constants
ri = bc.ri
rm = bc.rm
ro = bc.ro
nr = 1000
r = np.linspace(ri, ro, nr)

cp = bc.cp
Tm = bc.Tm
pm = bc.pm
rhom = bc.rhom
gam = bc.gamma
g = bc.G*bc.M/r**2
dgdr = -2.0*g/r

# First make a plot for different values of delta

colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y', 'b', 'r', 'g', 'm', 'c', 'k', 'y']

delta = 0.01*rm
kvals = np.linspace(0., 5., 7)
count = 0
firstplot = True
for k in kvals:
    s = np.zeros(nr)
    dsdr = np.zeros(nr)
    d2sdr2 = np.zeros(nr)
    for ir in range(nr):
        if r[ir] <= rm - delta:
            s[ir] = (8./15.)*k*cp*(delta/rm) + k*cp*(r[ir]/rm - 1.)
            dsdr[ir] = k*cp/rm
            d2sdr2[ir] = 0.
        elif r[ir] < rm:
            s[ir] = k*cp*(delta/rm)*((2./3.)*((r[ir] - rm)/delta)**3. -\
                    (1./5.)*((r[ir] - rm)/delta)**5.)
            dsdr[ir] = (k*cp/rm)*(1. - (1. - ((r[ir] - rm)/delta)**2.)**2.)
            d2sdr2[ir] = (4./delta)*(k*cp/rm)*\
                    (1. - ((r[ir] - rm)/delta)**2.)*((r[ir] - rm)/delta)
        else:
            s[ir] = 0.
            dsdr[ir] = 0.
            d2sdr2[ir] = 0.
    
    T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
        arbitrary_atmosphere(r, s, dsdr, d2sdr2, g, dgdr, rm, Tm, pm, cp,\
        gam)
    
    if firstplot:
        fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
            d2lnrho, label=r'$k=%.1f$' %k, color=colors[count])
        firstplot = False
    else:
        plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr, d2lnrho, \
                label=r'$k=%.1f$' %k, color=colors[count],\
                fig=fig, axs=axs)
    count += 1

plt.legend()
plt.tight_layout() 
    
axs[0,0].set_title('          CZ and RZ, delta/rm=0.01, quartic matching',\
        **csfont)
    
plt.savefig('figures/CZ_RZ_quartmatch_vs_k.pdf')
plt.close()
