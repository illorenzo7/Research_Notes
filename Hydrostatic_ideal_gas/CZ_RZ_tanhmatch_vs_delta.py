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
r = np.linspace(ri, ro, 1000)

cp = bc.cp
Tm = bc.Tm
pm = bc.pm
rhom = bc.rhom
gam = bc.gamma

# First make a plot for different values of delta

k = 1.
colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y']

deltas = rm*np.array([0.001, 0.01, 0.03, 0.05, 0.1, 0.15])
count = 0
firstplot = True
for delta in deltas:
    d2sdr2 = -k*cp/rm*(1./(2.*delta))*(1./np.cosh((r - rm)/delta))**2.
    dsdr = k*cp/rm*(1./2.)*(1. - np.tanh((r - rm)/delta))
    s = k*cp/rm*(1./2.)*((r - rm) - delta*np.log(np.cosh((r - rm)/delta)))
    g = bc.G*bc.M/r**2
    dgdr = -2.0*g/r
    
    T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
        arbitrary_atmosphere(r, s, dsdr, d2sdr2, g, dgdr, rm, Tm, pm,\
        cp, gam)
    
    if firstplot:
        fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
            d2lnrho, label=r'$\delta/r_{\rm{m}}=%.3f$' %(delta/rm),\
            color=colors[count])
        firstplot = False
    else:
        plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr, d2lnrho, \
                label=r'$\delta/r_{\rm{m}}=%.3f$' %(delta/rm),\
                color=colors[count], fig=fig, axs=axs)
    count += 1

plt.legend()
plt.tight_layout() 
    
axs[0,0].set_title('            CZ and RZ, k=1.0, tanh matching', **csfont)
    
plt.savefig('figures/CZ_RZ_tanhmatch_vs_delta.pdf')
plt.close()
