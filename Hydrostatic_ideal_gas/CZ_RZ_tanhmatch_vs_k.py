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
    d2sdr2 = -k*cp/rm*(1./(2.*delta))*(1./np.cosh((r - rm)/delta))**2.
    dsdr = k*cp/rm*(0.5*(1.0 - np.tanh((r - rm)/delta)))
    s = k*cp*(0.5*((r/rm - 1.0) -\
            (delta/rm)*np.log(np.cosh((r - rm)/delta))))
    
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
    
axs[0,0].set_title('          CZ and RZ, delta/rm=0.01, tanh matching',\
        **csfont)
    
plt.savefig('figures/CZ_RZ_tanhmatch_vs_k.pdf')
plt.close()
