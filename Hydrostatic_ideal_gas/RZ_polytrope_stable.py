import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['idref'])
sys.path.append(os.environ['co'])
sys.path.append(os.environ['pl'])
from plotref import plotref

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import basic_constants as bc

# Define radius of a hundred grid points
ri = bc.ri
rm = bc.rm
Tm = bc.Tm
pm = bc.pm
rhom = bc.rhom
r = np.linspace(ri, rm, 100)

colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y']
# Plot a range of polytropic indices
count = 0
firstplot = True
for n in np.arange(1.5, 20, 3):
    # Define the constant a
    a = bc.G*bc.M/((n + 1.)*bc.R*Tm*rm)
    zeta = a*rm/r + (1. - a)
    
    T = Tm*zeta
    p = pm*zeta**(n + 1.0)
    rho = rhom*zeta**n
    s = bc.cv*(n/bc.n0 - 1.)*(np.log(r/rm) - np.log(a + (1. - a)*r/rm))
    dsdr = (n/bc.n0 - 1.)*bc.cv/(r + (1. - a)*r**2./(a*rm))
    
    dzeta = -a*bc.rm/r**2.
    d2zeta = 2.*a*bc.rm/r**3.
    dlnzeta = dzeta/zeta
    d2lnzeta = -dzeta**2./zeta**2. + d2zeta/zeta
    dlnT = dlnzeta
    dlnp = (n + 1.)*dlnzeta
    dlnrho = n*dlnzeta
    d2lnrho = n*d2lnzeta

    # Plot the polytrope for this n-value
    if firstplot:
        fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
                d2lnrho, label=r'$n=%.1f$' %n, color=colors[count])
        firstplot = False
    else:
        plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr, d2lnrho,\
                label=r'$n=%.1f$' %n,fig=fig, axs=axs,\
                color=colors[count])
    count += 1

plt.legend()
plt.tight_layout() 
    
axs[0,0].set_title('    Stable polytropes ' + '(radiative zone)', **csfont)
plt.savefig('figures/RZ_polytrope_stable.pdf')
plt.close()
