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
r = np.linspace(bc.rm, bc.ro, 100)

# Plot a bunch of reference states corresponding to different values 
# of the polytropic index n
colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y', 'b', 'r', 'g', 'm', 'c', 'k',\
        'y']
nvals = np.array([1.5, 3., 5., 10., 100., 1000.])

count = 0
firstplot = True
for nval in nvals:
    # Define the constant a
    a = bc.G*bc.M/((nval + 1.)*bc.R*bc.Tm*bc.rm)
    zeta = a*bc.rm/r + (1. - a)

    T = bc.Tm*zeta
    p = bc.pm*zeta**(nval + 1.)
    rho = bc.rhom*zeta**nval

    dzeta = -a*bc.rm/r**2.
    d2zeta = 2.*a*bc.rm/r**3.
    dlnzeta = dzeta/zeta
    d2lnzeta = -dzeta**2./zeta**2. + d2zeta/zeta
    dlnT = dlnzeta
    dlnrho = nval*dlnzeta
    d2lnrho = nval*d2lnzeta

    dlnp = (nval + 1.)*dlnzeta

    s = bc.cv*(nval/bc.n0 - 1.)*(np.log(r/bc.rm) - \
            np.log(a + (1. - a)*r/bc.rm))
    dsdr = (nval/bc.n0 - 1.)*bc.cv/(r + (1. - a)*r**2./(a*bc.rm))

    # Plot this basic polytrope
    if firstplot:
        fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
                d2lnrho, label=r'$n=%.1f$' %nval, color=colors[count])
        firstplot = False
    else:
        plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr, d2lnrho,\
                label=r'$n=%.1f$' %nval,fig=fig, axs=axs,\
                color=colors[count])
    count += 1

axs[0,0].set_title('Convectively stable polytropes ' + r'$(N_\rho=3)$',\
        **csfont)
plt.tight_layout()
plt.legend()

plt.savefig('figures/CZ_polytrope_stable.pdf')
plt.close()
