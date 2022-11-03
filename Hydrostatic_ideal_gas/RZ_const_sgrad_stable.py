import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['idref'])
sys.path.append(os.environ['co'])
sys.path.append(os.environ['pl'])
from plotref import plotref
from scipy.special import expn

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import basic_constants as bc

# Define basic constants for a solar-like CZ
ri = bc.ri
rm = bc.rm
r = np.linspace(ri, rm, 1000)

Tm = bc.Tm
pm = bc.pm
rhom = bc.rhom
gam = bc.gamma
cp = bc.cp
cv = bc.cv
g = bc.G*bc.M/r**2
dgdr = -2.*g/r

colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y', 'b', 'r', 'g', 'm', 'c', 'k', 'y']
    
# Plot a range of polytropic indices
count = 0
firstplot = True

for k in np.linspace(0., 5., 7):
    # Here we START with the entropy gradient
    s = k*bc.cp*(r/rm - 1.)
    dsdr = k*bc.cp/rm*np.ones_like(r)
    d2sdr2 = np.zeros_like(r)
    expr = np.exp(-gam/(gam - 1.)*k*(r/rm - 1.))
    
    # Define the constant a0
    a0 = bc.G*bc.M/(bc.cp*Tm*rm)
    
    zeta = np.exp(k)*a0*(rm/r)*expn(2, k*r/rm) +\
            (1. - np.exp(k)*expn(2, k)*a0)
    zeta *= np.exp(k*(r/rm - 1.0))
    
    T = Tm*zeta
    p = pm*expr*zeta**(gam/(gam - 1.))
    rho = rhom*expr*zeta**(1./(gam - 1.))
    
    dlnT = dsdr/bc.cp - g/(bc.cp*T)
    dlnrho = 1./(gam - 1.)*(dlnT - dsdr/bc.cp)
    dlnp = dlnT + dlnrho
    
    d2lnT = d2sdr2/cp - dgdr/(cp*T) + g/(cp*T) * dlnT
    d2lnrho = (1./(gam - 1.))*(d2lnT - d2sdr2/cv)
    
#####    # Compute instead from arbitrary atmosphere
#    T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
#        arbitrary_atmosphere(r, s, dsdr, d2sdr2, g, dgdr, rm, Tm, pm, cp, gam)  

    # Plot the polytrope for this k-value  
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
    
axs[0,0].set_title('                  Constant entropy gradient ' + '(RZ)', **csfont)
    
plt.savefig('figures/RZ_const_sgrad.pdf')
plt.close()
