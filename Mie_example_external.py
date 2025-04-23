# -*- coding: utf-8 -*-
"""
@author: werts


Compare spectra calculated using original miniMie.Mie() function and
Chip Legett's clegget_mie.mie().


Remember:

# calculate extinction cross sections
sig_ext = Qext * pi * r**2

Extinction coefficient:
epsilon=1/(3.82e-25)*sig_ext

Scattering efficiency:
phisca=Qsca/Qext    
"""

import numpy as np
import matplotlib.pyplot as plt

from miniMie import Mie_spectrum
import clegett_mie

d_nm = 60.
wavelens =	np.linspace(380, 1000, 500)
Qext, Qsca = Mie_spectrum(wavelens,
                         d_nm, mat='gold', mfp=False,
                         n_medium=1.33)
Qext_cl, Qsca_cl = Mie_spectrum(wavelens,
                         d_nm, mat='gold', mfp=False,
                         n_medium=1.33,
                         MieFun = clegett_mie.mie)

# Test if results calculated by original miniMie and clegett_mie are
# identical.
assert np.allclose(Qext, Qext_cl)
assert np.allclose(Qsca, Qsca_cl)

print("Yes. miniMie.Mie() and clegett_mie.mie() give identical results with Mie_spectrum().")

plt.figure(1)
plt.clf()
plt.plot(wavelens, Qext, label='Q_ext')
plt.plot(wavelens, Qext_cl, ':', label=' clegett_mie')
plt.ylabel('Q')
plt.xlabel('wavelength / nm')
plt.legend()

plt.figure(2)
plt.clf()
plt.plot(wavelens, Qsca, label='Q_sca')
plt.plot(wavelens, Qsca_cl, ':', label=' clegett_mie')
plt.ylabel('Q')
plt.xlabel('wavelength / nm')
plt.legend()
 
plt.show()
