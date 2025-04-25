# -*- coding: utf-8 -*-
"""
@author: werts


Remember

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
from miniMie import JC_gold

d_nm = 60.
wavelens =	np.linspace(380, 1000, 500)
Qext,Qsca = Mie_spectrum(wavelens,
                         d_nm, material = JC_gold(),
                         n_medium=1.33)
Qext_mfp, Qsca_mfp = Mie_spectrum(wavelens,
                         d_nm, material = JC_gold(MFPradius_nm = d_nm/2),
                         n_medium=1.33)

plt.figure(1)
plt.clf()
plt.plot(wavelens, Qext, label='Q_ext')
plt.plot(wavelens, Qext_mfp, label=' mfp')
plt.ylabel('Q')
plt.xlabel('wavelength / nm')
plt.legend()

plt.figure(2)
plt.clf()
plt.plot(wavelens, Qsca, label='Q_sca')
plt.plot(wavelens, Qsca_mfp, label=' mfp')
plt.ylabel('Q')
plt.xlabel('wavelength / nm')
plt.legend()
 
plt.show()
