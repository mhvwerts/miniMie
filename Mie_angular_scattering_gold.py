# -*- coding: utf-8 -*-
"""
Angular distribution of unpolarized Mie scattering intensity
for a series of gold nanoparticle diameters


Reproduces the left panel of Figure 1 in

Shortell, M. P.; Hewins, R. A.; Fernando, J. F. S.; Walden, S. L.; 
Waclawik, E. R.; Jaatinen, E. A., "Multi-Angle Fluorometer Technique for the 
Determination of Absorption and Scattering Coefficients of Subwavelength
Nanoparticles.", Opt. Express 2016, 24, 17090. doi:10.1364/OE.24.017090.
"""

import numpy as np
from numpy import cos, pi
import matplotlib.pyplot as plt

from miniMie import Mie_tetascan
from miniMie import JC_gold



wvln_nm = 500.
d_nm_list = [30., 60., 90., 120., 150., 180.]
fmts_list = ['C3--', 'C4-', 'C2--', 'C9-', 'C5--', 'C0-']
lw = 2.0

plt.figure(1)
plt.clf()

for d_nm, fmts in zip(d_nm_list, fmts_list):
    teta, Ipar, Iperp, Iunpol = Mie_tetascan(wvln_nm, d_nm,
                                             material = JC_gold(),
                                             n_medium=1.33,
                                             normalize = True)
    plt.polar(teta, Iunpol, fmts, linewidth = lw,
              label=f'{d_nm:.0f}nm')
    plt.polar(-teta, Iunpol, fmts, linewidth = lw)

# normalized Rayleigh
K = 1/(3*pi) 
Icos2 = K*(1+cos(teta)**2)
plt.polar(teta, Icos2, 'k:', label='Rayleigh')
plt.polar(-teta, Icos2, 'k:')

plt.title('Gold nanoparticles (norm. unpol. scatt.)')
plt.legend()
