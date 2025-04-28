# -*- coding: utf-8 -*-
"""
Angular distribution of unpolarized Mie scattering intensity
for a series of gold nanoparticle diameters


Reproduces Figure 1 in

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

reldev_list = []
for d_nm, fmts in zip(d_nm_list, fmts_list):
    teta, Ipar, Iperp, Iunpol = Mie_tetascan(wvln_nm, d_nm,
                                             material = JC_gold(),
                                             n_medium = 1.33,
                                             normalize = True)
    # normalized Rayleigh, for relative deviation
    K1 = 16*pi/3 
    Icos2 = (1+cos(teta)**2)/K1
    reldev = Iunpol/Icos2 - 1
    reldev_list.append(reldev)
    
    plt.polar(teta, Iunpol, fmts, linewidth = lw,
              label=f'{d_nm:.0f}nm')
    plt.polar(-teta, Iunpol, fmts, linewidth = lw)

plt.polar(teta, Icos2, 'k:', label='Rayleigh')
plt.polar(-teta, Icos2, 'k:')

# plt.title('Gold nanoparticles (norm. unpol. scatt.)')
# plt.legend()



plt.figure(2)
plt.clf()
for d_nm, fmts, reldev in zip(d_nm_list, fmts_list, reldev_list):
    plt.plot(teta*180/pi, reldev*100.0, fmts, linewidth = lw,
             label=f'{d_nm:.0f}nm')
plt.xlabel('scattering angle (\u03b8)')
plt.ylabel('rel. deviation from Rayleigh (%)')
plt.xlim(0, 180)
plt.legend(frameon=False)

