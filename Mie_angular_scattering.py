# -*- coding: utf-8 -*-
"""
Angular distribution of Mie scattering intensity

@author: werts-moltech
"""



import numpy as np
from numpy import cos, pi
import matplotlib.pyplot as plt

from miniMie import Mie_tetascan
from miniMie import JC_gold, Material



wvln_nm = 530.



# gold nanoparticle, JC_gold()

d_nm = 150.

teta, Ipar, Iperp, Iunpol = Mie_tetascan(wvln_nm, d_nm,
                                         material = JC_gold(),
                                         n_medium = 1.33,
                                         normalize = True)

mu = np.cos(teta)
total_teta = 2 * np.pi * np.trapezoid(Iunpol*np.sin(teta), teta)
total_mu = -2 * np.pi * np.trapezoid(Iunpol, mu) # minus because of inversed limits
print('Full sphere integral: ', total_teta, total_mu )
print()

plt.figure(1)
plt.clf()

plt.polar(teta, Ipar, 'r', label='parallel')
plt.polar(-teta, Ipar, 'r')
plt.polar(teta, Iperp, 'b', label='perpendicular')
plt.polar(-teta, Iperp, 'b')
plt.polar(teta, Iunpol, 'g', label='unpolarized')
plt.polar(-teta, Iunpol, 'g')

plt.title(f'{d_nm:.0f}nm gold particles (normalized)')
plt.legend()



# Tiny silica particles, n=1.47-0.0j
# Mie and compare to Rayleigh scattering

d_nm = 10.


teta, Ipar, Iperp, Iunpol = Mie_tetascan(wvln_nm, d_nm,
                                         material = Material(1.47),
                                         n_medium=1.33,
                                         normalize = True)

mu = np.cos(teta)
total_teta = 2 * np.pi * np.trapezoid(Iunpol*np.sin(teta), teta)
total_mu = -2 * np.pi * np.trapezoid(Iunpol, mu) # minus because of inversed limits
print('Full sphere integral: ', total_teta, total_mu )
print()

K1 = 16*pi/3 # normalize Rayleigh integral to unity for full sphere
Icos2 = (1+cos(teta)**2)/K1

plt.figure(2)
plt.clf()

plt.polar(teta, Ipar, 'r', label='parallel')
plt.polar(-teta, Ipar, 'r')
plt.polar(teta, Iperp, 'b', label='perpendicular')
plt.polar(-teta, Iperp, 'b')
plt.polar(teta, Iunpol, 'g', label='unpolarized')
plt.polar(-teta, Iunpol, 'g')


plt.polar(teta, Icos2, 'k:', label='Rayleigh')
plt.polar(-teta, Icos2, 'k:')

plt.title(f'{d_nm:.0f}nm silica particles (normalized)')
plt.legend()





# Unnormalized calculation for comparison: the overall shape should not change

teta, Ipar, Iperp, Iunpol = Mie_tetascan(wvln_nm, d_nm,
                                         material = Material(1.47),
                                         n_medium=1.33,
                                         normalize = False)

plt.figure(3)
plt.clf()

plt.polar(teta, Ipar, 'r', label='parallel')
plt.polar(-teta, Ipar, 'r')
plt.polar(teta, Iperp, 'b', label='perpendicular')
plt.polar(-teta, Iperp, 'b')
plt.polar(teta, Iunpol, 'g', label='unpolarized')
plt.polar(-teta, Iunpol, 'g')

plt.title(f'{d_nm:.0f}nm silica particles (not normalized)')
plt.legend()


