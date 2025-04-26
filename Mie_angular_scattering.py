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
                                         n_medium=1.33,
                                         normalize = True)

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

K = 1/(3*pi) # normalize Rayleigh to integral = 0.5
Icos2 = K*(1+cos(teta)**2)

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



# Not normalized for comparison: the overall shape should not change

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


