# -*- coding: utf-8 -*-
"""
Angular distribution of Mie scattering intensity

@author: werts-moltech
"""

# TO DO: normalization option in Mie_tetascan


import numpy as np
import matplotlib.pyplot as plt

from miniMie import Mie_tetascan
from miniMie import JC_gold, Material



wvln_nm = 530.



# gold nanoparticle, JC_gold()

d_nm = 150.

teta, Ipar, Iperp, Iunpol = Mie_tetascan(wvln_nm, d_nm,
                                         material = JC_gold(),
                                         n_medium=1.33)

plt.figure(1)
plt.clf()
plt.polar(teta, Ipar, 'r', label='parallel')
plt.polar(-teta, Ipar, 'r')
plt.polar(teta, Iperp, 'b', label='perpendicular')
plt.polar(-teta, Iperp, 'b')
plt.polar(teta, Iunpol, 'k', label='unpolarized')
plt.polar(-teta, Iunpol, 'k')
plt.legend()



# Ludox, n=1.47-0.0j

d_nm = 30.

teta, Ipar, Iperp, Iunpol = Mie_tetascan(wvln_nm, d_nm,
                                         material = Material(1.47),
                                         n_medium=1.33)

plt.figure(2)
plt.clf()
plt.polar(teta, Ipar, 'r', label='parallel')
plt.polar(-teta, Ipar, 'r')
plt.polar(teta, Iperp, 'b', label='perpendicular')
plt.polar(-teta, Iperp, 'b')
plt.polar(teta, Iunpol, 'k', label='unpolarized')
plt.polar(-teta, Iunpol, 'k')
plt.legend()
