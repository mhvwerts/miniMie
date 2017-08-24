# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 14:13:42 2016

@author: werts
"""

import numpy as np
import matplotlib.pyplot as plt

from Mie_mfp import Mie_gold_water

d_nm=50.
wavelens =	np.linspace(400,1000,500)
Qext,Qsca = Mie_gold_water(d_nm, wavelens)

plt.clf()
plt.plot(wavelens,Qext,label='Q_ext')
plt.plot(wavelens,Qsca,label='Q_sca')
plt.ylabel('Q')
plt.xlabel('wavelength / nm')
plt.show()

# à toi de convertir Q_ext en coefficient d'extinction!
# la réponse est dans l'article Analyst 2013
