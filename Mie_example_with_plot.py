# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 14:13:42 2016

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

d_nm=90.
wavelens =	np.linspace(380,1000,500)
Qext,Qsca = Mie_spectrum(wavelens, d_nm, mfp=False,mat='gold')
Qext_mfp, Qsca_mfp = Mie_spectrum(wavelens, d_nm, mfp=True,mat='gold')

plt.figure(1)
plt.clf()
plt.plot(wavelens,Qext,label='Q_ext')
plt.plot(wavelens,Qext_mfp,label=' mfp')
plt.ylabel('Q')
plt.xlabel('wavelength / nm')
plt.legend()

plt.figure(2)
plt.clf()
plt.plot(wavelens,Qsca,label='Q_sca')
plt.plot(wavelens,Qsca_mfp,label=' mfp')
plt.ylabel('Q')
plt.xlabel('wavelength / nm')
plt.legend()
 
plt.show()
