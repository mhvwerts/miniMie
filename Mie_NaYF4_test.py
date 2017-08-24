# -*- coding: utf-8 -*-
"""
@author: werts


Remember

# calculate extinction cross sections
sig_ext = (r_part*r_part)*pi*Qext    
sig_ext = Qext * pi * r**2

Extinction coefficient:
epsilon=1/(3.82e-25)*sig_ext

Scattering efficiency:
phisca=Qsca/Qext    

"""


import numpy as np
import matplotlib.pyplot as plt

from MiePlouc import Mie_spectrum

# some gold to see that everything works fine
d_nm=90.
wavelens =	np.linspace(250,1000,500)
Qext,Qsca = Mie_spectrum(wavelens, d_nm, mfp=False,mat='gold')
Qext_mfp, Qsca_mfp = Mie_spectrum(wavelens, d_nm, mfp=True,mat='gold')

plt.figure(1)
plt.clf()
plt.plot(wavelens,Qext,label='Q_ext')
plt.plot(wavelens,Qext_mfp,label=' mfp')
plt.ylabel('Q')
plt.xlabel('wavelength / nm')
plt.legend()


# test NaYF4 "spheres" using constant refractive index
n_mat = 0.5*(1.4748+1.4778) # values from  Sokolov Optics&Spectroscopy 2015
#TODO implement NaYF4 wavelength-dependent refractive index in MiePlouc
# select material  mat = "NaYF4"
d_nm=150.
Qext,Qsca = Mie_spectrum(wavelens, d_nm, mat=n_mat, n_medium=1.3855, mfp=False) # n-heptane
plt.figure(2)
plt.clf()
plt.plot(wavelens,Qext,label='Q_ext')
plt.plot(wavelens,Qsca,label='Q_sca')
plt.plot(wavelens,2e8*pow(wavelens,-4),label='lambda-4')
plt.ylabel('Q')
plt.xlabel('wavelength / nm')
plt.legend()

plt.show()
