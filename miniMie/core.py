# -*- coding: utf-8 -*-
"""
Mie calculation routines:
output Qext and Qsca for input vector of wavelengths

M. H. V. Werts, CNRS, Université d'Angers, France.

Read the license text at the end of this file before using this software.



 Literature references

 (Bohren and Huffman 1983):
       C.F. Bohren and D.R. Huffman, "Absorption and Scattering of Light
       by Small Particles", Wiley Interscience: New York, 1983.

 (Mätzler 2002): 
       C. Mätzler, "MATLAB Functions for Mie Scattering and Absorption, 
       Version 2", Research Report 2002-11, Institut für Angewandte Physik, 
       Universität Bern, 2002
       https://doi.org/10.7892/boris.146550

 (Johnson and Christy 1972):
       P.B. Johnson and R.W. Christy, "Optical Constants of the Noble 
       Metals", Phys. Rev. B. 1972, (6), 4370

 (Haiss et al. 2007):
       W. Haiss, N.T.K. Thanh, J. Aveyard, D.G. Fernig, "Determination of
       size and concentration of gold nanoparticles from UV-vis spectra.",
       Anal. Chem. 2007, (79), 4215

 (Kreibig 1974):
       U. Kreibig, "Electronic properties of small silver particles: the
       optical constants and their temperature dependence",
       J. Phys. F: Metal Phys. 1974, 4, 999

 (Murata and Tanaka 2010):
       K.I. Murata and H. Tanaka, "Surface-wetting effects on the liquid-liquid
       transition of a single-component molecular liquid.",
       Nature Commun. 2010, 1, 16

"""


#TODO ADD BENCHMARKING RESULTS in particular for gold, silver particles

import numpy as np
from numpy import array, arange, zeros, concatenate
from numpy import sqrt, sin, cos, pi
from scipy import special
from scipy import interpolate

from .clegett_mie import mie
from .materials import Material




def Mie_spectrum(wvln_nm, d_nm, material=Material(1.5), n_medium=1.33):
    """generate extinction and scattering spectra 
    for a sphere of diameter d_nm in medium with refractive index n_medium
    sampled on the wavelengths specified in wvln_nm
    
    output: 2-tuple of numpy vectors (extinction and scattering)
    
    The kwarg `MieFun` enables to supply an external Mie calculation function,
    e.g. from another Mie library.
    """ 
    
    r_sphere=(d_nm*1e-9)/2 # allows use of both SI unit-based values
    wvln = wvln_nm*1e-9   # and simple floats (already divided by nm for example)

    # get dielectric function
    ncmplx_wvln = material.get_ncmplx_vector(wvln_nm)
    
    # CALCULATION of spectra
    Npts = len(wvln)
    Qext = zeros(Npts)
    Qsca = zeros(Npts)
    # not used:
    # Qabs = zeros(Npts)
    # asy = zeros(Npts)
    for idx in range(Npts):
        xco = (2*pi*n_medium*r_sphere)/wvln[idx]
        m = ncmplx_wvln[idx]/n_medium # use bulk dielectric function
        resulttuple = mie(m, xco)
        Qext[idx] = resulttuple[3]
        Qsca[idx] = resulttuple[4]
        # not used:
        # Qabs[idx] = resulttuple[5]
        # asy[idx]  = resulttuple[7]
    return (Qext, Qsca)
    



#
#Copyright M. H. V. Werts, 2013-2025
#
#martinus point werts à univ-angers point fr
#
#
#This software is a computer program whose purpose is to calculate 
#the optical cross sections of nanoparticles using Mie theory.
#
#This software is governed by the CeCILL  license under French law and
#abiding by the rules of distribution of free software.  You can  use, 
#modify and/ or redistribute the software under the terms of the CeCILL
#license as circulated by CEA, CNRS and INRIA at the following URL
#"https://cecill.info". See also the file 'LICENSE' distributed with this
#software.
#
#As a counterpart to the access to the source code and  rights to copy,
#modify and redistribute granted by the license, users are provided only
#with a limited warranty  and the software's author,  the holder of the
#economic rights,  and the successive licensors  have only  limited
#liability. 
#
#In this respect, the user's attention is drawn to the risks associated
#with loading,  using,  modifying and/or developing or reproducing the
#software by the user in light of its specific status of free software,
#that may mean  that it is complicated to manipulate,  and  that  also
#therefore means  that it is reserved for developers  and  experienced
#professionals having in-depth computer knowledge. Users are therefore
#encouraged to load and test the software's suitability as regards their
#requirements in conditions enabling the security of their systems and/or 
#data to be ensured and,  more generally, to use and operate it in the 
#same conditions as regards security. 
#
#The fact that you are presently reading this means that you have had
#knowledge of the CeCILL license and that you accept its terms.
#
#

