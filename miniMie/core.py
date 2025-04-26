# -*- coding: utf-8 -*-
"""
High-level Mie calculation routines


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


import numpy as np
from numpy import pi, cos

from .clegett_mie import mie
from .clegett_mie import mie_s12
from .materials import Material



def Mie_spectrum(wvln_nm, d_nm, material=Material(1.5), n_medium=1.33):
    """
    Generate extinction and scattering spectra using Mie theory

    Parameters
    ----------
    wvln_nm : np.ndarray(dtype=float)
        Vector of desired wavelengths [nm].
    d_nm : float
        Nanosphere diameter [nm].
    material : instance of Material, optional
        For evaluation of dielectric function. The default is Material(1.5).
    n_medium : float, optional
        Refractive index of medium. The default is 1.33.

    Returns
    -------
    Qext : np.ndarray(dtype=float)
        Vector of extinction efficiencies Q_ext
    Qsca : np.ndarray(dtype=float)
        Vector of scattering efficiencies Q_sca

    """
    
    r_sphere=(d_nm*1e-9)/2
    wvln = wvln_nm*1e-9

    # get dielectric function
    ncmplx_wvln = material.get_ncmplx_vector(wvln_nm)
    
    # CALCULATION of spectra
    Npts = len(wvln)
    Qext = np.zeros(Npts)
    Qsca = np.zeros(Npts)
    # not used:
    # Qabs = zeros(Npts)
    # asy = zeros(Npts)
    for idx in range(Npts):
        xco = (2*pi*n_medium*r_sphere)/wvln[idx]
        m = ncmplx_wvln[idx]/n_medium
        resulttuple = mie(m, xco)
        Qext[idx] = resulttuple[3]
        Qsca[idx] = resulttuple[4]
        # not used:
        # Qabs[idx] = resulttuple[5]
        # asy[idx]  = resulttuple[7]
    return (Qext, Qsca)
    


def Mie_tetascan(wvln_nm, d_nm, material=Material(1.5), n_medium=1.33,
                 Npts = 400):
    """
    Computation of Mie Power Scattering from 0° to 180°

    Parameters
    ----------
    wvln_nm : float
        Vacuum wavelength of incident light.
    d_nm : float
        Nanosphere diameter [nm].
    material : instance of Material, optional
        For evaluation of dielectric function. The default is Material(1.5).
    n_medium : float, optional
        Refractive index of medium. The default is 1.33.
    Npts : int, optional
        Number of points to be evaluated. The default is 400.

    Returns
    -------
    teta : np.ndarray(dtype=float)
        Scattering angle (in radians).
    SL : np.ndarray(dtype=float)
        Scattered intensity, parallel polarization.
    SR : np.ndarray(dtype=float)
        Scattered intensity, perpendicular polarization.
    SU : np.ndarray(dtype=float)
        Scattered intensity, unpolarized, 0.5*(SL+SR).

    """
    r_sphere=(d_nm*1e-9)/2 
    wvln = wvln_nm*1e-9   

    # get dielectric function
    ncmplx = material.get_ncmplx_vector(np.array([wvln_nm]))[0]
    
    x = (2*pi*n_medium*r_sphere)/wvln
    m = ncmplx/n_medium
    
    teta = np.linspace(0, pi, Npts)
    SL = np.zeros_like(teta)
    SR = np.zeros_like(teta)
    for j, tetaval in enumerate(teta):
        u = cos(tetaval)
        S1, S2 = mie_s12(m, x, u)
        SL[j] = (S1*S1.conj()).real
        SR[j] = (S2*S2.conj()).real
    SU = 0.5*(SL+SR)
    
    return (teta, SL, SR, SU)




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

