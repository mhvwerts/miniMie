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

 (Maetzler 2002): 
       C. Maetzler, "MATLAB Functions for Mie Scattering and Absorption", 
       Research Report 2002-08, Institut fuer Angewandte Physik, 
       Universitaet Bern, 2002
       downloaded from http://www.iap.unibe.ch/publications/download/201/en/

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

from .materials import Material


# definition of Mie routines Mie_abcd and Mie

def Mie_abcd(m, x, nmax):
    """based on the MATLAB code by C. Maetzler, 2002
    Ref.: (Maetzler 2002)
    """
    n = arange(1,(nmax+1))*1.0
    nu = n+0.5
    z = m*x
    m2 = m*m
    sqx = sqrt(0.5 * pi / x)
    sqz = sqrt(0.5 * pi / z)
    bx = special.jv(nu, x) * sqx
    bz = special.jv(nu, z) * sqz
    yx = special.yv(nu, x) * sqx
    hx = (complex(1,0)*bx + complex(0,1)*yx)
    b1x = concatenate((array([(sin(x)/x)]),bx[0:(nmax-1)]))
    b1z = concatenate((array([(sin(z)/z)]),bz[0:(nmax-1)]))
    y1x = concatenate((array([(-cos(x)/x)]),yx[0:(nmax-1)]))
    h1x = complex(1,0)*b1x + complex(0,1)*y1x
    ax = x*b1x - n*bx
    az = z*b1z - n*bz
    ahx = x*h1x - n*hx
    an = (m2*bz*ax - bx*az)/(m2*bz*ahx - hx*az)
    bn = (bz*ax - bx*az)/(bz*ahx - hx*az)
    cn = (bx*ahx - hx*ax)/(bz*ahx - hx*az)
    dn = m*(bx*ahx - hx*ax)/(m2*bz*ahx - hx*az)
    return (an,bn,cn,dn)



def Mie(m, x):
    """The Mie routine adapted from Maetzler MATLAB code (Maetzler 2002).
    It calculates extinction, scattering and absorption cross sections, as well as
     the asymmetry parameter (avg cos theta) (and more later)
    for a single value of x, based on the complex refractive index contrast m.
    See the Maetzler document for the definition of x and m.
    The result is returned in the form of a 'tuple'. 
    Not all properties calculated in the original code are calculated here 
    (partial implementation).
    This function needs the Mie_abcd function. 
    """
    # check x==0  and avoid singularity
    if x==0:
        return (m.real, m.imag, 0., 0., 0., 0., 0.)
    nmax = int(round(2.0+x+4.0*x**(1./3.)))
    n1 = nmax - 1
    n = arange(1,nmax+1)
    cn = 2.0*n + 1.0
    c1n = n*(n+2.0)/(n+1.0)
    c2n = cn/n/(n+1.0)
    x2 = x*x
    (Mie_an,Mie_bn,Mie_cn,Mie_dn)=Mie_abcd(m, x, nmax) 
    anp = Mie_an.real
    anpp = Mie_an.imag
    bnp = Mie_bn.real
    bnpp = Mie_bn.imag
    g1=zeros((4,nmax))
    g1[0,0:n1]=anp[1:nmax]
    g1[1,0:n1]=anpp[1:nmax]
    g1[2,0:n1]=bnp[1:nmax]
    g1[3,0:n1]=bnpp[1:nmax]
    dn = cn*(anp+bnp)
    q = sum(dn)
    Qext = 2*q/x2
    en = cn*(anp*anp + anpp*anpp + bnp*bnp + bnpp*bnpp)
    q = sum(en)
    Qsca = 2*q/x2
    Qabs = Qext-Qsca
    asy1 = c1n*(anp*g1[0,:]+anpp*g1[1,:]+bnp*g1[2,:]+bnpp*g1[3,:])
    asy2 = c2n*(anp*bnp + anpp*bnpp)
    asy = 4.0/x2 * sum(asy1+asy2)/Qsca

    # We do not yet calculate the following properties,
    # contrary to the original Maetzler code:
    # Qb (backscatter), Qratio
  
    # return results as a tuple
    # follow the same format as clegett_mie.mie()
    # Qb and Qratio are not calculated by miniMie, set to NaN
    Qb = np.nan
    Qratio = np.nan
    return (m.real, m.imag, x, Qext, Qsca, Qabs, Qb, asy, Qratio)





def Mie_spectrum(wvln_nm, d_nm, material=Material(1.5), n_medium=1.33,
                 MieFun = Mie):
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
        resulttuple = MieFun(m, xco)
        Qext[idx] = resulttuple[3]
        Qsca[idx] = resulttuple[4]
        # not used:
        # Qabs[idx] = resulttuple[5]
        # asy[idx]  = resulttuple[7]
    return (Qext, Qsca)
    


# Execute the following only if run as a script.
# Testing code goes here.
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    d_nm=50.
    wavelens =	np.linspace(380, 1000, 500)
    Qext,Qsca = Mie_spectrum(wavelens, d_nm, mfp=False)
    Qext_mfp, Qsca_mfp = Mie_spectrum(wavelens, d_nm, mfp=True)
    
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

