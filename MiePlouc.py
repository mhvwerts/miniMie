# -*- coding: utf-8 -*-
"""
Mie calculation routines:
output Qext and Qsca for input vector of wavelengths

M. H. V. Werts

"""

#TODO BENCHMARK? in particular for gold, silver particles

from numpy import sqrt,array,arange,pi,concatenate,sin,cos,zeros
from scipy import special
from scipy import interpolate
import numpy as np



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

    # we do not calculated the following, contrary to original code:
    # Qb backscatter, Qratio
  
    # return results as a tuple
    return (m.real, m.imag, x, Qext, Qsca, Qabs, asy)



def ncmplx_mfpcorr(ncmplx_bulk, radius, waveln, FV, OMP, OM0):
    """mean free path correction
    
    the input takes the "bulk" complex dielectric function
    (as a vector)
    together with a vector of the wavelengths
    and material parameters
    
    returns the MFP-corrected dielectric function
    
    
    adapted from Haiss FORTRAN code (Haiss et al. 2007),
    which is in cgs units, which was maintained here
    We used this code with minimal changes 
    (this is why there are UPPERCASE variable names...)
    radius, waveln in nanometers
    FV in cm/s, OMP, OM0 in 1E-14 Hz
    """
    rn = ncmplx_bulk.real
    rk = ncmplx_bulk.imag
    CL = 2.998E+10
    # 	Calculate EPS1 and EPS2 from rn and rk:
    EPS1 = rn*rn - rk*rk
    EPS2 = 2.*rn*rk
    # 	Calculate OM and A1 and A2: 
    #	CL speed of light in cm/s
    OM = (2.*pi*CL/(waveln*1.E-7))/1.E+14 
    #   OMP: bulk plasma frequency in Hz divided by 1E+14
    #  	OM0: collision frequency in Hz/1E+14
    A1 = 1.-(OMP*OMP/(OM*OM + OM0*OM0)) 
    A2 = OMP*OMP*OM0/(OM*(OM*OM + OM0*OM0))
    #  	Contribution of the bond electrons to n (B1) and k (B2): 
    B1 = EPS1 - A1
    B2 = EPS2 - A2
    #   Calculate R dependent OM0 (OM0R)
    OM0R = OM0 + (FV/(radius*1.E-7))/1.E+14
    #   Calculate R dependent contributions of the free electrons:
    A1R = 1.-(OMP*OMP/(OM*OM + OM0R*OM0R)) 
    A2R = OMP*OMP*OM0R/(OM*(OM*OM + OM0R*OM0R))
    #   Calculate R dependent EPS (EPS1R and EPS2R)
    EPS1R = A1R + B1
    EPS2R = A2R + B2
    #	Reconvert EPS1R and EPS2R back to n and k:
    rnr = sqrt((A1R + B1)/2. + \
        sqrt((A1R/2.+B1/2.)*(A1R/2.+B1/2.)+(A2R/2.+B2/2.)*(A2R/2.+B2/2.))) 
    rkr = sqrt(-(A1R + B1)/2. + \
    sqrt((A1R/2.+B1/2.)*(A1R/2.+B1/2.)+(A2R/2.+B2/2.)*(A2R/2.+B2/2.)))
    # reconstruct complex index   
    ncmplx_corr = complex(1,0) * rnr + complex(0,1) * rkr
    return ncmplx_corr



def get_ncmplx_vector(wvln_nm, mat, MFPradius_nm = None):
    """generate a vector of complex dielectric function of a material
    sampled to the wavelengths (nm) in the input vector
    
    In:
        wvln_nm    vector of wavelengths (float)
                   for which dielectric function
                   should be calculated
        mat        (string) describing the material
        MFPdiam_nm    (float) if specified it applies a mean free path
                    correction (only available for gold and silver,
                    ignored elsewhere). Diameter of the particle
                    is specified in nm
    """
    
    # complex dielectric functions of gold and silver by Johnson and Christy
    if mat == 'gold':
        E = array([0.64,0.77,0.89,1.02,1.14,1.26,1.39,1.51,1.64,1.76,1.88,
                   2.01,2.13,2.26,2.38,2.50,
                   2.63,2.75,2.88,3,3.12,3.25,3.37,3.5,3.62,3.74,3.87,3.99,
                   4.12,4.24,4.36,4.49,4.61,4.74,4.86,
                   4.98,5.11,5.23,5.36,5.48,5.6,5.73,5.85,5.98,6.1,6.22,
                   6.35,6.47,6.6])
        n = array([0.92,0.56,0.43,0.35,0.27,0.22,0.17,0.16,0.14,0.13,0.14,
                   0.21,0.29,0.43,0.62,1.04,
                   1.31,1.38,1.45,1.46,1.47,1.46,1.48,1.50,1.48,1.48,1.54,
                   1.53,1.53,1.49,1.47,1.43,1.38,1.35,
                   1.33,1.33,1.32,1.32,1.30,1.31,1.30,1.30,1.30,1.30,1.33,
                   1.33,1.34,1.32,1.28])
        k = array([13.78,11.21,9.519,8.145,7.150,6.350,5.663,5.083,4.542,
                   4.103,3.697,3.272,
                   2.863,2.455,2.081,1.833,1.849,1.914,1.948,1.958,1.952,
                   1.933,1.895,1.866,1.871,1.883,1.898,
                   1.893,1.889,1.878,1.869,1.847,1.803,1.749,1.688,1.631,
                   1.577,1.536,1.497,1.460,1.427,1.387,
                   1.350,1.304,1.277,1.251,1.226,1.203,1.188])
        # following values are from Haiss et al.
        FV = 1.4E8 # Fermi velocity in cm/s - needed by mean free path correction
        OMP = 138. # plasma frequency in Hz/1E+14
        OM0 = 0.333 # collision frequency in Hz/1E+14  
        ncmplx = complex(1,0) * n + complex(0,1) * k # construct complex vector
        ncmplx_interpol=interpolate.interp1d(E, ncmplx, kind='cubic')
        # wavelength-sampled dielectric function
        ncmplx_wvln0 = ncmplx_interpol(1240.0/(wvln_nm))
        if not MFPradius_nm == None:
            ncmplx_wvln = ncmplx_mfpcorr(ncmplx_wvln0, MFPradius_nm, wvln_nm,
                                         FV, OMP, OM0)   
        else:
            ncmplx_wvln = ncmplx_wvln0
    elif mat == 'silver':
        E = array([0.64,0.77,0.89,1.02,1.14,1.26,1.39,1.51,1.64,1.76,1.88,
                   2.01,2.13,
                   2.26,2.38,2.50,2.63,2.75,2.88,3,3.12,3.25,3.37,3.5,3.62,
                   3.74,3.87,3.99,4.12,
                   4.24,4.36,4.49,4.61,4.74,4.86,4.98,5.11,5.23,5.36,5.48,
                   5.6,5.73,5.85,5.98,6.1,6.22,6.35,6.47,6.6])
        n = array([0.24,0.15,0.13,0.09,0.04,0.04,0.04,0.04,0.03,0.04,0.05,
                   0.06,0.05,
                   0.06,0.05,0.05,0.05,0.04,0.04,0.05,0.05,0.05,0.07,0.1,0.14,
                   0.17,0.81,1.13,1.34,
                   1.39,1.41,1.41,1.38,1.35,1.33,1.31,1.3,1.28,1.28,1.26,
                   1.25,1.22,1.20,1.18,1.15,1.14,1.12,1.10,1.07])
        k = array([14.08,11.85,10.10,8.828,7.795,6.692,6.312,5.727,5.242,
                   4.838,4.483,
                   4.152,3.858,3.586,3.324,3.093,2.869,2.657,2.462,2.275,
                   2.07,1.864,1.657,1.419,
                   1.142,0.829,0.392,0.616,0.964,1.161,1.264,1.331,1.372,
                   1.387,1.393,1.389,1.378,
                   1.367,1.357,1.344,1.342,1.336,1.325,1.312,1.296,1.277,
                   1.255,1.232,1.212])
        # following values are from Kreibig 1974 + Murata, Tanaka 2010
        FV = 1.4E8 # Fermi velocity in cm/s - needed by mean free path correction
        OMP = 137. # plasma frequency in Hz/1E+14
        OM0 = 0.27 # collision frequency in Hz/1E+14
        ncmplx = complex(1,0) * n + complex(0,1) * k # construct complex vector
        ncmplx_interpol=interpolate.interp1d(E, ncmplx, kind='cubic')
        # wavelength-sampled dielectric function
        ncmplx_wvln0 = ncmplx_interpol(1240.0/(wvln_nm))
        if not MFPradius_nm == None:
            ncmplx_wvln = ncmplx_mfpcorr(ncmplx_wvln0, MFPradius_nm, wvln_nm,
                                         FV, OMP, OM0)   
        else:
            ncmplx_wvln = ncmplx_wvln0
    else:
        raise ValueError("material 'mat' not known (lowercase only!)")
    return ncmplx_wvln



def Mie_spectrum(wvln_nm, d_nm, mat="gold", n_medium=1.33, mfp = True):
    """generate extinction and scattering spectra 
    for a sphere of diameter d_nm in medium with refractive index n_medium
    sampled on the wavelengths specified in wvln_nm
    
    output: 2-tuple of numpy vectors (extinction and scattering)
    """ 
    
    r_sphere=(d_nm*1e-9)/2 # allows use of both SI unit-based values
    wvln = wvln_nm*1e-9   # and simple floats (already divided by nm for example)

    # get dielectric function
    if mfp:
        ncmplx_wvln = get_ncmplx_vector(wvln_nm, mat, MFPradius_nm = d_nm/2.)
    else:
        ncmplx_wvln = get_ncmplx_vector(wvln_nm, mat)
    
    # CALCULATION of spectra
    Npts = len(wvln)
    Qext = zeros(Npts)
    Qsca = zeros(Npts)
    Qabs = zeros(Npts)
    asy = zeros(Npts)
    for idx in range(Npts):
        xco = (2*pi*n_medium*r_sphere)/wvln[idx]
        m = ncmplx_wvln[idx]/n_medium # use bulk dielectric function
        resulttuple = Mie(m, xco)
        Qext[idx] = resulttuple[3]
        Qsca[idx] = resulttuple[4]
        Qabs[idx] = resulttuple[5]
        asy[idx]  = resulttuple[6]
    return (Qext,Qsca)
    

# Literature references
#
#TODO (Bohren and Huffman) <insert Ref>
#    
# (Maetzler 2002): 
#       C. Maetzler, "MATLAB Functions for Mie Scattering and Absorption", 
#       Research Report 2002-08, Institut fuer Angewandte Physik, 
#       Universitaet Bern, 2002
#       downloaded from http://www.iap.unibe.ch/publications/download/201/en/
#
# (Johnson and Christy 1972):
#       P.B. Johnson and R.W. Christy, "Optical Constants of the Noble 
#       Metals", Phys. Rev. B. 1972, (6), 4370
#
# (Haiss et al. 2007):
#       W. Haiss, N.T.K. Thanh, J. Aveyard, D.G. Fernig, "Determination of
#       size and concentration of gold nanoparticles from UV-vis spectra.",
#       Anal. Chem. 2007, (79), 4215
#
# (Kreibig 1974):
#       U. Kreibig, "Electronic properties of small silver particles: the
#       optical constants and their temperature dependence",
#       J. Phys. F: Metal Phys. 1974, 4, 999
#
# (Murata and Tanaka 2010):
#       K.I. Murata and H. Tanaka, "Surface-wetting effects on the liquid-liquid
#       transition of a single-component molecular liquid.",
#       Nature Commun. 2010, 1, 16




if __name__ == "__main__":
    # execute only if run as a script
    # testing code goes here

    #import numpy as np
    import matplotlib.pyplot as plt
    
    d_nm=50.
    wavelens =	np.linspace(380,1000,500)
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





    
    





