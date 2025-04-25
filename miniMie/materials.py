# -*- coding: utf-8 -*-
"""
Classes describing the various optical dielectric functions of materials

Part of miniMie by M. H. V. Werts

See 'core.py' for literature references
"""
import numpy as np
from numpy import pi, sqrt
from scipy import interpolate



class Material:
    """Base Material class
    
    Defines a material with a constant (complex or real) refractive index
        ncmplx
    """
    def __init__(self, ncmplx = 1.0):
        self.ncmplx = ncmplx


    def get_ncmplx_vector(self, wvln_nm):
        """generate a vector of complex dielectric function of a material
        sampled to the wavelengths (nm) in the input vector
        
        In:
            wvln_nm    vector of wavelengths (float)
                       for which dielectric function
                       should be calculated
        """
        ncmplx_wvln = np.ones_like(wvln_nm) * self.ncmplx
        return ncmplx_wvln
    


class JC_gold:
    """Gold as reported by Johnson & Christy (1972)
    
       An optional mean free-path correction can be applied as described by
       Haiss et al. (2007) by setting an effective radius of the particle
           MFPradius_nm
    """
    def __init__(self, MFPradius_nm = None):
        self.MFPradius_nm = MFPradius_nm
        
        # complex dielectric function of gold Johnson and Christy
        E = np.array([0.64,0.77,0.89,1.02,1.14,1.26,1.39,1.51,1.64,1.76,1.88,
                      2.01,2.13,2.26,2.38,2.50,
                      2.63,2.75,2.88,3,3.12,3.25,3.37,3.5,3.62,3.74,3.87,3.99,
                      4.12,4.24,4.36,4.49,4.61,4.74,4.86,
                      4.98,5.11,5.23,5.36,5.48,5.6,5.73,5.85,5.98,6.1,6.22,
                      6.35,6.47,6.6])
        n = np.array([0.92,0.56,0.43,0.35,0.27,0.22,0.17,0.16,0.14,0.13,0.14,
                      0.21,0.29,0.43,0.62,1.04,
                      1.31,1.38,1.45,1.46,1.47,1.46,1.48,1.50,1.48,1.48,1.54,
                      1.53,1.53,1.49,1.47,1.43,1.38,1.35,
                      1.33,1.33,1.32,1.32,1.30,1.31,1.30,1.30,1.30,1.30,1.33,
                      1.33,1.34,1.32,1.28])
        k = np.array([13.78,11.21,9.519,8.145,7.150,6.350,5.663,5.083,4.542,
                      4.103,3.697,3.272,
                      2.863,2.455,2.081,1.833,1.849,1.914,1.948,1.958,1.952,
                      1.933,1.895,1.866,1.871,1.883,1.898,
                      1.893,1.889,1.878,1.869,1.847,1.803,1.749,1.688,1.631,
                      1.577,1.536,1.497,1.460,1.427,1.387,
                      1.350,1.304,1.277,1.251,1.226,1.203,1.188])
        
        self.E = E
        self.n = n
        self.k = k
        
        # mean-free path parameters
        # following values are from Haiss et al.
        self.FV = 1.4E8 # Fermi velocity in cm/s - needed by mean free path correction
        self.OMP = 138. # plasma frequency in Hz/1E+14
        self.OM0 = 0.333 # collision frequency in Hz/1E+14  
        
    def ncmplx_mfpcorr(self, ncmplx_bulk, waveln):
        """Mean free path correction
        
        Not necessarily used by all materials
        
        The input takes the "bulk" complex dielectric function
        (as a vector)
        together with a vector of the wavelengths
        and material parameters
        
        Returns the MFP-corrected dielectric function
        
        
        Adapted from Haiss FORTRAN code (Haiss et al. 2007).
        This code is in cgs units, which was maintained here.
        We used the code with minimal changes in order to avoid errors;
        this is why there are UPPERCASE variable names...
        
        radius, waveln in nanometers
        FV in cm/s, OMP, OM0 in 1E-14 Hz
        """
        FV = self.FV
        OMP = self.OMP
        OM0 = self.OM0
        radius = self.MFPradius_nm
        
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

    def get_ncmplx_vector(self, wvln_nm):
        """
        Generate a vector of dielectric function values of the material

        Parameters
        ----------
        wvln_nm : np.ndarray(dtype=float)
            vector of wavelengths for which dielectric function
            should be evaluated.

        Returns
        -------
        ncmplx_wvln : np.ndarray(dtype=complex)
            vector of complex dielectric function of the material
            sampled to the wavelengths (nm) in the input vector.

        """
        ncmplx = complex(1,0) * self.n + complex(0,1) * self.k 
        ncmplx_interpol=interpolate.interp1d(self.E, ncmplx, kind='cubic')
        
        # wavelength-sampled dielectric function
        ncmplx_wvln0 = ncmplx_interpol(1240.0/(wvln_nm))
        
        # mean-free path correction if required
        if not self.MFPradius_nm == None:
            ncmplx_wvln = self.ncmplx_mfpcorr(ncmplx_wvln0, wvln_nm)   
        else:
            ncmplx_wvln = ncmplx_wvln0
        return ncmplx_wvln



class JC_silver(JC_gold):
    """Silver as reported by Johnson & Christy (1972)
    
    An optional mean free-path correction can be applied as described by
    Haiss et al. (2007) by setting an effective radius of the particle
        MFPradius_nm
    """
    def __init__(self, MFPradius_nm = None):
        self.MFPradius_nm = MFPradius_nm
        
        # complex dielectric function of silver Johnson and Christy
        E = np.array([0.64,0.77,0.89,1.02,1.14,1.26,1.39,1.51,1.64,1.76,1.88,
                      2.01,2.13,
                      2.26,2.38,2.50,2.63,2.75,2.88,3,3.12,3.25,3.37,3.5,3.62,
                      3.74,3.87,3.99,4.12,
                      4.24,4.36,4.49,4.61,4.74,4.86,4.98,5.11,5.23,5.36,5.48,
                      5.6,5.73,5.85,5.98,6.1,6.22,6.35,6.47,6.6])
        n = np.array([0.24,0.15,0.13,0.09,0.04,0.04,0.04,0.04,0.03,0.04,0.05,
                      0.06,0.05,
                      0.06,0.05,0.05,0.05,0.04,0.04,0.05,0.05,0.05,0.07,0.1,0.14,
                      0.17,0.81,1.13,1.34,
                      1.39,1.41,1.41,1.38,1.35,1.33,1.31,1.3,1.28,1.28,1.26,
                      1.25,1.22,1.20,1.18,1.15,1.14,1.12,1.10,1.07])
        k = np.array([14.08,11.85,10.10,8.828,7.795,6.692,6.312,5.727,5.242,
                      4.838,4.483,
                      4.152,3.858,3.586,3.324,3.093,2.869,2.657,2.462,2.275,
                      2.07,1.864,1.657,1.419,
                      1.142,0.829,0.392,0.616,0.964,1.161,1.264,1.331,1.372,
                      1.387,1.393,1.389,1.378,
                      1.367,1.357,1.344,1.342,1.336,1.325,1.312,1.296,1.277,
                      1.255,1.232,1.212])
        self.E = E
        self.n = n
        self.k = k
        
        # mean-free path parameters
        # following values are from Kreibig (1974) + Murata&Tanaka (2010)
        self.FV = 1.4E8 # Fermi velocity in cm/s - needed by mean free path correction
        self.OMP = 137. # plasma frequency in Hz/1E+14
        self.OM0 = 0.27 # collision frequency in Hz/1E+14
