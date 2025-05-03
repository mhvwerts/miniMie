# miniMie: calculate optical cross sections of spherical particles

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7657794.svg)](https://doi.org/10.5281/zenodo.7657794)

miniMie is a small scientific Python library for calculating optical cross sections (extinction, scattering, absorption) of spherical (nano)particles in a medium using Mie theory. It depends on `numpy` and `scipy`, and on `matplotlib` for plotting. 

miniMie calculates extinction, scattering and absorption spectra via `Mie_spectrum()`. The extinction cross sections calculated by miniMie can be readily converted into the molar extinction coefficients (see the example Jupyter Notebook). From the absorption and extinction cross sections, one may deduce the photothermal efficiency. miniMie can also calculate the angular distribution of scattering intensity via `Mie_tetascan()`.

The library finds its origins in the script that we used and published in:
J. R. G. Navarro and M. H. V. Werts, "Resonant light scattering spectroscopy of gold, silver and gold-silver alloy nanoparticles and optical detection in microfluidic channels", [*Analyst* **2013**, *138*, 583-592](https://doi.org/10.1039/c2an36135c).

The Mie calculation code is based on MATLAB code from the [report by C. Mätzler](https://doi.org/10.7892/boris.146550) through the [`mie`](https://github.com/clegett/mie) library by [C. Legett](https://github.com/clegett). The [`mie`](https://github.com/clegett/mie) library has been integrated and adapted into the miniMie code base, as it gives results identical to our initial Python port. 

Also included in miniMie are the dielectric functions of gold and silver from Johnson & Christy, [*Phys. Rev. B* **1972**, *6*, 4370](https://doi.org/10.1103/PhysRevB.6.4370), with the possibility of applying a mean-free path correction as described by Haiss et al, [*Anal. Chem* **2007**, *79*, 4215](https://doi.org/10.1021/ac0702084). The dielectric function is interpolated to find the values at specific, user-defined wavelengths.

This implementation is distributed under the CeCILL license (a GNU GPL-compatible license). See: [https://cecill.info/index.en.html](https://cecill.info/index.en.html)



## Examples

### Jupyter Notebooks

* [Calculation of molar extinction coefficients of nanospheres](https://github.com/mhvwerts/miniMie/blob/master/Example%20-%20Extinction%20coefficients%20of%20gold%20nanospheres.ipynb)


### Scripts
* [Calculate and plot extinction and scattering spectra of gold nanoparticles](https://github.com/mhvwerts/miniMie/blob/master/Mie_example_with_plot.py)
* [Calculate and plot angular scattering by nanoparticles](https://github.com/mhvwerts/miniMie/blob/master/Mie_angular_scattering.py)
* [Reproduce angular scattering calculation results for gold nanoparticles by Shortell et al. (Opt. Express 2016)](https://github.com/mhvwerts/miniMie/blob/master/Mie_angular_scattering_gold.py)



## Suggestions for future work

* Include further application examples, in particular as Jupyter Notebooks
* Add benchmark calculations for testing purposes, in particular spectral properties relying on the included dielectric functions
* Add more dielectric functions
* Include the Mie calculation code from [`miepython`](https://github.com/scottprahl/miepython), which is not based on Bohren & Huffman, but on Wiscombe, and may be numerically superior in certain cases. This mainly concerns [`miepython/mie_nojit.py`](https://github.com/scottprahl/miepython/blob/main/miepython/mie_nojit.py) and [`miepython/core.py`](https://github.com/scottprahl/miepython/blob/main/miepython/core.py) with some interfacing and clean-up.



## Links

* [`mie`](https://github.com/clegett/mie), a recent port of the Mätzler MATLAB code to Python, without dielectric functions for gold and silver.
* [`Mie-Simulation-Maetzler-MATLAB-code`](https://github.com/garatbeo/Mie-Simulation-Maetzler-MATLAB-code) contains a copy of Mätzler's original MATLAB code, together with the accompanying report.
* [scattport.org](https://scattport.org/index.php/programs-menu/mie-type-codes-menu/111-mie-matlab-maetzler) contains another copy of the MATLAB code.
* [`pyMieScatt`](https://github.com/bsumlin/PyMieScatt), a comprehensive Python library for Mie calculations, with code for solving the inverse Mie problem (obtaining the complex refractive index from absorption and scattering measurements). Dielectric functions for gold and silver are not included. The `pyMieScatt` library is described in [*J. Quant. Spectrosc. Radiat. Transf.* **2018**, *205*, 127](https://doi.org/10.1016/j.jqsrt.2017.10.012).
* [`miepython`](https://github.com/scottprahl/miepython) is a pure Python module to calculate light scattering for non-absorbing, partially-absorbing, or perfectly-conducting spheres. Mie computations are done with a procedure that differs from the Mätzler/Bohren&Huffman code (but gives similar results, of course). `miepython` has [good documentation](https://miepython.readthedocs.io/en/latest/), also describing different conventions used by different codes, and benchmarking against results from literature.
