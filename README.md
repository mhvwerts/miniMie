# miniMie: calculate optical cross sections of spherical particles

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7657794.svg)](https://doi.org/10.5281/zenodo.7657794)

`miniMie` is a small scientific Python module for calculating optical cross sections (extinction, scattering, absorption) of spherical (nano)particles in a medium using Mie theory. It depends on `numpy` and `scipy`. For plotting, `matplotlib` is needed. 

The module finds its origins in the script that we used and published in:
J. R. G. Navarro and M. H. V. Werts, "Resonant light scattering spectroscopy of gold, silver and gold-silver alloy nanoparticles and optical detection in microfluidic channels", [*Analyst* **2013**, *138*, 583-592](https://doi.org/10.1039/c2an36135c).

The Mie calculation code is based on MATLAB code from the [report by C. Mätzler](https://boris.unibe.ch/146551/1/201-1.pdf). Included in `miniMie` are the dielectric functions of gold and silver from Johnson & Christy, [*Phys. Rev. B* **1972**, *6*, 4370](https://doi.org/10.1103/PhysRevB.6.4370), with the possibility of applying a mean-free path correction as described by Haiss et al, [*Anal. Chem* **2007**, *79*, 4215](https://doi.org/10.1021/ac0702084).

The extinction cross section calculated by `miniMie` can be readily converted into the molar extinction coefficient (see the example Jupyter Notebook). From the absorption and extinction cross sections, one may deduce the photothermal efficiency.

This implementation is distributed under the CeCILL license (a GNU GPL-compatible license). See: [https://cecill.info/index.en.html](https://cecill.info/index.en.html)

The basic code for doing simple spectral Mie calculations functions correctly. Furthermore, we have included the functions of the [`mie`](https://github.com/clegett/mie) library by [C. Legett](https://github.com/clegett), which is a more recent and complete Python port of the Mätzler MATLAB code. In the demo spectrum calculations, it yields results identical to the original miniMie.


## Example notebooks

* [Calculation of molar extinction coefficients of nanospheres](https://github.com/mhvwerts/miniMie/blob/master/Example%20-%20Extinction%20coefficients%20of%20gold%20nanospheres.ipynb)


## To-do list

* Include further application examples, in particular as Jupyter Notebooks
* Add benchmark calculations for testing purposes (*e.g.*, see the Mätzler report, and other sources)
* Add more dielectric functions, in particular analytic models using a limited set of parameters to accurately describe the dielectric function over a range of wavelengths.
* Use more functionalities from the `miniMie.clegett_mie` submodule.


## Links to other programs of interest

* [`Mie-Simulation-Maetzler-MATLAB-code`](https://github.com/garatbeo/Mie-Simulation-Maetzler-MATLAB-code) contains a copy of Mätzler's original MATLAB code, together with the accompanying report.
* [scattport.org](https://scattport.org/index.php/programs-menu/mie-type-codes-menu/111-mie-matlab-maetzler) contains another copy of the MATLAB code.
* [`mie`](https://github.com/clegett/mie), a more recent and more complete port of the Mätzler MATLAB code to Python, without dielectric functions for gold and silver.
* [`pyMieScatt`](https://github.com/bsumlin/PyMieScatt), a recent, comprehensive Python library for Mie calculations, with code for solving the inverse Mie problem (obtaining the complex refractive index from absorption and scattering measurements). Dielectric functions for gold and silver are not included. The `pyMieScatt` library is described in [*J. Quant. Spectrosc. Radiat. Transf.* **2018**, *205*, 127](https://doi.org/10.1016/j.jqsrt.2017.10.012).
