# miniMie: calculate optical cross sections of spherical particles

This repository is very much a work in progress. It will be updated irregularly. The basic code for doing simple Mie calculations is functioning correctly, but it needs further documentation and illustration.

miniMie is a small scientific Python module for calculating optical cross sections (extinction, scattering, absorption) for spherical (nano)particles in a medium using Mie theory. It depends on numpy and scipy. For plotting, also matplotlib is needed.

The module is essentially the script that we used and published in:
J. R. G. Navarro and M. H. V. Werts, "Resonant light scattering spectroscopy of gold, silver and gold-silver alloy nanoparticles and optical detection in microfluidic channels", *Analyst* **2013**, *138*, 583-592.

The calculated extinction cross section can be readily converted into the molar extinction coefficient. From the absorption and extinction cross sections, one may deduce the photothermal efficiency.

This implementation is distributed under the CeCILL license (a GNU GPL-compatible license). See: [https://cecill.info/index.en.html](https://cecill.info/index.en.html)

## TO DO

* Include further application examples, in particular as Jupyter Notebooks
* Clean up code
* Add benchmarking results for testing purposes
