# PyLineProfile

Calculation of accurate undersampled (FWHM ~ 1 pixel) 1D Gaussian and Lorentzian
profiles. Useful in PSF fitting, etc. The profiles are calculated by integrating
them analyticallly over each pixel, and the main parts of the code are written
in Fortran for efficiency.

## Requirements
  - Fortran compiler
  - NumPy

## Installation

### From github

    git clone https://github.com/hpparvi/psf-centroid.git PyFC
    cd

    python setup.py config_fc --fcompiler=gnu95 --opt="-Ofast" --f90flags="-fopenmp -march=native" build
    python install --user
