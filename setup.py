from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration
import distutils.sysconfig as ds

setup(name='PyFastCentroid',
      version='0.1',
      description='Fast PSF centroiding routines for Python.',
      author='Hannu Parviainen',
      author_email='hpparvi@gmail.com',
      url='',
      package_dir={'pyfc':'src'},
      packages=['pyfc'],
      ext_modules=[Extension('pyfc.gaussianf', ['src/gaussian.f90'], libraries=['gomp','m']),
                   Extension('pyfc.lorentzianf', ['src/lorentzian.f90'], libraries=['gomp','m'])]
     )
