from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration
import distutils.sysconfig as ds

setup(name='PyCntr',
      version='0.1',
      description='Fast centroiding routines for Python.',
      author='Hannu Parviainen',
      author_email='hpparvi@gmail.com',
      url='',
      package_dir={'pycntr':'src'},
      packages=['pycntr'],
      ext_modules=[Extension('pycntr.cntrf', ['src/cntr.f90','src/cntr.pyf'], libraries=['gomp','m'])]
     )
