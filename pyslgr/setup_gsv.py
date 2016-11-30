from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

ext = Extension("GSV",
                sources=["GSV.pyx"],
		extra_objects = ['../slgr_engine/slgr_lib.a'],
                language="c++"
                ,include_dirs = ['../slgr_engine/',numpy.get_include()]
                ,extra_compile_args = ["-std=c++0x","-fpermissive"]
                ,extra_link_args = ['-lfftw3','-lstdc++']
                )

setup(name="GSV",
      ext_modules=cythonize(ext))
