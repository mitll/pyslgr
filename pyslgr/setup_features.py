from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

ext = Extension("LLFeatures",
                sources=["LLFeatures.pyx"],
		extra_objects = ['../slgr_engine/slgr_lib.a','../get_f0_lib/get_f0.a'],
                language="c++"
                ,include_dirs = ['../slgr_engine/',numpy.get_include()]
                ,extra_compile_args = ["-std=c++0x","-fpermissive"]
                ,extra_link_args = ['-lfftw3','-lstdc++']
                )

setup(name="LLFeatures",
      ext_modules=cythonize(ext))
