from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("LDPC",
                 sources=["LDPC.pyx", "decode_BEC.c"],
                 include_dirs=[numpy.get_include()],extra_compile_args=['-std=c99'])],
)
