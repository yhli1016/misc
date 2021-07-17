from setuptools import Extension, setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np


ext_modules = [
    Extension(
        "core",
        ["core.pyx"],
        include_dirs=[np.get_include()],
    )
]

setup(
    name="core",
    cmdclass={"build_ext": build_ext},
    ext_modules=cythonize(ext_modules)
)
