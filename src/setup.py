from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        'engine',  # Name of the module
        ['./c_backend/engine.cpp'],  # Source file(s)
        include_dirs=[pybind11.get_include()],  # Include pybind11 headers
        language='c++',
        extra_compile_args=['-std=c++17', '-fopenmp', '-O3'],  # Compiler flags
        extra_link_args=['-fopenmp', '-O3'],  # Linker flags
    ),
]

setup(
    name='example',
    version='0.0.1',
    ext_modules=ext_modules,
)